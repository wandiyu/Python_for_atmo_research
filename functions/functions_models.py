def mvr_trend_model(var,index,**kw):
    '''
    model to calculate trend using multivariate regression
    input:
        var: monthly mean data, 1d (time), 2d (time,level) or 3d (time,level,latitude) array
            if input an array without a lev, must first reshape the array to have a level demension
        index: indexes that needed to be regressed out, format: 2-d [time, n_indexes],
            if index==False, then do not regress out anything
        **kw: 
            var_ref: reference variable field, must be the same dimension as var, default var
            select_month: call select month function, input target_mon
            sliding_window: calculate sliding window mean, and kw is window length (year)
            silent: do not print model information 
    output: 
        result: a dictionary containing trend, relative trend, confidence_level,
            relative_confidence_level,p_value, in */decade 
    '''
    import numpy as np
    from datetime import datetime
    from functions_basics import cal_anomaly,select_month
    from functions_statistics import linear_regress,mvr_filter
    if 'silent' not in kw:
        print ('Starting the MVR model, time now: {}'.format(datetime.now()))

#-----initiate values----------------------    
    s = var.shape
    if len(s)==1:
        var = var.reshape([s[0],1,1])
    elif len(s)==2:
        var = var.reshape([s[0],s[1],1])
    s = var.shape
    n_time,n_lev,n_lat = s
    if len(index) == 1:
        index = np.ones([s[0],1])
    if 'sliding_window' in kw:
        length_window = kw['sliding_window']
        if 'silent' not in kw:
            print ('sliding window mode: window length: {}'.format(length_window))
        n_year = n_time//12
        n_window = n_year-length_window+1
    else:
        n_window = 1
    if ('silent' not in kw) & ('select_month' in kw):
        print ('select month, keyword: {}'.format(kw['select_month']))
    trend, relative_trend, confidence_level,\
            relative_confidence_level,p_value = [np.zeros([n_window,n_lev,n_lat]) for i in range(5)]
#----------check var_ref--------------------------
    if 'var_ref' in kw:
        var_ref = np.nanmean(kw['var_ref'],axis=0)
    else:
        var_ref = np.nanmean(var,axis=0)
#---------------calculation---------------------
    for i_window in range(n_window):
        if 'sliding_window' in kw:
            window_start = i_window*12
            window_end = (i_window+length_window)*12
            window_range = range(window_start,window_end)
        else:
            window_range = range(n_time)
        for i_lev in range(n_lev):
            for i_lat in range(n_lat):
                y = cal_anomaly(var[window_range,i_lev,i_lat])  # remove seasonal cycle 
                x = np.arange(len(y)) # month numbers 
                X = index[window_range]
                if 'select_month' in kw:
                    x,y,X = [select_month(xx,kw['select_month']) for xx in [x,y,X]]
                    #here month number x is also selected, so this is the 'monthly trend'
                i_finite = np.where((np.isfinite(y)) & (np.isfinite(X.sum(axis=-1))))[0]
                if len(i_finite)>2:
                    y = y[i_finite]  # target array 
                    X = X[i_finite] # array to be filtered out 
                    y_filtered = mvr_filter(y,X) # targeted array, with indexes filtered out  
                    x = x[i_finite] # month numbers 
                    regress_result = linear_regress(x,y_filtered)
#--------------wrap up values----------------------------------                
                    trend[i_window,i_lev,i_lat] = regress_result['slope']
                    relative_trend[i_window,i_lev,i_lat] = regress_result['slope']/var_ref[i_lev,i_lat]
                    confidence_level[i_window,i_lev,i_lat] = regress_result['confidence_level']
                    relative_confidence_level[i_window,i_lev,i_lat]= \
                             regress_result['confidence_level']/var_ref[i_lev,i_lat]
                    p_value[i_window,i_lev,i_lat] = regress_result['p_value']
                else:
                    trend[i_window,i_lev,i_lat],relative_trend[i_window,i_lev,i_lat],confidence_level[i_window,i_lev,i_lat],\
                    relative_confidence_level[i_window,i_lev,i_lat]\
                    ,p_value[i_window,i_lev,i_lat] = [np.nan for i in range(5)]
                
    result = {
        'trend':trend.squeeze()*120,
        'relative_trend':relative_trend.squeeze()*120,
        'confidence_level':confidence_level.squeeze()*120,
        'relative_confidence_level':relative_confidence_level.squeeze()*120,
        'p_value':p_value.squeeze()
    }
    return result           



def pmc_0d_pressure_model(p,t,h2o,vpop=1):
    '''
    IDL code from Mark Hervig, GATS Inc., and converted to python 
;  input:
;     z.......altitude (km), fltarr(nz)
;     t.......temp (K), fltarr(nz)
;     p.......pressure (mb), fltarr(nz)
;     h2o.....h2o mixing ratio (ppmv), fltarr(nz)
    output:
;     v_ice......ice volume density (um3 / cm3), fltarr(nz)
;     m_ice......ice mass density (ng / m3), fltarr(nz)
;     h2o_ice....gas phase equivalent H2O in ice (ppmv), fltarr(nz)
    '''
    import numpy as np
    nt,nz,nlat,nlon = h2o.shape
    p = p.reshape([1,nz,1,1])

    Mww = 18.0        # molec wt. h2o, g/mol
    di  = 0.93        # density of ice, g/cm3
    R   = 8.314       # J/mol/K
    Sti = 0.12        # surface tension of ice in the presence of air, J/m2

    if vpop == 1:
        p_ice = 0.01*np.exp(9.550426-5723.265/t+3.53068*np.log(t)-0.00728332*t)    
    if vpop == 2: 
        p_ice = 0.01*10**(14.88-3059.0/t)
    if vpop == 3: 
        p_ice = 0.01*np.exp(28.868 - 6132.935 / t) 

    h2o_sat  = 1e6 *p_ice / p  # saturation mixing ratio, ppmv
    s_ice = h2o / h2o_sat      # saturation ratio

    #-   Equilibrium ice properties

    q_xs = h2o - h2o_sat    # excess h2o mix ratio, ppmv 
    q_xs[q_xs<0] = 0


    n_xs = p*1e2 *q_xs*1e-6/(R*t)  # excess mols h2o per m3 air
    v_ice = 1e6 * n_xs * Mww / di   # ice volume density, microns3 / cm3     
    m_ice = 1e9 * n_xs * Mww        # ice mass density, ng/m3       
    h2o_ice = q_xs                    # H2O(ice), ppmv
    
    result = {'h2o_ice':h2o_ice,'m_ice':m_ice,'v_ice':v_ice
        
    }
    return result

