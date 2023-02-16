def mvr_filter(y,X):
    '''
    multivariate regression filter 
    input: 
        y: target array
        X: array to filter out
    output: 
        filtered y array
    '''
    from sklearn.linear_model import LinearRegression
    import numpy as np
    reg = LinearRegression().fit(X, y)
    
    return y - np.sum(reg.coef_*X,axis=1)

def corr_nan(x,y):
    '''
    calculate correlation coefficient considering nans 
    output:
        a dictionary containing: correlation coefficient (corr), p_value, and number of invalid values (n)
    '''
    from scipy import stats
    import numpy as np
    i_finite = np.where((np.isfinite(x)) & (np.isfinite(y)))[0]
    print ('Correlation total values: {}, Valid values: {}'.format(len(x),len(i_finite)))
    slope, intercept, r_value, p_value, std_err = stats.linregress(x[i_finite],y[i_finite])
    result = {'corr':r_value,'p_value':p_value,'n':len(i_finite)}
    return result 



def linear_regress(x,y, sl=0.95,auto_correlation=True,print_info=False):
    '''
    linear regression 
    input: x, y: 1-d array 
        sl: significant level, default 0.95
        auto_correlation: consider auto_correlation or not
    output: 
        result dictionary, including slope, intercept, r_value, p_value, std_err, confidence_level
    '''
    from scipy import stats
    import numpy as np
    i_finite = np.where((np.isfinite(x)) & (np.isfinite(y)))[0]
    if print_info==True:
        print ('Linear regression total values: {}, Valid values: {}'.format(len(x),len(i_finite)))
    if len(x)>2*len(i_finite):
        print ('WARNING: half of the data is not valid')
    if len(i_finite):
        y = y[i_finite]
        x = x[i_finite]
        n = len(y)
        slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
#------calculate confidence level considering auto_correlation--------            
        if auto_correlation:
            aa,bb = slope,intercept
            residual = (y-aa*x-bb)
            xx = residual
            z = np.correlate(xx,xx,mode='full')
            z = z[len(z)//2:]
            lag1 = z[1]/z.max()
            df = (n*(1-lag1)/(1+lag1))-2
        else:
            df = n-2
        confidence_level = stats.t.ppf(1-(1-sl)/2,df)*std_err
        result = {'slope':slope,'intercept':intercept,
            'r_value':r_value,'p_value':p_value,'std_err':std_err,
            'df':df, 'confidence_level':confidence_level,
                  'sig_level':sl,'auto_correlation':auto_correlation
        }
    else:
        result = {'slope':np.nan,'intercept':np.nan,
            'r_value':np.nan,'p_value':np.nan,'std_err':np.nan,
            'df':np.nan, 'confidence_level':np.nan,
                 'sig_level':sl,'auto_correlation':auto_correlation}
    return result
    
def detrend(x):
    '''
    function to calculate a detrended time series 
    x: time series of 1-4d, with time to be the first dimension
    return: detrended x 
    '''
    import numpy as np
    s = x.shape
    trend = np.zeros(s)
    if len(s) == 1:
        trend = linear_regress(np.arange(s[0]),x)['slope']*np.arange(s[0])   
    elif len(s) == 2:
        for j in range(s[1]):
            trend[:,j] = detrend(x[:,j])
    elif len(s) == 3:
        for j in range(s[2]):
            trend[:,:,j] = detrend(x[:,:,j])
    elif len(s) == 4:
        for j in range(s[3]):
            trend[:,:,:,j] = detrend(x[:,:,:,j])
    return x-trend
        
