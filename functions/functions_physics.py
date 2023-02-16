def SWV(temp,pr): 
    '''
    input: 
        temp: temperature in K
        pr: pressure in hPa
    output:
        saturated water vapor mixing ratio in ppmv
    '''
    esi=np.zeros(np.shape(temp),dtype='f')
    i_h2o=np.where(temp>273.15)
    n_h2o=np.size(i_h2o)
    if n_h2o>0:
        esi[i_h2o]=np.exp(54.842763-6763.22/temp[i_h2o]-4.210*np.log(temp[i_h2o])\
        +0.000367*temp[i_h2o]+np.tanh(0.0415*(temp[i_h2o]-218.8))\
        *(53.878-1331.22/temp[i_h2o])-9.44523*np.log(temp[i_h2o]+0.014025*temp[i_h2o]))/100.
    i_ice=np.where((temp>0.0) & (temp<273.15))
    n_ice=np.size(i_ice)
    if n_ice>0:
        esi[i_ice]=np.exp(9.550426-5723.265/temp[i_ice]+3.53068*np.log(temp[i_ice])\
        -0.00728332*temp[i_ice])/100.0
    
    esi[esi>=pr] = np.nan #for ext
    h2o=esi/(pr-esi)
    
    return h2o*1.0e6



def number_density(lev,T,mr):
    '''
    convert mixing ratio to number density 
    input: 
         lev: pressure in hPa
         T: temperature in K
         mr: mixing ratio in mol/mol
    output: 
         nd: neutral density in molecular/cm-3
    '''

    K = 1.3806504e-23
    n = lev*100/K/T
    nd = n*mr*1e-6
    return nd 

def get_cooling_power(ffc):
    '''
    calculate cooling power
    input: 
        ffc: xarray file containing T and QNO
    oputut:
        cooling power (watt/m^3)
    '''
    R = 287.05
    cp = 1.006*1000
    rho = ffc.lev*100/R/ffc.T
    cooling_power = ffc.QNO*cp*rho
    return -cooling_power


def power_spectrum(y):
    '''
    calculate the power spectrum 
    input:
        y: 1d array 
    output:
        power spectrum of y in different period
    '''
    x = np.arange(len(y))
    yi = np.arange(len(y))

    mask = np.isfinite(y)
    yfiltered = np.interp(yi, yi[mask], y[mask])

    ft = np.fft.fft(yfiltered)
    xf = np.fft.fftfreq(len(x),1/len(x))
    xp = np.where(xf>=0)[0]
    f_ps = np.sqrt(ft.real[xp]**2+ft.imag[xp]**2)

    xp_sort=np.argsort(xf[xp])
    return f_ps[xp_sort]

def cold_point_tropopause(t,lev,method='CubicSpline'):
    '''
    calculate the cold point tropopause
    input:
        t: temperature, 1d (lev), 2d(time-lev), or 3d(time-lev-lat), or 4d(time-lev-lat-lon)
        lev: pressure levels
        method: CubicSpline or interp1d
    output: 
        tropopause pressure 
        tropopause temperature 
        tropopause_relative_coordinates: the same size as t, relative coordinate of the tropopause pressure, positive --> over the tropopause
    '''
    def cpt_recursive(t,lev,method):
        from scipy.interpolate import CubicSpline,interp1d
        lev_fine = np.arange(50,110)
        s = t.shape
        methods = {'CubicSpline':CubicSpline,'interp1d':interp1d}
        if len(s) == 1:
            f = methods[method](lev,t)
            t_fine = f(lev_fine)
            i_tropopause = np.argmin(t_fine)
            tropopause_p = lev_fine[i_tropopause]
            tropopause_t = t_fine[i_tropopause]
            tropopause_relative_coordinates = tropopause_p-lev
        
        else:
            tropopause_p = np.zeros([s[0],*s[2:]])
            tropopause_t = np.zeros([s[0],*s[2:]])
            tropopause_relative_coordinates = np.zeros(s)
            if len(s) == 2:
                for i in range(s[0]):
                    tropopause_p[i],tropopause_t[i],tropopause_relative_coordinates[i] = cpt_recursive(t[i],lev,method=method)
            elif len(s) == 3:
                for j in range(s[2]):
                    tropopause_p[:,j],tropopause_t[:,j],tropopause_relative_coordinates[:,:,j] = cpt_recursive(t[:,:,j],lev,method=method)
            else:
                for k in range(s[3]):
                    tropopause_p[:,:,k],tropopause_t[:,:,k],tropopause_relative_coordinates[:,:,:,k] = cpt_recursive(t[:,:,:,k],lev,method=method)
        return tropopause_p,tropopause_t,tropopause_relative_coordinates
    
    tropopause_p,tropopause_t,tropopause_relative_coordinates = cpt_recursive(t,lev,method=method)
    
    result = {
        'tropopause_p':tropopause_p,
        'tropopause_t':tropopause_t,
        'tropopause_relative_coordinates':tropopause_relative_coordinates
    }
    return result

def get_dtdy_term(t,vstar,lat):
    '''
    function to calculate the dtdy term in thermal dynamic equation  
    '''
    ae = 6370000 
    return np.gradient(t,(lat/180*np.pi),axis=2)*vstar/ae