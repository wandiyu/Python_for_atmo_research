def read_nc_data_universal(filename,varlist,input_dimension, output_dimension,**kw):
    '''
    read data 
    input: 
        filename: filename 
        varlist: a list of variable names 
        input_dimension: input dimension, 1,2,3, or 4
        output_dimension: output dimension, 1,2,3, or 4
        **kw:
            lat_range_n: range of averaging latitude
            lat_range_s: range of averaging latitude,lat_range_s<lat_range_n
            lev_high: higher bound of pressure
            lev_low: lower bound of pressure 
                    if output_dimension == 1, use lev_high == lev_low to specify a level for 1-d output
            silent: do not print out file information
            dim_names: dimension names, format: {'Time':'Time','Lev':'Lev','Lat':'Lat','Lon':'Lon'}
    '''
    import numpy as np
    from datetime import datetime
    from functions_basics import geo_avg,find_index
    from functions_plot import print_lat
    from netCDF4 import Dataset
    if 'silent' not in kw:
        print ('now reading file: {}\n variables: {} \n time now: {}'.format(filename, ','.join(varlist), datetime.now()))
    if output_dimension>input_dimension:
        raise ValueError('output dimension should not be larger than input dimension') 
    file = Dataset(filename)
    var_output = {}
    dim_names = {'Time':'Time','Lev':'Lev','Lat':'Lat','Lon':'Lon'}
    if 'dim_names' in kw:
        dim_names = {dim:kw['dim_names'][dim] for dim in dim_names}

    for var_name in varlist:
        var = file.variables[var_name][:].data
        if 'lev_high' in kw:
            lev = file.variables[dim_names['Lev']][:].data
            i_lev_low = find_index(lev,kw['lev_low'],method='lower')
            i_lev_high = find_index(lev,kw['lev_high'],method='higher')
            unit = file.variables['Lev'].unit
            if 'silent' not in kw:
                print (' Choosing levels between {:.3f} {} and {:.3f} {}'.format(lev[i_lev_low], unit,lev[i_lev_high],unit))
            if i_lev_low>i_lev_high:
                i_lev_high, i_lev_low = i_lev_low, i_lev_high
            lev_out = lev[i_lev_low:i_lev_high+1]
            var = var[:,i_lev_low:i_lev_high+1]
        if output_dimension == input_dimension:
            var_output[var_name] = var
        else:
            if output_dimension == 3: # output_dimension == 3, input_dimension = 4
                var_output[var_name] = np.nanmean(var,axis=-1)
            else:
                if input_dimension>=3: # output_dimension == 1 or 2, input_dimension = 3 or 4 
                    lat = file.variables[dim_names['Lat']][:].data
                    if 'lat_range_n' in kw:
                        i_lat = np.where((lat>=kw['lat_range_s']) & (lat<=kw['lat_range_n']))[0]
                        if 'silent' not in kw:
                            print (' Averaging between {} and {}'.format(print_lat(kw['lat_range_s']),print_lat(kw['lat_range_n'])))
                        lat = lat[i_lat]
                        var = var[:,:,i_lat]
                    var = geo_avg(var,lat,dim = 2)
                if output_dimension == 2:
                    var_output[var_name] = var
                elif output_dimension == 1:
                    if kw['lev_high']==kw['lev_low']:
                        var_output[var_name] = var[:,0]
                    else:
                        raise SyntaxError('Please use lev_high == lev_low to specify a level for 1-d output') 
    axes = ['Time','Lev','Lat','Lon'][:output_dimension]
    for axis in axes:
        var_output[axis] = file.variables[dim_names[axis]][:].data
    if 'lev_high' in kw:
        var_output['Lev'] = lev_out
    return var_output

def add_dimension_to_filename(f):
    '''
    read and add a dimension as a suffix as file name
    '''
    import os
    from netCDF4 import Dataset
    i_fname_end = f.find('.nc')
    os.rename(f,f[:i_fname_end]+'_{}d'.format(len(Dataset(f).dimensions))+f[i_fname_end:])
    
    
def save_to_2d(fout,variables,mons,lev):
    '''
    save variables to 2-d .nc file (time-lev)
    input:
        fout: output file name
        variables: a dictionary of data, format: {varname:[unit,longname,values]}
        mons: time dimension
        lev: level dimension
    '''
    from netCDF4 import Dataset
    f2 = Dataset(fout,'w')

    Time=f2.createDimension('Time',None)
    Lev=f2.createDimension('Lev',len(lev))

    Time=f2.createVariable('Time','f',('Time',),zlib=True)
    Time.units="month"
    Time.long_name="Time"

    Lev=f2.createVariable('Lev','f',('Lev',),zlib=True)
    Lev.unit="hPa"
    Lev.long_name="Level"
    Lev.min=lev.min()
    Lev.max=lev.max()

    Time[:]=mons
    Lev[:]=lev

    for var in variables:
        var1=f2.createVariable(var,'f',('Time','Lev',),zlib=True)
        var1.unit=variables[var][0]
        var1.long_name=variables[var][1]
        var1[:]=variables[var][2]
    f2.close()
    
    
def save_to_3d(fout,variables,mons,lev,lat):
    '''
    save variables to 3-d .nc file (time-lev-lat)
    input:
        fout: output file name
        variables: a dictionary of data, format: {varname:[unit,longname,values]}
        mons: time dimension
        lev: level dimension
        lat: latitude dimension
    '''
    from netCDF4 import Dataset
    f2 = Dataset(fout,'w')

    Time=f2.createDimension('Time',None)
    Lev=f2.createDimension('Lev',len(lev))
    Lat=f2.createDimension('Lat',len(lat))

    Time=f2.createVariable('Time','f',('Time',),zlib=True)
    Time.units="month"
    Time.long_name="Time"
    
    Lat=f2.createVariable('Lat','f',('Lat',),zlib=True)
    Lat.units="degree"
    Lat.long_name="Latitude"
    Lat.min=lat.min()
    Lat.max=lat.max()

    Lev=f2.createVariable('Lev','f',('Lev',),zlib=True)
    Lev.unit="hPa"
    Lev.long_name="Level"
    Lev.min=lev.min()
    Lev.max=lev.max()

    Time[:]=mons
    Lev[:]=lev
    Lat[:]=lat

    for var in variables:
        var1=f2.createVariable(var,'f',('Time','Lev','Lat',),zlib=True)
        var1.unit=variables[var][0]
        var1.long_name=variables[var][1]
        var1[:]=variables[var][2]
    f2.close()
    
def save_to_4d(fout,variables,mons,lev,lat,lon):
    '''
    save variables to 4-d .nc file (time-lev-lat-lon)
    input:
        fout: output file name
        variables: a dictionary of data, format: {varname:[unit,longname,values]}
        mons: time dimension
        lev: level dimension
        lat: latitude dimension
        lon: longitude dimension 
    '''
    from netCDF4 import Dataset
    f2 = Dataset(fout,'w')
    
    Time=f2.createDimension('Time',None)
    Lat=f2.createDimension('Lat',len(lat))
    Lon=f2.createDimension('Lon',len(lon))
    Lev=f2.createDimension('Lev',len(lev))
    
    Time=f2.createVariable('Time','f',('Time',),zlib=True)
    Time.units="month"
    Time.long_name="Time"

    Lon=f2.createVariable('Lon','f',('Lon',),zlib=True)
    Lon.units="degree"
    Lon.long_name="Longitude"
    Lon.min=lon.min()
    Lon.max=lon.max()

    Lat=f2.createVariable('Lat','f',('Lat',),zlib=True)
    Lat.units="degree"
    Lat.long_name="Latitude"
    Lat.min=lat.min()
    Lat.max=lat.max()

    Lev=f2.createVariable('Lev','f',('Lev',),zlib=True)
    Lev.unit="hPa"
    Lev.long_name="Level"
    Lev.min=lev.min()
    Lev.max=lev.max()
    
    Time[:]=mons
    Lon[:]=lon
    Lat[:]=lat
    Lev[:]=lev

    for var in variables:
        var1=f2.createVariable(var,'f',('Time','Lev','Lat','Lon',),zlib=True)
        var1.unit=variables[var][0]
        var1.long_name=variables[var][1]
        var1[:]=variables[var][2]
    f2.close()
    
def read_mls_h2o(FILE_NAME):
    import h5py
    import numpy as np
    with h5py.File(FILE_NAME, mode='r') as f:

        name = '/HDFEOS/SWATHS/H2O/Data Fields/L2gpValue'
        pname = '/HDFEOS/SWATHS/H2O/Geolocation Fields/Pressure'
        subset = 0
        data = f[name][:]
        pres = f[pname][:]
        punits = f[pname].attrs['Units'].decode() 
        units = f[name].attrs['Units'].decode() 
        longname = f[name].attrs['Title'].decode()

        # Get the geolocation data
        latitude = f['/HDFEOS/SWATHS/H2O/Geolocation Fields/Latitude'][:]
        longitude = f['/HDFEOS/SWATHS/H2O/Geolocation Fields/Longitude'][:]
    
        # get the uncertainty 
        precision = f['/HDFEOS/SWATHS/H2O/Data Fields/L2gpPrecision'][:]
        status = f['/HDFEOS/SWATHS/H2O/Data Fields/Status'][:]
        quality = f['/HDFEOS/SWATHS/H2O/Data Fields/Quality'][:]
        convergence = f['/HDFEOS/SWATHS/H2O/Data Fields/Convergence'][:]

        # screening 
        data[:,:37][data[:,:37]<0.101/1.5e6] = np.nan
        data[precision<=0] = np.nan
        data = data[(status % 2 == 0) & ((status<16) | (status >32)) & (quality > 0.7) & (convergence<2)] 
        latitude = latitude[(status % 2 == 0) & ((status<16) | (status >32)) & (quality > 0.7) & (convergence<2)] 
        longitude = longitude[(status % 2 == 0) & ((status<16) | (status >32)) & (quality > 0.7) & (convergence<2)]
        
    return data,pres,latitude,longitude

def read_mls_t(FILE_NAME):
    import h5py
    with h5py.File(FILE_NAME, mode='r') as f:

        name = '/HDFEOS/SWATHS/Temperature/Data Fields/L2gpValue'
        pname = '/HDFEOS/SWATHS/Temperature/Geolocation Fields/Pressure'
        subset = 0
        data = f[name][:]
        pres = f[pname][:]
        punits = f[pname].attrs['Units'].decode() 
        units = f[name].attrs['Units'].decode() 
        longname = f[name].attrs['Title'].decode()

        # Get the geolocation data
        latitude = f['/HDFEOS/SWATHS/Temperature/Geolocation Fields/Latitude'][:]
        longitude = f['/HDFEOS/SWATHS/Temperature/Geolocation Fields/Longitude'][:]
    return data,pres,latitude,longitude

def read_obs_index(start_year,end_year,lag=0):
    '''
    read indexes to be filtered out for mvr_trend_model in functions_models.py
    indexes include: 
        ENSO3.4: https://www.ncdc.noaa.gov/teleconnections/enso/indicators/sst/
        F10.7: https://lasp.colorado.edu/lisird/data/noaa_radio_flux/
        QBO1&2, Singapore winds 30 hPa and 50 hPa: https://psl.noaa.gov/data/correlation/qbo.data, 
                    https://www.cpc.ncep.noaa.gov/data/indices/qbo.u50.index
        starting from ~1948
    input:
        start_year, end_year: the year to start and end the data
        lag: lag month for ENSO
    output:
        index ready for mvr_trend_model
    '''
    import numpy as np
    from datetime import datetime, timedelta
    with open('/glade/work/yuwandi/index/enso34.txt') as f1:
        t = f1.readlines()

    year_enso = []
    enso = []
    for i in range(0,len(t)):
        tt = t[i].split()
        year_enso.append(int(tt[0]))
        enso+=list(map(float,tt[1:]))  
        
    def julian_to_date(x):
        return datetime(1,1,1)+timedelta(days=x-1721425.5)
    with open('/glade/work/yuwandi/index/penticton_radio_flux.csv') as f:
        t = f.readlines()
    header = t[0]
    julian_date = []
    flux1 = []
    flux2 = []
    for i in range(1,len(t)):
        tt = list(map(float,t[i].strip().split(',')))
        julian_date.append(tt[0])
        flux1.append(tt[1])
        flux2.append(tt[2])
    date = list(map(julian_to_date,julian_date))
    year = np.array([d.year for d in date])
    month = np.array([d.month for d in date])
    flux1 = np.array(flux1)
    flux2 = np.array(flux2)
    solar1=[]
    solar2=[]
    year_solar=list(range(1947,2021))
    for y in year_solar:
        for m in range(1,13):
            i_month = np.where((year==y) & (month==m))[0]
            solar1.append(np.nanmean(flux1[i_month]))
            solar2.append(np.nanmean(flux2[i_month]))
    
#    year_qbo = []
#    qbo = []
#     with open('/glade/work/yuwandi/index/qbo.ucar.txt') as f:
#         t = f.readlines()
#     for i in range(len(t)):
#         tt = t[i].split()
#         year_qbo.append(int(tt[0]))
#         qbo += list(map(float,tt[1:]))
    year_qbo30 = []
    qbo30 = []
    with open('/glade/work/yuwandi/index/qbo.u30.index') as f:
        t = f.readlines()
    for i in range(4,46):
        tt = t[i].split()
        year_qbo30.append(int(tt[0]))
        qbo30 += list(map(float,tt[1:]))
        
    year_qbo50 = []
    qbo50 = []
    with open('/glade/work/yuwandi/index/qbo.u50.index') as f:
        t = f.readlines()
    for i in range(4,46):
        tt = t[i].split()
        year_qbo50.append(int(tt[0]))
        qbo50 += list(map(float,tt[1:]))
        
    years = [year_enso,year_solar,year_qbo30,year_qbo50]
    series = [enso,solar1,qbo30,qbo50]

# slice the regression indexes from start_year to end_year

    x = []
    for i,yr,sr in zip(range(4),years,series):
        i_start = len([y for y in yr if y<start_year])
        i_end = len([y for y in yr if y<=end_year])
        if i == 0:
            x.append(sr[i_start*12-lag:i_end*12-lag])
        else:
            x.append(sr[i_start*12:i_end*12])
    return np.array(x).T


def read_cmip6_index(scenario,run_num):
    '''
    read indexes to be filtered out for mvr_trend_model in functions_models.py
    read indexes according to SSP future scenario and run number (1-5)
    if run_num == 0, then read index from observations
    output: [time,n] array for ENSO,QBO1&2, and f10.7
    '''
    import numpy as np
    if run_num == 0:
        print (1)
    else:
        with open('/glade/u/home/yuwandi/h2o_21st_century/indexes/f107.txt') as f107:
            x = f107.readlines()
            f107 = [float(i) for i in x[0].split()]
        with open('/glade/u/home/yuwandi/h2o_21st_century/indexes/QBO30_{}-WACCM.00{}.txt'.format(scenario, str(run_num))) as fqbo:
            x = fqbo.readlines()
            qbo1 = [float(i) for i in x[0].split()]
        with open('/glade/u/home/yuwandi/h2o_21st_century/indexes/QBO50_{}-WACCM.00{}.txt'.format(scenario, str(run_num))) as fqbo:
            x = fqbo.readlines()
            qbo2 = [float(i) for i in x[0].split()]
        with open('/glade/u/home/yuwandi/h2o_21st_century/indexes/ENSO_{}-WACCM.00{}.txt'.format(scenario, str(run_num))) as fenso:
            x = fenso.readlines()
            enso = [float(i) for i in x[0].split()]
        X = np.zeros([4,len(enso)])
    X[0] = enso
    X[1] = qbo1
    X[2] = f107
    X[3] = qbo2
    return X.T
