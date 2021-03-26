import xarray as xr, pandas as pd
import numpy as np    
from netCDF4 import Dataset
from glob import glob
from datetime import datetime,timedelta

startyear=2006
startmonth=8
year=startyear
month=startmonth

def read_he5(fin,var):
    import h5py
    import numpy as np
    f1 = h5py.File(fin)
    loc = f1.require_group('/HDFEOS/SWATHS/'+var+'/Geolocation Fields')
    lon = loc['Longitude'][:]
    lat = loc['Latitude'][:]
    pre = loc['Pressure'][:]
    t   = loc['Time'][:]

    var1  = f1.require_group('/HDFEOS/SWATHS/'+var+'/Data Fields')
    value = var1['L2gpValue'][:]
    miss  = var1['L2gpValue'].attrs['MissingValue']
    prec  = var1['L2gpPrecision'][:]
    qual  = var1['Quality'][:]
    conv  = var1['Convergence'][:]
    stat  = var1['Status'][:]

    ind = np.where(np.isnan(conv) | np.isnan(qual) | np.isnan(stat))
    conv[ind[0]] = var1['Convergence'].attrs['MissingValue']
    qual[ind[0]] = var1['Quality'].attrs['MissingValue']
    stat[ind[0]] = var1['Status'].attrs['MissingValue']

    return dict(lon=lon,lat=lat,pr=pre,t=t,value=value,miss=miss,prec=prec,qual=qual,conv=conv,stat=stat)

years = [2005+i for i in range(12)]
months = [6+i for i in range(3)]
for yearin in years:
    for monthin in months:
        
        #f1 = glob('/co2/hao/download_mls/raw_data/h2o/MLS-Aura_L2GP-H2O_v04_'+str(year).zfill(4)+'d???.he5')
        f1 = glob('/mnt/data/ice1/wyu/download_mls/MLS-Aura_L2GP-O3_v04_'+str(yearin).zfill(4)+'d???.he5')
        f1.sort()
#         f2 = glob('/co2/hao/download_mls/raw_data/T/MLS-Aura_L2GP-Temperature_v04_'+str(year).zfill(4)+'d???.he5')
#         f2.sort()
        month_start = (datetime(yearin,monthin,1)-datetime(yearin,1,1)).days+1
        month_end = (datetime(yearin,monthin+1,1)-timedelta(days=1)-datetime(yearin,1,1)).days+1

        f1start = glob('/mnt/data/ice1/wyu/download_mls/MLS-Aura_L2GP-O3_v04_'+str(yearin).zfill(4)+'d'+str(month_start).zfill(3)+'.he5')
        f1end = glob('/mnt/data/ice1/wyu/download_mls/MLS-Aura_L2GP-O3_v04_'+str(yearin).zfill(4)+'d'+str(month_end).zfill(3)+'.he5')
        index_start = f1.index(f1start[0])
        index_end = f1.index(f1end[0])
        f1 = f1[index_start:index_end+1]

# f2start = glob('/co2/hao/download_mls/raw_data/T/MLS-Aura_L2GP-Temperature_v04_'+str(year).zfill(4)+'d'+str(month_start).zfill(3)+'.he5')
# f2end = glob('/co2/hao/download_mls/raw_data/T/MLS-Aura_L2GP-Temperature_v04_'+str(year).zfill(4)+'d'+str(month_end).zfill(3)+'.he5')
# index_start = f2.index(f2start[0])
# index_end = f2.index(f2end[0])
# f2 = f2[index_start:index_end+1]


# readdata 

        mlsh2o = np.zeros([1,55])
        mlstemp = np.zeros([1,55])
        lat = np.zeros(1)
        lon = np.zeros(1)
        t = np.zeros(1)

        for i in range(len(f1)):
            #mlsdata = read_he5(f1[i],'H2O')
            mlsdata = read_he5(f1[i],'O3')
            mlsh2o = np.append(mlsh2o,np.array(mlsdata['value']),axis=0)
            lat = np.append(lat,np.array(mlsdata['lat']))
            lon = np.append(lon,np.array(mlsdata['lon']))
            t = np.append(t,np.array(mlsdata['t']))
            #mlsdata = read_he5(f2[i],'Temperature')
            #mlstemp = np.append(mlstemp,np.array(mlsdata['value']),axis=0)

        mlsh2o = mlsh2o[1:]
        #mlstemp = mlstemp[1:]
        lat = lat[1:]
        lon = lon[1:]
        t = t[1:]
        level = np.array(mlsdata['pr'])

        # process time 
        mlshour = np.zeros(len(t))
        mlsyear = np.zeros(len(t))
        mlsmonth = np.zeros(len(t))
        mlsday = np.zeros(len(t))
        for i in range(len(t)):
            mlstime = (datetime(1993, 1, 1)+timedelta(seconds=t[i]))
            mlsyear[i] = mlstime.year
            mlsmonth[i] = mlstime.month
            mlsday[i] = mlstime.day
            mlshour[i] = mlstime.hour
        # convert -999 to np.nan 
        lon[np.where(lon<-800)[0]] = np.nan
        lat[np.where(lat<-800)[0]] = np.nan
        mlsh2o[np.where(mlsh2o<-800)[0]] = np.nan
        lon[np.where(lon<0)[0]] = lon[np.where(lon<0)[0]] +360
        
        #mlstheta = mlstemp*(1000/level[np.newaxis,:])**(2/7.)
    '''
        # save the data 
        outfile='/sn2/wyu/mls/satellite/MLS_L2_O3_'+str(yearin).zfill(4)+str(monthin).zfill(2)+'_global.nc'
        f=Dataset(outfile,'w',format='NETCDF4')

        index=f.createDimension('index',len(t))
        altitude=f.createDimension('altitude',len(level))

        index=f.createVariable('index','f',('index',))
        altitude=f.createVariable('altitude','f',('altitude',))

        #h2o=f.createVariable('h2o','f',('index','altitude',))
        o3 =f.createVariable('o3','f',('index','altitude',))
        #temperature=f.createVariable('temperature','f',('index','altitude',))
        #theta=f.createVariable('theta','f',('index','altitude',))
        year=f.createVariable('year','f',('index',))
        month=f.createVariable('month','f',('index',))
        day=f.createVariable('day','f',('index',))
        hour=f.createVariable('hour','f',('index',))
        longitude=f.createVariable('longitude','f',('index',))
        latitude=f.createVariable('latitude','f',('index',))


        index[:]=np.arange(len(t))
        altitude[:]=level

        year[:]=mlsyear
        month[:]=mlsmonth
        day[:]=mlsday
        hour[:]=mlshour
        longitude[:]=lon
        latitude[:]=lat
        #h2o[:]=mlsh2o*1e6
        o3[:]=mlsh2o*1e6
        #temperature[:] = mlstemp
        #theta[:] = mlstheta

        f.close()
        '''
        outfile='/sn2/wyu/mls/satellite/MLS_L2_O3_'+str(yearin).zfill(4)+str(monthin).zfill(2)+'_global.nc'
        fMLS = xr.open_dataset(outfile)
        i_na = np.where(np.all([fMLS.latitude<=49,fMLS.latitude>=10,fMLS.longitude>=230,fMLS.longitude<=290],axis=0))[0]

        # save NA data 
        outfile='/sn2/wyu/mls/satellite/MLS_L2_O3_'+str(yearin).zfill(4)+str(monthin).zfill(2)+'_NA.nc'
        #os.remove(outfile)
        f=Dataset(outfile,'w',format='NETCDF4')

        index=f.createDimension('index',len(t[i_na]))
        altitude=f.createDimension('altitude',len(level))

        index=f.createVariable('index','f',('index',))
        altitude=f.createVariable('altitude','f',('altitude',))

        o3 =f.createVariable('o3','f',('index','altitude',))
        # h2o=f.createVariable('h2o','f',('index','altitude',))
        # temperature=f.createVariable('temperature','f',('index','altitude',))
        # theta=f.createVariable('theta','f',('index','altitude',))
        z=f.createVariable('z','f',('index','altitude',))
        year=f.createVariable('year','f',('index',))
        month=f.createVariable('month','f',('index',))
        day=f.createVariable('day','f',('index',))
        hour=f.createVariable('hour','f',('index',))
        longitude=f.createVariable('longitude','f',('index',))
        latitude=f.createVariable('latitude','f',('index',))


        index[:]=fMLS.index[i_na]
        altitude[:]=level

        year[:]=mlsyear[i_na]
        month[:]=mlsmonth[i_na]
        day[:]=mlsday[i_na]
        hour[:]=mlshour[i_na]
        longitude[:]=lon[i_na]
        latitude[:]=lat[i_na]
        o3[:]=mlsh2o[i_na]*1e6
        # h2o[:]=mlsh2o[i_na]*1e6
        # temperature[:]=mlstemp[i_na]
        # theta[:]=mlstheta[i_na]

        f.close()
