import xarray as xr, pandas as pd
import numpy as np    
from netCDF4 import Dataset
from datetime import datetime,timedelta
from scipy import stats
from glob import glob
from scipy.io import readsav

months = [str(i+6).zfill(2) for i in range(3)]
years = [str(i+2005) for i in range(12)]
lonGrid = np.arange(230,291,2)
latGrid = np.arange(25,51,2)

def grid_and_average(lon, lat, lonGrid, latGrid, varin):
    import numpy as np
    nlatGrid  = len(latGrid)                           # the number of latitude grids
    nlonGrid  = len(lonGrid)                           # the number of longitude grids
    latStep   = latGrid[1] - latGrid[0]                # calculate the latstep for the regular latitudes
    lonStep   = lonGrid[1] - lonGrid[0]                # calculate the lonstep for the regular longitude

    #s  = varin.squeeze().shape
    varout = np.ones([np.size(lonGrid),np.size(latGrid)])*np.nan    
    for i in range(nlonGrid):
        for j in range(nlatGrid):
            lon_low = lonGrid[i]-lonStep/2.0 ; lon_high = lonGrid[i]+lonStep/2.0
            lat_low = latGrid[j]-latStep/2.0 ; lat_high = latGrid[j]+latStep/2.0
            iH = np.where(((lon >= lon_low) & (lon < lon_high)) & ((lat >= lat_low) & (lat < lat_high)))
            if len(iH[0]) != 0:
                varout[i, j] = np.nanmean(varin[iH[0]])
    return varout.T

def get_parcel_history(year,month):
    lonGrid = np.arange(230,291,2)
    latGrid = np.arange(25,51,2)
    days = [str(i+1).zfill(2) for i in range(30)]
    fs = glob('/mnt/data/ice2/wyu/ttl_traj_out/200223_s100_i100_ERAi_MLS_6hrly_back_'+year+month+\
                  '_NA/traj_s100_i100_'+year+'_????_I??????.sav')
    parcel_history1 = np.zeros(len(fs))
    init_lon = np.zeros(len(fs))
    init_lat = np.zeros(len(fs))
    for i in range(len(fs)):
        f = readsav(fs[i])
        init_lon[i] = f.lonrec[0]
        init_lat[i] = f.latrec[0]
        if len(np.where(np.any([f.latrec<25,f.latrec>50,f.lonrec<230,f.lonrec>290],axis=0))[0])==0:
            parcel_history1[i] = 480
        else:
            parcel_history1[i] = np.where(np.any([f.latrec<25,f.latrec>50,f.lonrec<230,f.lonrec>290],axis=0))[0][0]
    print (len(fs))
    parcel_history = grid_and_average(init_lon, init_lat, lonGrid, latGrid, parcel_history1)
    return parcel_history

for i,year in enumerate(years):
    for j,month in enumerate(months):
        print (year,month)
        p1 = get_parcel_history(year,month)
        fout = '/mnt/data/sn2/wyu/mls/parcel_history/parcel_history'+year+month+'_grided_center.nc'
        if os.path.exists(fout):
            os.remove(fout)
        f2 = Dataset(fout,'w')

        #Time=f2.createDimension('Time',None)
        Lat=f2.createDimension('Lat',13)
        Lon=f2.createDimension('Lon',31)

#         Time=f2.createVariable('Time','int',('Time',))
#         Time.units="day"
#         Time.long_name="Time"
#         Time[:]= int(month)+int(year)*100

        Lon=f2.createVariable('Lon','f',('Lon',))
        Lon.units="degree"
        Lon.long_name="Longitude"
        Lon[:]= lonGrid

        Lat=f2.createVariable('Lat','f',('Lat',))
        Lat.units="degree"
        Lat.long_name="Latitude"
        Lat[:]= latGrid
        
        Parcel_history = f2.createVariable('Parcel_history','f',('Lat','Lon',))
        Parcel_history.unit="hours"
        Parcel_history.long_name="parcel history"
        Parcel_history[:]=p1
        
        f2.close()
