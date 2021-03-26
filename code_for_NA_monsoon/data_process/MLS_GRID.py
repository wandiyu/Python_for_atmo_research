import xarray as xr, numpy as np, pandas as pd
from netCDF4 import Dataset
from glob import glob
from datetime import datetime,timedelta
from scipy.io import readsav

year = 2016
loninterval = 8
latinterval = 4
longrid = np.arange(4,364,loninterval)
latgrid = np.arange(-88,92,latinterval)

itime = datetime(year,1,1,0)
endtime = datetime(year+1,1,1,0)
days = (endtime-itime).days

# choose the level: 146.78 hPa to 68.13 hPa
lev_start =int(10)
lev_end = int(15)
f = glob('/co2/hao/download_mls/nc_file/h2o/MLS-Aura_L2GP-H2O_v04_20110101.nc')
lev = Dataset(f[0]).variables['Pressure'][lev_start:lev_end]


def GRID_ACE_VARS(lon, lat, lonGrid, latGrid, varin):
    import numpy as np
    nlatGrid  = len(latGrid)                           # the number of latitude grids
    nlonGrid  = len(lonGrid)                           # the number of longitude grids
    latStep   = latGrid[1] - latGrid[0]                # calculate the latstep for the regular latitudes
    lonStep   = lonGrid[1] - lonGrid[0]                # calculate the lonstep for the regular longitude

    s  = varin.squeeze().shape
    varout = np.ones([np.size(lonGrid),np.size(latGrid),s[1]])*np.nan
    for i in range(nlonGrid):
        for j in range(nlatGrid):
            lon_low = lonGrid[i]-lonStep/2.0 ; lon_high = lonGrid[i]+lonStep/2.0
            lat_low = latGrid[j]-latStep/2.0 ; lat_high = latGrid[j]+latStep/2.0
            iH = np.where(((lon >= lon_low) & (lon < lon_high)) & ((lat >= lat_low) & (lat < lat_high)))
            if len(iH[0]) != 0:

                for k in range(s[1]):
                    varout[i, j, k] = np.nanmean(varin[iH[0],k])#alon[iH[iV]],alat[iH[iV]],ah2om[iH[iV]]     
    return varout.T

h2o = np.ones([days,5,45,45])*np.nan
temp = np.ones([days,5,45,45])*np.nan
theta = np.ones([days,5,45,45])*np.nan
rh = np.ones([days,5,45,45])*np.nan
z = np.ones([days,5,45,45])*np.nan
mlstime = np.ones(days)*np.nan

for i in range(days):
    ftime = itime+timedelta(days=i)
    mlstime[i] = float(ftime.strftime('%Y%m%d'))
    f = glob('/co2/hao/download_mls/nc_file/h2o/MLS-Aura_L2GP-H2O_v04_'+ftime.strftime('%Y%m%d')+'.nc')
    if len(f)!=0:
        ff = Dataset(f[0])
        h2oday = ff.variables['h2o_mix'][:,lev_start:lev_end]
        lon = ff.variables['lon'][:]
        lat = ff.variables['lat'][:]
        lev = ff.variables['Pressure'][lev_start:lev_end]
        lon[lon<0] = lon[lon<0]+360
        h2o [i] = GRID_ACE_VARS(lon,lat,longrid,latgrid,h2oday)
    f = glob('/co2/hao/download_mls/nc_file/T/MLS-Aura_L2GP-temperature_v04-20_'+ftime.strftime('%Y%m%d')+'.nc')
    if len(f)!=0:
        ff = Dataset(f[0])
        tempday = ff.variables['temperature'][:,lev_start:lev_end]
        lon = ff.variables['lon'][:]
        lat = ff.variables['lat'][:]
        lon[lon<0] = lon[lon<0]+360
        thetaday = tempday*(1000/lev[np.newaxis,:])**(2/7.)
        temp [i] = GRID_ACE_VARS(lon,lat,longrid,latgrid,tempday)
        theta [i] = GRID_ACE_VARS(lon,lat,longrid,latgrid,thetaday)
    f = glob('/co2/hao/download_mls/nc_file/RH/MLS-Aura_L2GP-RHI_v04_'+ftime.strftime('%Y%m%d')+'.nc')
    if len(f)!=0:
        ff = Dataset(f[0])
        rhday = ff.variables['rh'][:,lev_start:lev_end]
        lon = ff.variables['lon'][:]
        lat = ff.variables['lat'][:]
        lon[lon<0] = lon[lon<0]+360
        rh [i] = GRID_ACE_VARS(lon,lat,longrid,latgrid,rhday)
    else:
        print (ftime)
for i in range(5):
    z[:,i] = -7*np.log(lev[i])
import os
output = '/sn2/wyu/mls/satellite/MLS_h2o_grided_daily_'+str(year).zfill(4)+'.nc'
if os.path.exists(output):
    os.remove(output)
f2=Dataset(output,'w')
Time=f2.createDimension('Time',len(mlstime))
Lat=f2.createDimension('Lat',45)
Lon=f2.createDimension('Lon',45)
Lev=f2.createDimension('Lev',5)

Time=f2.createVariable('Time','int',('Time',))
Time.units="day"
Time.long_name="Time"
Time[:]= mlstime

Lon=f2.createVariable('Lon','f',('Lon',))
Lon.units="degree"
Lon.long_name="Longitude"
Lon[:]=longrid

Lat=f2.createVariable('Lat','f',('Lat',))
Lat.units="degree"
Lat.long_name="Latitude"
Lat[:]=latgrid

Lev=f2.createVariable('Lev','f',('Lev',))
Lev.unit="hPa"
Lev.long_name="Level"
Lev[:]=lev

H2o_mix=f2.createVariable('H2o_mix','f',('Time','Lev','Lat','Lon',))
H2o_mix.unit="ppmv"
H2o_mix.long_name="h2o mixing ratio"
H2o_mix[:]=h2o*1e6

RH=f2.createVariable('RH','f',('Time','Lev','Lat','Lon',))
RH.unit="K"
RH.long_name="relative humidity "
RH[:]=rh

Temperature=f2.createVariable('Temperature','f',('Time','Lev','Lat','Lon',))
Temperature.unit="K"
Temperature.long_name="Temperature"
Temperature[:]=temp

Theta=f2.createVariable('Theta','f',('Time','Lev','Lat','Lon',))
Theta.unit="K"
Theta.long_name="Potential Temperature"
Theta[:]=theta

Z=f2.createVariable('Z','f',('Time','Lev','Lat','Lon',))
Z.unit="K"
Z.long_name="altitude"
Z[:]=z

f2.close()
