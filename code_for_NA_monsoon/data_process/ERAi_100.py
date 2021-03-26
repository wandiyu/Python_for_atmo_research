import xarray as xr, numpy as np, pandas as pd
from netCDF4 import Dataset
from glob import glob
from datetime import datetime,timedelta
from scipy.interpolate import CubicSpline,interp1d

years = [2005+i for i in range(12)]

months = [str(year)+str(j+6).zfill(2) for year in years for j in range(3)]

def get_uv_100hPa(month):
    print (month)
    windFolder = '/co1/hao/reanlys_raw2nc/ECMWF_ERAi/winds/'
    fs = glob(windFolder+'ECMWF_ERAi_hybrid_6hr_UVQ_'+month+'??.nc')
    tempfolder = '/co1/hao/reanlys_raw2nc/ECMWF_ERAi/temperature/'
    fs2 = glob(tempfolder+'ECMWF_ERAi_hybrid_6hr_t_'+month+'??.nc')
    u100 = np.zeros([len(fs2),26,61])
    v100 = np.zeros([len(fs2),26,61])
    for i in range(len(fs2)):
        fwind = xr.open_dataset(fs[i]).sel(latitude=slice(25,50),longitude=slice(230,290)).isel(levels=slice(20,41)).mean(dim='time')
        ftemp = xr.open_dataset(fs2[i]).sel(latitude=slice(25,50),longitude=slice(230,290)).isel(levels=slice(20,41)).mean(dim='time')
        pressure = 1000.*(ftemp.t/ftemp.theta)**(7./2.)
        u = fwind.u
        v = fwind.v
        s = u100.shape
        for i_lat in range(s[1]):
            for i_lon in range(s[2]):
                cs = CubicSpline(pressure[:,i_lat,i_lon],u[:,i_lat,i_lon])
                u100[i,i_lat,i_lon] = cs(100)
                cs = CubicSpline(pressure[:,i_lat,i_lon],v[:,i_lat,i_lon])
                v100[i,i_lat,i_lon] = cs(100)
    x = fwind.longitude.values
    y = fwind.latitude.values
    
    output = ('/mnt/data/ice2/wyu/wind_100/wind_100hPa_{}_NA.nc'.format(month))
    #os.remove(output)
    f2=Dataset(output,'w')
    Time=f2.createDimension('Time',None)
    Lat=f2.createDimension('Lat',26)
    Lon=f2.createDimension('Lon',61)

    Time=f2.createVariable('Time','int',('Time',))
    Time.units="day"
    Time.long_name="Time"
    Time[:]= np.array([int(month+str(i+1).zfill(2)) for i in range(len(fs2))])


    Lon=f2.createVariable('Lon','f',('Lon',))
    Lon.units="degree"
    Lon.long_name="Longitude"
    Lon[:]=x

    Lat=f2.createVariable('Lat','f',('Lat',))
    Lat.units="degree"
    Lat.long_name="Latitude"
    Lat[:]=y

    U=f2.createVariable('U','f',('Time','Lat','Lon',))
    U.unit="m/s"
    U.long_name="U"
    U[:]=u100

    V=f2.createVariable('V','f',('Time','Lat','Lon',))
    V.unit="m/s"
    V.long_name="V"
    V[:]=v100
    
    f2.close()

def get_t_100hPa(month):
    print (month)
    tempfolder = '/co1/hao/reanlys_raw2nc/ECMWF_ERAi/temperature/'
    fs2 = glob(tempfolder+'ECMWF_ERAi_hybrid_6hr_t_'+month+'??.nc')
    t100 = np.zeros([len(fs2),21,61])
    for i in range(len(fs2)):
        ftemp = xr.open_dataset(fs2[i]).sel(latitude=slice(0,20),longitude=slice(230,290)).isel(levels=slice(20,41)).mean(dim='time')
        pressure = 1000.*(ftemp.t/ftemp.theta)**(7./2.)
        t = ftemp.t
        s = t100.shape
        for i_lat in range(s[1]):
            for i_lon in range(s[2]):
                cs = CubicSpline(pressure[:,i_lat,i_lon],t[:,i_lat,i_lon])
                t100[i,i_lat,i_lon] = cs(100)
    x = ftemp.longitude.values
    y = ftemp.latitude.values
    
    output = ('/mnt/data/ice2/wyu/ERAi_100/temp_100hPa_{}_tro.nc'.format(month))
    #os.remove(output)
    f2=Dataset(output,'w')
    Time=f2.createDimension('Time',None)
    Lat=f2.createDimension('Lat',21)
    Lon=f2.createDimension('Lon',61)

    Time=f2.createVariable('Time','int',('Time',))
    Time.units="day"
    Time.long_name="Time"
    Time[:]= np.array([int(month+str(i+1).zfill(2)) for i in range(len(fs2))])


    Lon=f2.createVariable('Lon','f',('Lon',))
    Lon.units="degree"
    Lon.long_name="Longitude"
    Lon[:]=x

    Lat=f2.createVariable('Lat','f',('Lat',))
    Lat.units="degree"
    Lat.long_name="Latitude"
    Lat[:]=y

    T=f2.createVariable('T','f',('Time','Lat','Lon',))
    T.unit="K"
    T.long_name="Temperature"
    T[:]=t100

    
    f2.close()

for month in months:
    get_t_100hPa(month)

for month in months:
    get_uv_100hPa(month)
