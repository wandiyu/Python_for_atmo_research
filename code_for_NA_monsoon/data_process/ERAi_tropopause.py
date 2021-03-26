import xarray as xr, numpy as np, pandas as pd
from netCDF4 import Dataset
from glob import glob
from datetime import datetime,timedelta
from scipy.io import readsav
from scipy.interpolate import CubicSpline,interp1d

years = [2005+i for i in range(12)]

months = [str(year)+str(j+6).zfill(2) for year in years for j in range(3)]


def get_tropopause(month):
    print (month)
    tempfolder = '/co1/hao/reanlys_raw2nc/ECMWF_ERAi/temperature/'
    fs = glob(tempfolder+'ECMWF_ERAi_hybrid_6hr_t_'+month+'??.nc')
    tropopause_z = np.zeros([len(fs),181,360])
    tropopause_pr = np.zeros([len(fs),181,360])
    pressure_new = np.arange(60,300,5)
    znew = -7*np.log(pressure_new/1000)
    for i in range(len(fs)):
        ftemp = xr.open_dataset(fs[i])
        #ftemp = ftemp.isel(levels=slice(0,41)).sel(longitude=slice(240,290)\
        #       ,latitude=slice(0,60)).mean(dim='time')
        ftemp = ftemp.isel(levels=slice(20,41)).mean(dim='time')
        vertical_temp = ftemp.t.values
        vertical_theta = ftemp.theta.values
        vertical_pressure = 1000.*(vertical_temp/vertical_theta)**(7./2.)
        #vertical_z = -7*np.log(vertical_pressure/1000)
        #lapse_rate = ((vertical_temp[:-1]-vertical_temp[1:])/(vertical_z[1:]-vertical_z[:-1]))
        s = tropopause_z.shape
        for i_lat in range(s[1]):
            for i_lon in range(s[2]):
        #for i_lat in [90]:
        #    for i_lon in range(s[2]):
                cs = interp1d(vertical_pressure[:,i_lat,i_lon],vertical_temp[:,i_lat,i_lon])
                tnew = cs(pressure_new)
                lapse_rate = ((tnew[:-1]-tnew[1:])/(znew[1:]-znew[:-1]))
                if len(np.where(lapse_rate<2)[0]):
                    tropopause_z[i,i_lat,i_lon] = (znew[np.where(lapse_rate<2)[0][-1]+1])
                    tropopause_pr[i,i_lat,i_lon] = (pressure_new[np.where(lapse_rate<2)[0][-1]+1])
                else:
                    tropopause_z[i,i_lat,i_lon] = (znew[np.where(tnew==np.nanmin(tnew))])
                    tropopause_pr[i,i_lat,i_lon] = (pressure_new[np.where(tnew==np.nanmin(tnew))])
    output = ('/mnt/data/ice2/wyu/tropopause_height_fine_grid/tropopause_height_{}_linearreg_global.nc'.format(month))
    #os.remove(output)
    f2=Dataset(output,'w')
    Time=f2.createDimension('Time',None)
    Lat=f2.createDimension('Lat',181)
    Lon=f2.createDimension('Lon',360)

    Time=f2.createVariable('Time','int',('Time',))
    Time.units="day"
    Time.long_name="Time"
    Time[:]= np.array([int(month+str(i+1).zfill(2)) for i in range(len(fs))])


    Lon=f2.createVariable('Lon','f',('Lon',))
    Lon.units="degree"
    Lon.long_name="Longitude"
    Lon[:]=np.arange(360)

    Lat=f2.createVariable('Lat','f',('Lat',))
    Lat.units="degree"
    Lat.long_name="Latitude"
    Lat[:]=np.arange(181)-90

    Tropopause_z=f2.createVariable('Tropopause_z','f',('Time','Lat','Lon',))
    Tropopause_z.unit="km"
    Tropopause_z.long_name="tropopause height"
    Tropopause_z[:]=tropopause_z

    Tropopause_pr=f2.createVariable('Tropopause_pr','f',('Time','Lat','Lon',))
    Tropopause_pr.unit="hPa"
    Tropopause_pr.long_name="tropopause pressure"
    Tropopause_pr[:]=tropopause_pr
    f2.close()



for month in months:
    get_tropopause(month)
