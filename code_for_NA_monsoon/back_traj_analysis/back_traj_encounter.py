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

def SAT(temp,pr): #calculate the saturation h2o mixing ratio 
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
    h2o=esi/(pr-esi)
    high=np.where((esi*100-pr)<0)
    nh=np.size(high)
    if nh>0:
        h2o[high]=(esi/pr)[high]
    return h2o*1.0e6

def back_traj_encounter(fMLS,year,month):
    level = 2
    lev_diff = 0
    i_backtime = 0
#----------read data--------------------    
    #nexrad_conv = xr.open_dataset('/sn2/wyu/mls/convection/NEXRAD_convection_occurance_grided\
    #_hourly_{}{}_pressure_level.nc'.format(year,month))
    #mls=xr.open_dataset('/sn2/wyu/research/mls_mon/20170419_mls_h2o_monthly.nc')
    #mls['Lev']=np.floor(mls.Lev)
    month2 = str(int(month)-1).zfill(2)
    ftropopause0 = xr.open_dataset('/mnt/data/ice2/wyu/tropopause_height_fine_grid/tropopause_height_{}{}_linearreg_global.nc'.format(year,month2))
    ftropopause1 = xr.open_dataset('/mnt/data/ice2/wyu/tropopause_height_fine_grid/tropopause_height_{}{}_linearreg_global.nc'.format(year,month))
    ftropopause = xr.concat([ftropopause0,ftropopause1],dim='Time')
    fbackfile = '/mnt/data/ice2/wyu/back_traj_result_test/MLS_backtrajectories_init100_NEXRAD_ERAi_{}{}_'.format(year,month)
    fbackfiles = glob(fbackfile+'10days_25_E*.nc')
    year,month = int(year),int(month)
    tempfolder = '/co1/hao/reanlys_raw2nc/ECMWF_ERAi/temperature/'
    
#------------initiate data-----------------------
    itime = datetime(year,7,1,0,0)-timedelta(days=10)
    levs = np.array([146,121,100,82])[level]
    i_overs = []
    levs = []
    lats = []
    lons = []
    temps = []
    tempmins = []
    h2omins = []
    convtimes = []
    tropopauses = []
    convlats = []
#----------find out encounter events ------------------
    for i_file in range(27):
        fback = xr.open_dataset(fbackfiles[i_file])
        #fback = fback.isel(Index=np.where(fMLS.longitude>=240)[0])
        MLS_n_convz = fback.Conv_z.values
        MLS_n_convz[fback.Lat>49] = np.nan
        MLS_n_convz[fback.Lat<25] = np.nan
        MLS_n_convz[fback.Lon<245] = np.nan
        MLS_n_convz[fback.Lat>290] = np.nan
        z = -7*np.log(fback.Lev.values/1000)
        #z = take_half(z)
        s = z.shape
        for i in range(s[0]):
            i_start = np.where(np.isfinite(z[i]))[0][-1]
            i_end = np.where(np.isfinite(z[i]))[0][0]
            z[i][i_end:i_end+24*i_backtime] = np.nan
        i_over = np.where(np.nanmax((MLS_n_convz-z+lev_diff)>=0,axis=1))[0]
        #print (len(set(i_over)))
#-----------make sure it is above the tropopause
        #print (len(i_over))
        if len(i_over)!=0:
            tropopause_pressure = np.zeros(len(i_over))
            backtrajtime = fback.Time.values
            lev = np.zeros(len(i_over))
            lat = np.zeros(len(i_over))
            lon = np.zeros(len(i_over))
            temp = np.zeros(len(i_over))
            tempmin = np.zeros(len(i_over))
            h2omin = np.zeros(len(i_over))
            convtime = np.zeros(len(i_over))
            convlat = np.zeros(len(i_over))
            for j in range(len(i_over)):
#------------calculate tropopause height----------------------
                i=i_over[j]
                encounter_time = backtrajtime[np.where((fback.Conv_z.values-z+lev_diff)[i]>=0)[0][-1]]
                tropopause_pressure[j] = ftropopause.Tropopause_pr.sel(Time=encounter_time//100,\
                            Lon=fback.Lon.values[i,np.where((fback.Conv_z.values-z+lev_diff)[i]>=0)[0][-1]]\
               ,Lat=fback.Lat.values[i,np.where((fback.Conv_z.values-z+lev_diff)[i]>=0)[0][-1]],method='nearest')
    #----------------compare tropopause height with parcel position --------------------
                lev[j] = fback.Lev.values[i,np.where((MLS_n_convz-z+lev_diff)[i]>=0)[0][-1]]
                lat[j] = fback.Lat.values[i,np.where((MLS_n_convz-z+lev_diff)[i]>=0)[0][-1]]
                lon[j] = fback.Lon.values[i,np.where((MLS_n_convz-z+lev_diff)[i]>=0)[0][-1]]
                i_mintemp = [np.where(fback.Temp[i].values==np.nanmin(fback.Temp.values[i,np.where((MLS_n_convz-z+lev_diff)[i]>=0)[0][-1]:]))[0][0]]
                i_observ = np.where(np.isfinite(fback.Temp[i].values))[0][-1]
                temp[j] = i_observ-i_mintemp
                tempmin[j] = fback.Temp[i,i_mintemp].values
                h2omin[j] = SAT(fback.Temp[i,i_mintemp].values,fback.Lev[i,i_mintemp].values)
                convlat[j] = fback.Lat[i,i_mintemp].values
                convtime[j] =np.where(np.isfinite(z[i]))[0][-1] -np.where((fback.Conv_z.values-z+lev_diff)[i]>=0)[0][-1]
            i_over_trop = np.where((tropopause_pressure-lev)>0)[0]
            #i_over_trop = np.arange(len(i_over))
            i_over = i_over[i_over_trop]
            levs = np.append(levs,lev[i_over_trop])
            lons = np.append(lons,lon[i_over_trop])
            lats = np.append(lats,lat[i_over_trop])
            temps = np.append(temps,temp[i_over_trop])
            tempmins = np.append(tempmins,tempmin[i_over_trop])
            h2omins = np.append(h2omins,h2omin[i_over_trop])
            tropopauses = np.append(tropopauses,tropopause_pressure[i_over_trop])
            convtimes  = np.append(convtimes,convtime[i_over_trop])
            convlats = np.append(convlats,convlat[i_over_trop])
        i_overs = np.append(i_overs,i_over)
    i_over = i_overs.astype(int)

    return (i_over,lats,lons,levs,temps,tempmins,h2omins,tropopauses,convtimes,convlats)  

months = [str(i+6).zfill(2) for i in range(3)]
years = [str(i+2005) for i in range(12)]
lonGrid = np.arange(230,291,2)
latGrid = np.arange(25,51,2)
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
