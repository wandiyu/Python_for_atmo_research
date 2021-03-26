import xarray as xr, pandas as pd
import numpy as np    
from netCDF4 import Dataset
from glob import glob
from datetime import datetime,timedelta
from scipy.io import readsav

# --------------parameter input:----------------------------- 
# level: 1 = 121hPa, 2 = 100 hPa, 3 = 82 hPa
level=2
# start year and start month 
startyear = 2010
startmonth = 8
# date of back trajectory 
back_date = 10
# directories 
datatype = 'NEXRAD'
reanalysis = 'ERAi'
kw = '_l'
ftraj = '/mnt/data/ice2/wyu/ttl_traj_out/190810_s100_i100_ERAi_MLS_6hrly_back_201008_NA_multipoint'+kw+'/traj_s100_i100_2010_'

def read_file(infile):
	
	# Import python libraries
	import sys
	import os
	import numpy as np
	import netCDF4
	
	# Check to see if file exists
	if not os.path.isfile(infile):
		print('File "' + infile + '" does not exist.  Returning -2.')
		return -2
	
	# Check to see if file has size of zero
	if os.stat(infile).st_size == 0:
		print('File "' + infile + '" contains no valid data.  Returning -1.')
		return -1
	
	from netCDF4 import Dataset
	from netCDF4 import Variable
	# Open GridRad netCDF file
	id = Dataset(infile, "r", format="NETCDF4")
	
	# Read global attributes
	Analysis_time           = str(id.getncattr('Analysis_time'          ))
	Analysis_time_window    = str(id.getncattr('Analysis_time_window'   ))
	File_creation_date      = str(id.getncattr('File_creation_date'     ))
	Grid_scheme             = str(id.getncattr('Grid_scheme'            ))
	Algorithm_version       = str(id.getncattr('Algorithm_version'      ))
	Algorithm_description   = str(id.getncattr('Algorithm_description'  ))
	Data_source             = str(id.getncattr('Data_source'            ))
	Data_source_URL         = str(id.getncattr('Data_source_URL'        ))
	NOAA_wct_export_Version = str(id.getncattr('NOAA_wct-export_Version'))
	Authors                 = str(id.getncattr('Authors'                ))
	Project_sponsor         = str(id.getncattr('Project_sponsor'        ))
	Project_name            = str(id.getncattr('Project_name'           ))
	
	# Read list of merged files
	file_list    = (id.variables['files_merged'])[:]
	files_merged = ['']*(id.dimensions['File'].size)
	for i in range(0,id.dimensions['File'].size):
		for j in range(0,id.dimensions['FileRef'].size):
			files_merged[i] += str(file_list[i,j])
	
	# Read longitude dimension
	x = id.variables['Longitude']
	x = {'values'    : x[:],             \
		  'long_name' : str(x.long_name), \
		  'units'     : str(x.units),     \
		  'delta'     : str(x.delta),     \
		  'n'         : len(x[:])}
	
	# Read latitude dimension
	y = id.variables['Latitude']
	y = {'values'    : y[:],             \
		  'long_name' : str(y.long_name), \
		  'units'     : str(y.units),     \
		  'delta'     : str(y.delta),     \
		  'n'         : len(y[:])}
	
	# Read altitude dimension
	z = id.variables['Altitude']
	z = {'values'    : z[:],             \
		  'long_name' : str(z.long_name), \
		  'units'     : str(z.units),     \
		  'delta'     : str(z.delta),     \
		  'n'         : len(z[:])}
	
	# Read observation and echo counts
	nobs  = (id.variables['Nradobs' ])[:]
	necho = (id.variables['Nradecho'])[:]
	index = (id.variables['index'   ])[:]
	
	# Read reflectivity variables	
	Z_H  = id.variables['Reflectivity' ]
	wZ_H = id.variables['wReflectivity']

	# Create arrays to store binned values	
	values    = np.zeros(x['n']*y['n']*z['n'])
	wvalues   = np.zeros(x['n']*y['n']*z['n'])
	values[:] = float('nan')

	# Add values to arrays
	values[index[:]]  =  (Z_H)[:]
	wvalues[index[:]] = (wZ_H)[:]
	
	# Reshape arrays to 3-D GridRad domain
	values  =  values.reshape((z['n'], y['n'] ,x['n']))
	wvalues = wvalues.reshape((z['n'], y['n'] ,x['n']))

	Z_H = {'values'     : values,              \
			 'long_name'  : str(Z_H.long_name),  \
			 'units'      : str(Z_H.units),      \
			 'missing'    : float('nan'),        \
			 'wvalues'    : wvalues,             \
			 'wlong_name' : str(wZ_H.long_name), \
			 'wunits'     : str(wZ_H.units),     \
			 'wmissing'   : wZ_H.missing_value,  \
			 'n'          : values.size}
	
	# Close netCDF4 file
	id.close()
	
	# Return data dictionary	
	return {'name'                    : 'GridRad analysis for ' + Analysis_time, \
			  'x'                       : x, \
			  'y'                       : y, \
			  'z'                       : z, \
			  'Z_H'                     : Z_H, \
			  'nobs'                    : nobs, \
			  'necho'                   : necho, \
			  'file'                    : infile, \
			  'files_merged'            : files_merged, \
			  'Analysis_time'           : Analysis_time, \
			  'Analysis_time_window'    : Analysis_time_window, \
			  'File_creation_date'      : File_creation_date, \
			  'Grid_scheme'             : Grid_scheme, \
			  'Algorithm_version'       : Algorithm_version, \
			  'Algorithm_description'   : Algorithm_description, \
			  'Data_source'             : Data_source, \
			  'Data_source_URL'         : Data_source_URL, \
			  'NOAA_wct_export_Version' : NOAA_wct_export_Version, \
			  'Authors'                 : Authors, \
			  'Project_sponsor'         : Project_sponsor, \
			  'Project_name'            : Project_name}


# GridRad filter routine
def filter_data(data0):
	
	# Import python libraries
	import sys
	import os
	import numpy as np	
	
	#Extract year from GridRad analysis time string
	year = int((data0['Analysis_time'])[0:4])

	wmin        = 0.1												# Set absolute minimum weight threshold for an observation (dimensionless)
	wthresh     = 1.33 - 1.0*(year < 2009)					# Set default bin weight threshold for filtering by year (dimensionless)
	freq_thresh = 0.6												# Set echo frequency threshold (dimensionless)
	Z_H_thresh  = 18.5											# Reflectivity threshold (dBZ)
	nobs_thresh = 2												# Number of observations threshold
	
	# Extract dimension sizes
	nx = (data0['x'])['n']
	ny = (data0['y'])['n']
	nz = (data0['z'])['n']
	
	echo_frequency = np.zeros((nz,ny,nx))					# Create array to compute frequency of radar obs in grid volume with echo

	ipos = np.where(data0['nobs'] > 0)						# Find bins with obs 
	npos = len(ipos[0])											# Count number of bins with obs

	if (npos > 0):
		echo_frequency[ipos] = (data0['necho'])[ipos]/(data0['nobs'])[ipos]		# Compute echo frequency (number of scans with echo out of total number of scans)

	inan = np.where(np.isnan((data0['Z_H'])['values']))				# Find bins with NaNs 
	nnan = len(inan[0])														# Count number of bins with NaNs
	
	if (nnan > 0): ((data0['Z_H'])['values'])[inan] = 0.0

	# Find observations with low weight
	ifilter = np.where( ((data0['Z_H'])['wvalues'] < wmin       )                                            | \
							 (((data0['Z_H'])['wvalues'] < wthresh    ) & ((data0['Z_H'])['values'] <= Z_H_thresh)) |
							  ((echo_frequency           < freq_thresh) &  (data0['nobs'] > nobs_thresh)))
	
	nfilter = len(ifilter[0])									# Count number of bins that need to be removed
	
	# Remove low confidence observations
	if (nfilter > 0): ((data0['Z_H'])['values'])[ifilter] = float('nan')
	
	# Replace NaNs that were previously removed
	if (nnan > 0): ((data0['Z_H'])['values'])[inan] = float('nan')
	
	# Return filtered data0
	return data0

	
def remove_clutter(data0, **kwargs):

	# Set defaults for optional parameters
	if ('skip_weak_ll_echo' not in kwargs): skip_weak_ll_echo = 0
	
	# Import python libraries
	import sys
	import os
	import numpy as np	
	
	# Set fractional areal coverage threshold for speckle identification
	areal_coverage_thresh = 0.32
	
	# Extract dimension sizes
	nx = (data0['x'])['n']
	ny = (data0['y'])['n']
	nz = (data0['z'])['n']
	
	# Copy altitude array to 3 dimensions
	zzz = ((((data0['z'])['values']).reshape(nz,1,1)).repeat(ny, axis = 1)).repeat(nx, axis = 2)

	# First pass at removing speckles
	fin = np.isfinite((data0['Z_H'])['values'])
	
	# Compute fraction of neighboring points with echo
	cover = np.zeros((nz,ny,nx))
	for i in range(-2,3):
		for j in range(-2,3):
			cover += np.roll(np.roll(fin, i, axis=2), j, axis=1)
	cover = cover/25.0
	
	# Find bins with low nearby areal echo coverage (i.e., speckles) and remove (set to NaN).
	ibad = np.where(cover <= areal_coverage_thresh)
	nbad = len(ibad[0])
	if (nbad > 0): ((data0['Z_H'])['values'])[ibad] = float('nan')

	# Attempts to mitigate ground clutter and biological scatterers
	if (skip_weak_ll_echo == 0):
		# First check for weak, low-level echo
		inan = np.where(np.isnan((data0['Z_H'])['values']))				# Find bins with NaNs 
		nnan = len(inan[0])															# Count number of bins with NaNs
	
		if (nnan > 0): ((data0['Z_H'])['values'])[inan] = 0.0

		# Find weak low-level echo and remove (set to NaN)
		ibad = np.where(((data0['Z_H'])['values'] < 10.0) & (zzz <= 4.0))
		nbad = len(ibad[0])
		if (nbad > 0): ((data0['Z_H'])['values'])[ibad] = float('nan')
		
		# Replace NaNs that were removed
		if (nnan > 0): ((data0['Z_H'])['values'])[inan] = float('nan')

		# Second check for weak, low-level echo
		inan = np.where(np.isnan((data0['Z_H'])['values']))				# Find bins with NaNs 
		nnan = len(inan[0])															# Count number of bins with NaNs
	
		if (nnan > 0): ((data0['Z_H'])['values'])[inan] = 0.0

		refl_max   = np.nanmax( (data0['Z_H'])['values'],             axis=0)
		echo0_max  = np.nanmax(((data0['Z_H'])['values'] >  0.0)*zzz, axis=0)
		echo0_min  = np.nanmin(((data0['Z_H'])['values'] >  0.0)*zzz, axis=0)
		echo5_max  = np.nanmax(((data0['Z_H'])['values'] >  5.0)*zzz, axis=0)
		echo15_max = np.nanmax(((data0['Z_H'])['values'] > 15.0)*zzz, axis=0)

		# Replace NaNs that were removed
		if (nnan > 0): ((data0['Z_H'])['values'])[inan] = float('nan')
		
		# Find weak and/or shallow echo
		ibad = np.where(((refl_max   <  20.0) & (echo0_max  <= 4.0) & (echo0_min  <= 3.0)) | \
							 ((refl_max   <  10.0) & (echo0_max  <= 5.0) & (echo0_min  <= 3.0)) | \
							 ((echo5_max  <=  5.0) & (echo5_max  >  0.0) & (echo15_max <= 3.0)) | \
							 ((echo15_max <   2.0) & (echo15_max >  0.0)))
		nbad = len(ibad[0])
		if (nbad > 0):
			kbad = (np.zeros((nbad))).astype(int)
			for k in range(0,nz):
				((data0['Z_H'])['values'])[(k+kbad),ibad[0],ibad[1]] = float('nan')


	# Find clutter below convective anvils
	k4km = ((np.where((data0['z'])['values'] >= 4.0))[0])[0]
	fin  = np.isfinite((data0['Z_H'])['values'])
	ibad = np.where((          fin[k4km         ,:,:]          == 0) & \
							 (np.sum(fin[k4km:(nz  -1),:,:], axis=0) >  0) & \
							 (np.sum(fin[   0:(k4km-1),:,:], axis=0) >  0))
	nbad = len(ibad[0])
	if (nbad > 0):
		kbad = (np.zeros((nbad))).astype(int)
		for k in range(0,k4km+1):
			((data0['Z_H'])['values'])[(k+kbad),ibad[0],ibad[1]] = float('nan')
	
	# Second pass at removing speckles
	fin = np.isfinite((data0['Z_H'])['values'])
	
	# Compute fraction of neighboring points with echo
	cover = np.zeros((nz,ny,nx))
	for i in range(-2,3):
		for j in range(-2,3):
			cover += np.roll(np.roll(fin, i, axis=2), j, axis=1)
	cover = cover/25.0
			
	# Find bins with low nearby areal echo coverage (i.e., speckles) and remove (set to NaN).
	ibad = np.where(cover <= areal_coverage_thresh)
	nbad = len(ibad[0])
	if (nbad > 0): ((data0['Z_H'])['values'])[ibad] = float('nan')
	
	return data0

def encountering_calculation(level,startyear,startmonth,back_date,datatype,ftraj,kw):
# -------------------read data -------------------------
    fMLS = xr.open_dataset('/mnt/data/ice2/wyu/mls_l2/MLS_L2_'+str(startyear).zfill(4)+str(startmonth).zfill(2)+'_NA.nc')
# ----------------initiate values------------------------------ 
# level 
    levs = np.array([146,121,100,82])
# time parameters 
    itime = datetime(startyear,startmonth,1,0)
    dateinit = itime-timedelta(days=back_date)
    dateend = datetime(startyear,startmonth+1,1,0)
    days = (dateend-dateinit).days

# outputs
    parcel_lons = np.ones([len(fMLS.index)*9,days*24])*np.nan
    parcel_lats = np.ones([len(fMLS.index)*9,days*24])*np.nan 
    parcel_levs = np.ones([len(fMLS.index)*9,days*24])*np.nan 
    parcel_temps = np.ones([len(fMLS.index)*9,days*24])*np.nan 
    MLS_convz = np.ones([len(fMLS.index)*9,days*24])*np.nan
    MLS_pfister_convth = np.ones([len(fMLS.index)*9,days*24])*np.nan
    backtrajtime = np.ones(days*24)
    for i in range(len(backtrajtime)):
        ftime = dateinit+timedelta(hours=i)
        backtrajtime[i] = float(ftime.strftime('%Y%m%d%H'))
    
# -----------------calculating: ---------------------------------------
# the position of all parcels at each time 
    f = glob(ftraj+'*.sav')
    for i in range(len(fMLS.index)):
        ftime = datetime(int(fMLS.year[i].values),int(fMLS.month[i].values),\
                int(fMLS.day[i].values), int(np.round(fMLS.hour[i].values)))
        f = glob(ftraj+\
             ftime.strftime("%d%H")+'_I'+str(int(fMLS.index[i].values)).zfill(6)+'_E?.sav')
        for i_file in range(9):
            ff = readsav(f[i_file])
    # Delta time is the difference between 2011 07 01 00 and the observation time 
            Deltatime = (ftime - itime).days*24+(ftime - itime).seconds//3600
            parcel_lons[i*9+i_file,Deltatime:Deltatime+back_date*24] = ff.lonrec[:back_date*24][::-1]
            parcel_lats[i*9+i_file,Deltatime:Deltatime+back_date*24] = ff.latrec[:back_date*24][::-1]
            parcel_levs[i*9+i_file,Deltatime:Deltatime+back_date*24] = ff.prrec[:back_date*24][::-1]
            parcel_temps[i*9+i_file,Deltatime:Deltatime+back_date*24] = ff.temprec[:back_date*24][::-1]
# convective height at each postion 
    if datatype == 'NEXRAD':
        convection_dir = '/NEXRAD/v3_1/level3/3d_filtered/'
        for j in range(24*days):
            i_finite = np.where(np.isfinite(parcel_lats[:,j]))[0]
            if len(i_finite)>0:
                datef = dateinit+timedelta(hours=j)
                print (datef,datetime.now())
                fconv = glob(convection_dir+str(datef.year).zfill(4)+'/'+datef.strftime('%Y%m')+'/nexrad_3d_filtered_v3_1_'+\
                     datef.strftime('%Y%m%d')+'T'+str(datef.hour).zfill(2)+'0000Z.nc')
                if len(fconv)>0:
                    out1 = read_file(fconv[0])
                    if out1 != -1:
                        out1 = filter_data(out1)
                        out1 = remove_clutter(out1)
                        convz = xr.DataArray(out1['Z_H']['values'], dims = ['alt','lat','lon'], coords = [out1['z']['values'],\
                                                                        out1['y']['values'],out1['x']['values']])
                        for i_file in range(len(i_finite)):
                            #conv_ver = convz.sel(lat = parcel_lats[i_finite[i_file],j], lon = parcel_lons[i_finite[i_file],j]
                            #         ,method = 'nearest').values
                            conv_ver = convz.sel(lat=slice(parcel_lats[i_finite[i_file],j]-0.125,parcel_lats[i_finite[i_file],j]+0.125),\
                                    lon=slice(parcel_lons[i_finite[i_file],j]-0.125,parcel_lons[i_finite[i_file],j]+0.125))
                            if len(conv_ver)!=0:
                                conv_ver = conv_ver.mean(dim=['lat','lon'])
                                if (len(np.where(conv_ver>10)[0]))>0:
                                    MLS_convz[i_finite[i_file],j] = convz.alt[np.where(conv_ver>10)[0][-1]]
    if datatype == 'Pfister':
        file_interval = 3
        convection_dir = '/mnt/data/co3/schoeberl/pfister_conv/'
        f0 = readsav(convection_dir+'cloudaltthetkm_2001010100_reanmix_noice_ocean_trmm.15.9_srch.3_offs1.00_tropw0.70.sav')
        convlat = f0.rainlat
        convlon = f0.rainlon
            
        for j in range(24*days//file_interval):
            i_finite = np.where(np.isfinite(parcel_lats[:,j*file_interval]))[0]
            if len(i_finite)>0:
                datef = dateinit+timedelta(hours=j*file_interval)
                fconv = glob(convection_dir+'cloudaltthetkm_'+datef.strftime('%Y%m%d%H')+'_reanmix_noice_ocean_trmm.15.9_srch.3_offs1.00_tropw0.70.sav')
                if len(fconv)>0:
                    ff = readsav(f[0])
                    fconv = readsav(fconv[0])
                    convth = xr.DataArray(fconv.rainthet, dims = ['lat','lon'], coords = [convlat,convlon])
                    convz = xr.DataArray(fconv.rainalt, dims = ['lat','lon'], coords = [convlat,convlon])
                    for i_file in range(len(i_finite)):
                        MLS_convz[i_finite[i_file],j*file_interval] = (convz.sel(lat = parcel_lats\
                            [i_finite[i_file],j*file_interval], lon = parcel_lons[i_finite[i_file],j*file_interval]
                                      ,method = 'nearest').values)
                        MLS_pfister_convth[i_finite[i_file],j*file_interval] = (convth.sel(lat = parcel_lats\
                         [i_finite[i_file],j*file_interval], lon = parcel_lons[i_finite[i_file],j*file_interval]
              	      ,method = 'nearest').values)
# ----------------------Data saving-------------------
    for i_file in range(9):
        f2=Dataset('/mnt/data/ice2/wyu/back_traj_result_test/MLS_backtrajectories_init'+str(int(levs[level]))+\
               '_'+datatype+'_'+reanalysis+'_'+str(startyear).zfill(4)\
               +str(startmonth).zfill(2)+'_10days_25_'+'E'+str(i_file+1).zfill(1)+kw+'.nc','w')
        Time=f2.createDimension('Time',None)
        Index = f2.createDimension('Index',len(fMLS.index))

        Time=f2.createVariable('Time','int',('Time',))
        Time.units="hour"
        Time.long_name="Time"
        Time[:]= backtrajtime

        Index=f2.createVariable('Index','f',('Index',))
        Index.units=""
        Index.long_name="Longitude"
        Index[:]=fMLS.index


        Lon=f2.createVariable('Lon','f',('Index','Time',))
        Lon.unit="degree"
        Lon.long_name="longitudes"
        Lon[:]=parcel_lons[i_file::9]

        Lat=f2.createVariable('Lat','f',('Index','Time',))
        Lat.unit="degree"
        Lat.long_name="latitude"
        Lat[:]=parcel_lats[i_file::9]
    
        Lev=f2.createVariable('Lev','f',('Index','Time',))
        Lev.unit="degree"
        Lev.long_name="altitude"
        Lev[:]=parcel_levs[i_file::9]
        
        Temp=f2.createVariable('Temp','f',('Index','Time',))
        Temp.unit="degree"
        Temp.long_name="temperature"
        Temp[:]=parcel_temps[i_file::9]

        Conv_z=f2.createVariable('Conv_z','f',('Index','Time',))
        Conv_z.unit="km"
        Conv_z.long_name="convection altitude of each position"
        Conv_z[:]=MLS_convz[i_file::9]
    
        if datatype == 'Pfister':
            Conv_th=f2.createVariable('Conv_th','f',('Index','Time',))
            Conv_th.unit="km"
            Conv_th.long_name="convection potential temperature of each position"
            Conv_th[:]=MLS_pfister_convth

    fMLS.close()
    f2.close()

encountering_calculation(level,startyear,startmonth,back_date,datatype,ftraj,kw)
