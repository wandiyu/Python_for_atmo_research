PRO CALL_MODEL2
	year = 2016
	month = 8
	pr1x = 100
	savFolder = '191005_s100_i'+STRING(pr1x, format='(I3)')+'_ERAi_MLS_6hrly_back_201608_NA/'
	lonRange=[240.0, 290.0] & lonStep = 2.0
    latRange=[10.0, 50.0] & latStep =  1.0
    nlon = UINT((lonRange[1] - lonRange[0])/lonStep)
    nlat = UINT((latRange[1] - latRange[0])/latStep)
    lonx=FINDGEN(nlon)*lonStep+lonRange(0)+lonStep/2.
    latx=FINDGEN(nlat)*latStep+latRange(0)+latStep/2.

	; injection grid of lat in equal area
    latStep=(SIN(latRange(1)*!PI/180.0)-SIN(latRange(0)*!PI/180.0))/nlat
    latx=SIN(latRange(0)*!PI/180.0)+FINDGEN(nlat)*latStep+latStep/2.0
    latx=ASIN(latx)*180.0/!PI
	
	; now mash the grid
    INDEXTO2D,lonx,latx,olon,olat  ;; lat oriented
    ;indexto2d,latx,lonx,olon,olat    ;; lon oriented
	findex = findgen(N_ELEMENTS(lonx)*N_ELEMENTS(latx))
	print, size(olon)
	print, N_ELEMENTS(lonx)*N_ELEMENTS(latx)
	FOR iday = 1,31 DO BEGIN 
		traj_erai_back_multiple,savFolder,pr1x, year, month, iday,0, olon,olat,findex  
	ENDFOR
END 
