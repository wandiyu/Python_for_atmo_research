PRO CALL_MODEL
  pr1x = 100
  infile = '/mnt/data/ice2/wyu/mls_l2/MLS_L2_201608_NA.nc'
  savFolder = '191029_s100_i'+STRING(pr1x, format='(I3)')+'_ERAi_MLS_6hrly_back_201608_NA_multipoint/'
  id = NCDF_OPEN(infile)
  NCDF_VARGET, id, 'year', year
  NCDF_VARGET, id, 'month', month
  NCDF_VARGET, id, 'day', day
  NCDF_VARGET, id, 'hour', hour
  NCDF_VARGET, id, 'longitude', longitude
  NCDF_VARGET, id, 'latitude', latitude
  NCDF_VARGET, id, 'index', findex
  NCDF_CLOSE, id
  s=SIZE(year)
  nn = 0
  FOR iday = 1,31 DO BEGIN 
    For ihour = 0,23 DO BEGIN
	idx = WHERE(ROUND(day) EQ iday AND ROUND(hour) EQ ihour, nidx)
    nn=nn+nidx
	PRINT, iday,ihour,nidx,nn
	IF nidx GE 1 THEN BEGIN 	 
	traj_erai_back_multiple,savFolder,pr1x, year[idx[0]],month[idx[0]],day[idx[0]],ROUND(hour[idx[0]]),$
[longitude[idx],longitude[idx],longitude[idx],longitude[idx]-0.25,longitude[idx]+0.25,longitude[idx]-0.25,longitude[idx]+0.25,longitude[idx]-0.25,longitude[idx]+0.25],$
[latitude[idx],latitude[idx]+0.25,latitude[idx]-0.25,latitude[idx],latitude[idx],latitude[idx]+0.25,latitude[idx]+0.25,latitude[idx]-0.25,latitude[idx]-0.25],$
[findex[idx],findex[idx],findex[idx],findex[idx],findex[idx],findex[idx],findex[idx],findex[idx],findex[idx]],$
[findex[idx]*0+1,findex[idx]*0+2,findex[idx]*0+3,findex[idx]*0+4,findex[idx]*0+5,findex[idx]*0+6,findex[idx]*0+7,findex[idx]*0+8,findex[idx]*0+9] 
	ENDIF 
	ENDFOR 
 ENDFOR
END
