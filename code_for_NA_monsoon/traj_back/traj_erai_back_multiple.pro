;;; info on variables in save files
;; lat/lon: latitude/longitude of parcel
;; zp: pressure altitude (pr=1000*exp(-zp/7))
;; h2o_mix: mixing ratio for water w/o methane oxidation
;; h2o_meth_mix: mixing ratio for water w/ methane oxidation
;; ch4f: methane mixing ratio for the h2o_meth_mix parcel
;; tem: temperature of parcel
;; index: number of days between start of run and day parcel was initialized
;; parcel_unique: a unique identifier for each parcel
;; parcel_un_lat,parcel_un_lon,parcel_un_press: lat/lon/pressure of initial location of the parcel
;; lat380/lon380: latitude and longitude of where the parcel ascended through the 380 K surface
;; convlat/convlon/convpress: latitude/longitude/pressure of the last convective intersection of the trajectory
;; if convpress is negative, then it means that the h2o_meth parcel was supersaturated when convection hit it
;; convindex: time of convective event in days since run started

;; fdhm_cnvSat: if conv_sat is negative, then it means that the h2o_meth parcel was supersaturated when convection hit it
;; fdh_cnvSat: if conv_sat is negative, then it means that the h2o parcel was supersaturated when convection hit it

;; for h2o_mix parcel:
;; fdh_lat/fdh_lon/fdh_press/fdh_temp/fdh_h2o: lat/lon/pressure/temp/dehydrated_h2o of last dehydration location
;; if pressure is negative, then it means that the event was convection into supersaturated air
;; fdh_index: time of last dehydration event, in days since run started
;; for h2o_meth_mix parcel: same variables, but with fdhm instead of fdh

;; if record all dehydration events (b_rec_dhALL=1)
;; for h2o_mix parcel:
;; dh_lat/dh_lon/dh_press/dh_temp/dh_xh2o: lat/lon/pressure/temp/dehydrated(excess)_h2o of all dehydration events wh
;; dh_h2o: after dehydration, the air parcel's water vapor amount
;; dh_un_lat/lat/press: where those dehydrated parcels originally initiated at
;; for h2o_meth_mix parcel: same variables, but with dhm instead of dh for parcels w/ methane

;; if track dhy fro each parcel (b_track_pacl_dh=1)
;; for each parcel, their dehydration history
;; pacl.dh_lon/dh_lat/dh_press/dh_temp/dh_xh2o/dh_h2o

;; Modify history:
;;  Updating the tropospheric water vapor using the reanalysis data. 
;;  11.06.2017 Wandi Yu

PRO TRAJ_ERAi_BACK_MULTIPLE,savFolder,pr1x, year1x,month1x,day1x,hour1x,lon1x,lat1x,findex,poindex
    SPAWN, 'hostname'
    PRINT, SYSTIME()
    PRINT, '******************** running start here ***************'
    PRINT
    
    sep = PATH_SEP()
    recday = 10
    ;year1x = 2004
    ;month1x = 03
    ;day1x = 29
    ;hour1x = 23
    ;lat1x = 46.85
    ;lon1x = -86.78+360
    ;pr1x = 100
	b_MLS= 1     
   	 
    supersat=1.00; 1.00  ; saturation level
    
    th_top=5000.  ;; remove parcels above this
    th_bottom=300.;; remove parcels below this
 
    ntdays  =  24             ;Number of timesteps per day; 86400s/ntdays=timeStep in seconds
    ntsave  =  ntdays         ;Number of saves per day
    ntplot  =  1 
    injecting_intvl = 2
    stopinjecting = 1
    firststep = 1
    initwater=200.
    b_MERRA_run=0 & b_ERAi_run=1 & b_CFSR_run=0 & b_MERmrgERAi_run=0 & b_JRA25_run=0
    b_GEOSCCM_run=0 & b_GEOSCCM_anvil_run=0
    b_WAC_ctl_run=0 & b_WAC_4xco2_run=0 & b_WAC_sst4k_run=0 & b_WAC_run=0 & b_WAC_ice_run=0
    

    satname=STRCOMPRESS('s'+STRING(FIX(supersat*100),'(I4)')+'_'+'i', /remove_all)
    
    satname=satname+STRING(pr1x, format='(I3)')

    
    ; the folder that will hold this run's results. This folder is always under '...ANDY_DIR/ttl_traj_out/'
    ;savFolder='190425_'+satname+'_ERAi_MLS_6hrly_back_201107_NA/' ; make sure to have a '/' after the name

 
    mdl_aux_data_dir='/atmomounts/home/grad/wyu/idl/traj3D_2_3/traj_mdl_aux_data/'
    cnv_freq_file=mdl_aux_data_dir+'new.convection.freq'
    ch4_oxid_para=mdl_aux_data_dir+'prodloss2d.sav'
    OLR_file=mdl_aux_data_dir+'olr.day.mean_20121231.nc'

    
    b_daily=0              ; 1=turn on daily
    b_6hrly=1             ; 1=turn on 6hrly
    
    IF b_daily EQ 1 AND b_6hrly EQ 1 THEN BEGIN
        PRINT, 'WRONG! YOU CANNOT RUN DAILY also 6HRLY! ' & STOP
    ENDIF
    IF b_daily EQ 0 AND b_6hrly EQ 0 THEN BEGIN
        PRINT, 'WRONG! YOU SHOULD RUN either DAILY or 6HRLY! ' & STOP
    ENDIF
    
    judgement=[b_MERRA_run, b_ERAi_run, b_CFSR_run, b_MERmrgERAi_run, b_GEOSCCM_run,b_GEOSCCM_anvil_run, b_JRA25_run, $
        b_WAC_ctl_run, b_WAC_4xco2_run, b_WAC_sst4k_run, b_WAC_run, b_WAC_ice_run]
        
    idx=WHERE(judgement EQ 1, nidx)
    IF nidx GT 1 THEN BEGIN
        PRINT, 'you cannot run the model with different wind and T setup at the same time'
        STOP
    ENDIF
    
    
    ; ERAi run
    IF b_ERAi_run EQ 1 THEN BEGIN
        IF b_daily EQ 1 THEN windDir='/co2/reanlys_wind_isen/ECMWF_ERAi/daily/'             ; only u, v
        IF b_6hrly EQ 1 THEN windDir='/co2/reanlys_wind_isen/ECMWF_ERAi/6hrly/'             ; only u, v
        tempDir = '/co1/hao/reanlys_raw2nc/ECMWF_ERAi/temperature/'
    ENDIF
 
	IF b_MERRA_run EQ 1 THEN BEGIN
        IF b_daily EQ 1 THEN windDir='/co2/reanlys_wind_isen/ECMWF_ERAi/daily/'             ; only u, v
        IF b_6hrly EQ 1 THEN windDir='/mnt/data/ice1/reanlys_wind_isen/MERRA2_official/6hrly_hybrid/'
        tempDir = '/mnt/data/ice1/reanlys_raw2nc/MERRA2_official/T_6hrly_hybrid/'
    ENDIF
   
    outDir= '/mnt/data/ice2/wyu/ttl_traj_out/'+savFolder
    IF DIR_EXIST(outDir) NE 1 THEN FILE_MKDIR, outDir
    
    
    ch4_delta=(1.8-1.54)/31.  ; from 1979 to 2012, the increment
    n_ch4_yr=2100L-1979L+1L
    ch4_1979_2012=1.54+DINDGEN(n_ch4_yr)*ch4_delta
    ch4_actl=MAKE_ARRAY(1979L+n_ch4_yr, /DOUBLE)
    ch4_actl[1979:2100]=ch4_1979_2012
   
	;load convection frequency
    RESTORE,cnv_freq_file
    ;return latH[36]: -87.5, -82.5, ... 87.5,           5 degree interval
    ;           olrH[26]: 95, 105, 115, 125, ..., 345,    10 K interval
    ;               pH[6]: 215, 178, 146, 121, 100, 83   hPa, 12 levels/decade change
    ;            totV[26,36,6,12]
    dph=(ph-SHIFT(ph,1))(1:*) ; -37,-32,-25,-21,-17
    
    ; load NOAA OLR data
    IF N_ELEMENTS(olr) LT 100 THEN READDAILYOLR,lono,lato,olrday,olr, fname=OLR_file
    ; return   lono[144]: 0,2.5,5.0,7.5,...,357.5
    ;                lato[73]: 90,87.5,...,2.5,0,-2.5,...,-87.5, -90
    ;           day[13xxx]: 17298624.0, 17298648.0, ...
    ;olr[144,73,13xxx]: lon x lat x day
    
    ;Constants for TRAJ3D programs
    DEFSYSV, '!xymax',   2.0, 1  ;
    DEFSYSV, '!nlimit',  89.0, 1  ;
    DEFSYSV, '!slimit', -89.0, 1  ;
    
   IF hour1x EQ 24 THEN BEGIN
		datename = MAKE_DATE(year1x, month1x, day1x, 0, 0, 0)
		datename = TIME_INC(datename, !TauDSolar)
		datef = datename
	print, datename
	ENDIF ELSE BEGIN
	datename  = MAKE_DATE(year1x, month1x, day1x, hour1x, 0, 0)
	starthour = LONG(hour1x/6)*6
	datef = MAKE_DATE(year1x,month1x,day1x,starthour,0,0)
	IF hour1x MOD 6 NE 0 THEN datef = TIME_INC(datef, !TauDSolar/4.)  ;start from the first 6 hourly time after hour1x
	ENDELSE
    
    run_initDate=datename  ; 05/21/2013  add

    ;; create tempfile name with random string
    ran1x=year1x
    dk = outDir+'scratch'+sep
    IF FILE_TEST(dk) EQ 0 THEN FILE_MKDIR, dk
    fn = dk+'tempfile'+STRING(LONG(RANDOMU(ran1x,1)*100000)/2+$
        LONG(RANDOMU(seed,1)*100000)/2,format='(I05)')
        
    ; load methane oxidation parameters
    ; return PROD&LOSS[12(months),4(species),45(latitude),76(altitude)]
    ; lat=-88:4:88, pres[76]=[973.167,817.610,...,0.0022], zz[76]=[0.5,1.5,2.5,3.5,...,91 km]
    ; M(12,45,76): # of densities,
    ; C[12(months),4(species),45(latitude),76(altitude)]:constituent mixing ratios for 4 species
    RESTORE, ch4_oxid_para ;'pres' is only used for calculating dlnp, and then used when adding methane
    ch4loss=REFORM(loss(*,1,*,*)/c(*,1,*,*)) * 86400.   ; units of 1/day [month, lat, pres] = [12,45,76]
    ch4prod=REFORM(prod(*,1,*,*)) *1.e6 * 86400.       ; units of ppmv/day [month, lat, pres] = [12,45,76]
    lat2d=lat & dlnp=(ALOG(pres)-SHIFT(ALOG(pres),1))(1:*) & p2d=pres
    ; lat2d[45]:[-88, -84, ..., 88] p2d[76]: [943.167,817.610,708.768,...,0.00228971]
	IF b_MLS EQ 1 THEN BEGIN 
	ddi=datef.year*10000L+datef.month*100L+datef.day*1L
        IF b_6hrly EQ 1 THEN BEGIN
            ;; added by Hao Ye, 12/11/2014
            ;; read temperature and theta for T interpolation
                READ_VAR_UNIVERSAL, ddi, datef, temperature, 't', longrid, latgrid, zprgrid, dir=tempDir, hrFlag=1
                READ_VAR_UNIVERSAL, ddi, datef, thetax, 'theta', longrid, latgrid, zprgrid, dir=tempDir, hrFlag=1
        ENDIF

        prx = 1000.*(temperature/thetax)^(7./2.)

        prgrid = [3.831187e+02, 3.162278e+02, 2.610157e+02, 2.154435e+02, 1.778279e+02, $
       1.467799e+02, 1.211528e+02, 1.000000e+02, 8.254042e+01, 6.812920e+01,$
       5.623413e+01, 4.641589e+01, 3.831187e+01, 3.162278e+01, 2.610157e+01,$
       2.154435e+01, 1.778279e+01, 1.467799e+01, 1.211528e+01, 1.000000e+01]
        npr=N_ELEMENTS(prgrid)
        dpr=(prgrid-SHIFT(prgrid,1))(1:*)

        ;interpolate theta[540,359,24] according to the intendedThetaGrid[64]
        ; then according to T[540,359,24] get new temperature in finer theta grid: tnew[540,359,64]
        ; 3. get temperature on our intended theta levels (thgrid): 'tnew' is still of reanalysis grid
        thetanew=SF_ZTOTHETA(thetax,prx,prgrid)
        xi=(lon1x-longrid(0))/(longrid(1)-longrid(0))     ;[1350,33]  0~600?
        yi=(lat1x-latgrid(0))/(latgrid(1)-latgrid(0))     ;[1350,33]  50~330?
        s = SIZE(lon1x)
        i_pre_over= WHERE(prgrid GT pr1x-1)
        zi = i_pre_over(-1)*(FLTARR(s[1])+1)
		theta1x=INTERPOLATE(thetanew,xi,yi,zi) ;temf[] , MISSING=!VALUES.F_NAN
	ENDIF
 
    olon = lon1x
    olat = lat1x
    oth = theta1x

   
    lon=[0.] & lat=[0.] & th=[0.] & pr=[0.] &  h2o_mix=[0.] 


 
    jd=0L
    
    ; start loop ---------------------
    firstDay = datef.year*10000L+datef.month*100L+datef.day*1L
    preDay = firstDay
    first_enter = 1  ;first time enter, injecting parcels; also, remove the first elements of all variables
    
    WHILE 1 EQ 1 DO BEGIN
    
        PRINT, '########################### a new loop starts ######################'
        
        datei=datef ;datei: Initial date and time (being overwritten after each day)    datef: Final date and time
        
        IF b_daily EQ 1 THEN datef  = TIME_INC(datef, -!TauDSolar) ; only one day time step(dt: time increment in seconds (positive or negative))
        IF b_6hrly EQ 1 THEN datef  = TIME_INC(datef, -!TauDSolar/4.) ; only one day time step(dt: time increment in seconds (positive or negative))
        
        currDay=datei.year*10000L+datei.month*100L+datei.day*1L
        
        IF currDay NE preDay THEN BEGIN
            PRINT, '$$$ day accumulating, jd=jd+1 $$$  ', jd
                jd=jd+1 ; day accumulating
            convcnt=0 ; reset convective count at the beginning of each day
        ENDIF
        
        ch4_currYr = ch4_actl[datei.year]
        
        ; if using WACCM run, then read WACCM CH4 amount on daily basis
        ddi_str=STRING(currDay,format='(I08)')
        
    
        ;; injecting parcels into air: only inject "stopinjecting" years!
        ; first time entering, must inject parcels into the air!!!!
        ; old injecting scheme: inject every day for the first year, so (currDay NE preDay) and (stopinjecting GT 0)
        ; IF (first_enter EQ 1) OR (currDay NE preDay) THEN BEGIN
        IF (first_enter EQ 1) OR ((currDay NE preDay) AND ((datei.day MOD injecting_intvl) EQ 0)) THEN BEGIN
        
            IF (stopinjecting GT 0)  THEN BEGIN ;;OR (inject_once EQ 1)

                ;color those new parcels
 
                lon=[lon,olon] & lat=[lat,olat]
                th=[th,oth] ; parcel initiated at isentropic surfaces

                
                h2o_mix=[h2o_mix,olat*0+initWater]                        ;h2o           initiated with 200 ppmv

                
                PRINT, '*** step into a new day, injecting parcels, now have  ', N_ELEMENTS(lon), 'parcels ***'
                stopinjecting = 0 
                inject_once=0
                
            ENDIF ;stopinjecting > 0: only inject 3 years to the air
        ENDIF
        
        ; for the first time of running, remove the first elements of value 0 that initiated at each variables
        IF (first_enter EQ 1) AND jd EQ 0 THEN BEGIN ; will never enter this loop after first enter!!!!
            lon=lon(1:*) & lat=lat(1:*) & th=th(1:*) 
            h2o_mix=h2o_mix(1:*) 
            lonrecs=0 & latrecs=0 & threcs=0 & prrecs=0 &  h2o_mix_recs=h2o_mix

        ENDIF ;(first_enter EQ 1) AND jd EQ 0
        
        ; remove parcels that are too low or too high
        PRINT,STRCOMPRESS('@@@'+MAKE_ISO_DATE_STRING(datei)+': Total number of parcels index'+STRING(N_ELEMENTS(index)))+ '   lon'+STRING(N_ELEMENTS(lon))
     

        ;; define starting position for parcels
        ;; construct input thrcel structure
        x        = {name         : 'Longitude', values   : lon, units    : 'Degrees'}
        y        = {name         : 'Latitude',  values   : lat, units    : 'Degrees'}
        z        = {name         : 'Theta',     values   : REFORM(th, N_ELEMENTS(th)), units    : 'K'}
        my_seed=0
        
        ri={x: x, y: y, z: z, seed: my_seed} ;Initial particle positions; Return thrcel initial positions
        yearxx=datei.year
        
        yearxx=UINT(datei.year)*100L+UINT(datei.month)     ;monthly written wind
        
        ; this is the input wind file
        windFile=FILE_SEARCH(windDir+'*_'+STRING(yearxx, format='(I6)')+'.nc')  ;;  str routine
        
        ;yearxx=UINT(datei.year)*10000L+UINT(datei.month)*100L+UINT(datei.day)*1L      ; test purpose
        TRAJ3D_2_3, windFile, datei, datef, ntdays, ntsave, ntplot, ri, /CLOBBER, OUTFILE = fn+'.tmp', TW_DEBUG=tw_debug

        ;reads traj run into the structure outs (outs is a structure includes lonf, latf, thf, dayno)
        READ_TRAJ,fn+'.tmp',outs, lon=lonf,lat=latf,pres=thf,dayno=dayno ;here pres is not used
        ; return lonf, latf, thf, dayno    temperature is later interpolated from MERRA data into trj grid
        ; trajectory output is the position of parcels: lonf, latf, thf, dayno
        SPAWN,'rm -f '+fn+'.tmp'
        PRINT, 'ch4 current year:', ch4_currYr
   
        ddi=datei.year*10000L+datei.month*100L+datei.day*1L ;20050101
        
        PRINT,'date  ',ddi
        mon=datei.month
        year=datei.year
        
        PRINT,' theta ',MAX(thf, /NAN),' min theta ',MIN(thf, /NAN)
        
        ; $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        ; $$$$$$$$$$$$$$$$$$$$$ this part can be replaced by reading already interpolated T $$$$$$$$$$$$$$$$$$
        ; read reanalysis temperature and convert to theta: 'temperature' is of reanalysis grids
        ; MERRA return: tempearture(3D), longrid, latgrid, zprgrid (24 levels in heights from H(2hPa) to H(500hPa), top-->bottom
        ; read_temp_merra_day_tw, ddi, temperature,longrid,latgrid,zprgrid, dir=tempDir ;/sn1/dessler/merra/temperatures/
         
        IF b_6hrly EQ 1 THEN BEGIN
            ;; added by Hao Ye, 12/11/2014
            ;; read temperature and theta for T interpolation
                READ_VAR_UNIVERSAL, ddi, datei, temperature, 't', longrid, latgrid, zprgrid, dir=tempDir, hrFlag=1 
                READ_VAR_UNIVERSAL, ddi, datei, thetax, 'theta', longrid, latgrid, zprgrid, dir=tempDir, hrFlag=1
         

        ENDIF  
       
       
        ; now create the theta grid field to get the temperatures
        ; get temperature on our intended theta levels (thgrid): 'tnew' is still of reanalysis grid
        ; theta[64]: 330:10:360K, 362:2:440K,450:20:570K, 600:100:1800K
        thgrid=[[300.,310.,320.,330., 335., 340., 345., 350., 352.5, 355., 357.5, 360.], $ ;352.5,  357.5,
            FINDGEN(40)*2+362, $   ;362:2:440
            [450,470,490,510,530,550,570,600,700,800,900,1000,$
            1100,1200,1300,1400,1600,1800, 2000, 2400, th_top]] ; intended 64 levels

        nth=N_ELEMENTS(thgrid)
        dth=(thgrid-SHIFT(thgrid,1))(1:*)
        
        ;interpolate theta[540,359,24] according to the intendedThetaGrid[64]
        ; then according to T[540,359,24] get new temperature in finer theta grid: tnew[540,359,64]
        ; 3. get temperature on our intended theta levels (thgrid): 'tnew' is still of reanalysis grid
        tnew=SF_ZTOTHETA(temperature,thetax,thgrid) ;this has cubSpl interp when b_cub_spl is set.
        ;now we have temperature on desired theta grid (thgrid) -- ultimate goal !!!
        
       
        
        ; interpolate tnew to our model grids: 'temf' is new temperature at our model grids
        ;; interpolate temperature(desired theta grid) to trajectories: get temperature at parcel locations
        ; lonf: trajectory [1350, 33]   0~360°* ntdays    longrid:MERRA temperature [540] within 0~360°
        ; latf: trajectory [1350, 33]  -60°~60°* ntdays    latgrid:MERRA temperature [540] within -89.5°~89.5°
        ; thf: trajectory [1350, 33]  -60°~60°* ntdays    latgrid:MERRA temperature [540] within -89.5°~89.5°
        ; s: size info of trajectory results!!!!!
        s=SIZE(latf)
        xi=(lonf-longrid(0))/(longrid(1)-longrid(0))     ;[1350,33]  0~600?
        yi=(latf-latgrid(0))/(latgrid(1)-latgrid(0))     ;[1350,33]  50~330?
        ind=VALUE_LOCATE(thgrid,thf)                  ;thgrid[64](std theta levels)   thf[1350, 33](from trajectory)
        zi=ind+(thf-thgrid(ind))/dth(ind)                  ;[1350,33]   7.5-8.5?
        nlen=s(1)*s(2) ;total length: lonNum*latNum from trajectory, i.e., wind grid: 1350*33=[(x*y)*tStepNum]
        ; tnew[540,359,64] --> temf[nlen]  tnew in theta grid into trajectory grid
        
        temf=INTERPOLATE(tnew,REFORM(xi,nlen),REFORM(yi,nlen),REFORM(zi,nlen)) ;temf[] , MISSING=!VALUES.F_NAN
        
        temf=REFORM(temf,s(1),s(2))                               ;[1350, 33]
        ;n=WHERE(temf LT 100 OR temf GT 350,cnt)    ;unreasonable temperature: <-173 or >70
        ;IF cnt GT 0 THEN temf(n)=initTemp                      ;adjust unreasonable temperature to 300K

        prf=1000.*(temf/thf)^(7./2.)

        PRINT,STRCOMPRESS('OF ALL STEPS, '+'maxTEMPF:'+STRING(MAX(temf, /NAN))+', minTEMPF:'+STRING(MIN(temf, /NAN)) +$
            ', maxTHF:'+STRING(MAX(thf, /NAN))+', minTHF:'+STRING(MIN(thf, /NAN))+$
            ', maxPRF:'+STRING(MAX(prf, /NAN))+', minPRF:'+STRING(MIN(prf, /NAN)))

        ;we also need the last step information [1350x1]
        ; last step of lat - latl is used to calculating methane oxidation
        
        ; if is the firststep, only keep first (6-n mod 6 +1) elemtents
        
        IF firststep EQ  1 THEN BEGIN
         IF hour1x MOD 6 NE 0 THEN BEGIN
		  lonf = lonf[*,0:(hour1x mod 6)]
          latf = latf[*,0:(hour1x mod 6)]
          prf = prf[*,0:(hour1x mod 6)]
          thf = thf[*,0:(hour1x mod 6)]
          temf = temf[*,0:(hour1x mod 6)]
		ENDIF
          firststep = 0
        ENDIF
        lonl=lonf[*, -1] &  latl=latf[*, -1] &  prl=prf[*, -1] & thl=thf[*, -1] & teml=temf[*,-1]
  
        PRINT,STRCOMPRESS('For LAST STEP, '+'maxTEMPL:'+STRING(MAX(teml, /NAN))+', minTEMPL:'+STRING(MIN(teml, /NAN)) +$
            ', maxTHL:'+STRING(MAX(thl, /NAN))+', minTH:'+STRING(MIN(thl, /NAN))+$
            ', maxPRL:'+STRING(MAX(prl, /NAN))+', minPRL:'+STRING(MIN(prl, /NAN)))
 
        ; add gravity wave temperature perturbations to parcel temperatures
        s=SIZE(latf)	
        time_step=FINDGEN(s(2))/s(2) ; time steps in days
        s=SIZE(prrecs)
        IF s[0] EQ 0 THEN BEGIN
          lonrecs=lonf
          latrecs=latf
          threcs=thf
		  prrecs=prf
          temprecs=temf
        ENDIF ELSE BEGIN
		  lonrecs=[[lonrecs],[lonf(*,1:*)]]
		  latrecs=[[latrecs],[latf(*,1:*)]]
          threcs=[[threcs],[thf(*,1:*)]]
          prrecs=[[prrecs],[prf(*,1:*)]]
          temprecs=[[temprecs],[temf(*,1:*)]]
        ENDELSE

        ; now get the water field - this has to be done with full time resolution temperature
        ; time step through the trajectories
  ;      FOR jh=0,s(2)-1 DO BEGIN ;s(2)=33=ntdays+1 times recording of trajectory results
  ;          lathr=latf(*, jh)     ;each time of trajectory lat; total ntdays+1 times
  ;          lonhr=lonf(*, jh)     ;each time of trajectory lon; total ntdays+1 times
  ;          thhr=thf(*, jh)      ;each time of trajectory theta ; total ntdays+1 times
  ;          temhr=temf(*, jh) ;each time of trajectory temperature(interpolated from MERRA T)
  ;          prhr=prf(*, jh)        ; prhr can be converted from zphr, or originally from temhr and thhr pr=1000*exp(-zp/7)
            ; 07/09/2013, change, let it suitable for subsaturated case!
  ;          h2olimit=SAT_MIX_MK(prhr, temhr)*supersat ;<initwater      ;from pr, temp --> thearetical saturation mixing ratio
          
            ; ######### non convective influenced parcels here ##############
            
            ; first water parcels without methane
            ;n=WHERE(h2o_mix/h2olimit GT supersat,cnt)
 ;           n=WHERE(h2o_mix/h2olimit GT 1.0, cnt) ;07/09/2013, change, let it suitable for subsaturated case!
 ;           temp_h2o_mix=h2o_mix;!!!!!!!!!delete after use    ;2004
 ;           IF cnt GT 0 THEN BEGIN    

 ;               ;07/09/2013, change, let it suitable for subsaturated case!
 ;               IF supersat GT 1.0 THEN BEGIN
 ;                   h2o_mix(n)=(SAT_MIX_MK(prhr, temhr))[n] ;for supersat>100%, always dehydrate to 100% RH
 ;               ENDIF ELSE BEGIN
 ;                   h2o_mix(n)=h2olimit(n) ;for supersat <=100%, always dehydrate to current supersat level

;                ENDELSE
                
              
               
;            ENDIF
            

   
;        ENDFOR ; jh steps loop end

        
        ;last position/temperatures of parcels
        lon=TEMPORARY(lonl)
        lat=TEMPORARY(latl)
        th=TEMPORARY(thl)
        tem=TEMPORARY(teml)
        pr=TEMPORARY(prl)
    
        n=WHERE(tem EQ 0,cnt)
        IF cnt GT 0 THEN BEGIN
            PRINT, 'Temperature error - continuing,',cnt
            tem(n)=initTemp
        ENDIF
        

        ; save data and restart file section --------------------
		s = SIZE(lonrecs)
        IF s[2] GE 24*recday  THEN BEGIN;; AND (datei.year GE 2004L)

            PRINT,'$$$$$$$$$$$ ', MAKE_ISO_DATE_STRING(datef), '   writting data into files $$$$'
                strout=STRCOMPRESS(STRING(datename.year,format='(I04)')+'_'+ $
                STRING(datename.Day,format='(I02)')+STRING(datename.hour,format='(I02)')+ $
             '_I'+ STRING(findex,format='(I06)'), /REMOVE_ALL)+'_E'+STRING(poindex,format='(I01)')
            datestr=STRCOMPRESS(STRING(datei.year)+STRING(datei.month,format='(I02)')+ $
                STRING(datei.DAY,format='(I02)'), /REMOVE_ALL)
            file_names='traj_'+satname+'_'+strout+'.sav'
            
			lonrecs = lonrecs(*,0:24*recday-1)
            latrecs = latrecs(*,0:24*recday-1)
            prrecs = prrecs(*,0:24*recday-1)
            temprecs = temprecs(*,0:24*recday-1)
		 
		FOR i_file = 0,s[1]-1 DO BEGIN 
			file_name = file_names[i_file] 
			lonrec = REFORM(lonrecs[i_file,*])
			latrec = REFORM(latrecs[i_file,*])
			threc = REFORM(threcs[i_file,*])
			prrec = REFORM(prrecs[i_file,*])
			temprec = REFORM(temprecs[i_file,*])
			h2o_mix_rec = REFORM(h2o_mix_recs[i_file,*])
			PRINT,'saving file '+outDir+file_name
            ; add header in file
            
			header={run_initDate:run_initDate, $
                b_daily:b_daily, b_6hrly:b_6hrly, year1x:year1x, $
                supersat:supersat, $
                windDir:windDir, windFile:windFile, tempDir:tempDir, $               
                $
                savFolder:savFolder, outDir:outDir}
		SAVE,FILENAME=outDir+file_name, jd,datestr,header,datei,datef, ddi, $
                lonrec, latrec, threc, prrec, temprec, h2o_mix_rec
    ENDFOR 
            PRINT, 'latRange: ', MIN(lat, /NAN), MAX(lat, /NAN)
            HELP, lon
            GOTO,skipp ;
            
        ENDIF
        
        PRINT, '*** ', MAKE_ISO_DATE_STRING(datei), ' to ', MAKE_ISO_DATE_STRING(datef), ', finished ***'
        
        PRINT
        
        IF (first_enter EQ 1) THEN BEGIN
            first_enter = 0
            PRINT, 'This is the first enter, and it finished!'
        ENDIF
        
        
        IF currDay NE preDay THEN BEGIN
            PRINT, 'chaning day: ', 'preday-', preDay, 'currDay-', currDay
            preDay = currDay
        ENDIF
        
    ENDWHILE
    
    skipp:
    PRINT,'program terminated normally'
    
    SPAWN,'rm -f '+fn+'.tmp'
    
END
