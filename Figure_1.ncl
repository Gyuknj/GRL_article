begin

;======================================
  yrStrt = 1982
  yrLast = 2022
  
  latS   = -60.
  latN   =  60.
  lonL   = -180.
  lonR   =  180.
  
   nmos = 12	;months in a year

  mxlag = 7
  lag = 1	;first-order autocorrelation
  
  thres = 1.0	;set a threshold
  
  xlabel = ispan(0,nmos-1,1)
  ;print(xlabel)

  p90 = 0.10
  
  xy_colors_dheating        = (/"red","deepskyblue","darksalmon","black","deeppink4","yellow2","blue","seagreen3"/)
  xy_colors_dheating1        = (/"deepskyblue","deepskyblue","darksalmon","black","deeppink4","yellow2","blue","seagreen3"/)
  xy_colors_dheating2        = (/"blue","deepskyblue","darksalmon","black","deeppink4","yellow2","blue","seagreen3"/)
  
;======================================================================================================================================
  f      = addfile("/home/gyuk/Downloads/GYUK/GEOMAR/Final_work/HadISST_sst.nc", "r")
  TIME   = f->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  hadi   = f->sst(iYYYY,:,:)
  hadi!1 = "lat"
  hadi!2 = "lon"

  delete(f)		;	Not needed anymore
  delete(TIME)
  delete(YYYY)
  delete(iYYYY)

  hadi = where(hadi.gt.-30.,hadi,hadi@_FillValue)	;removes padded values
  ;printVarSummary(hadi)       ; variable overview
 ; printMinMax(hadi,True)
;===================================================================================================================================
  hadiClm      = clmMonTLL(hadi)
  ;printVarSummary(hadiClm)                   ; (12,nlat,nlon)   
  
  hadi         = calcMonAnomTLL(hadi,hadiClm)    ; replace with anonamlies
  hadi@long_name = "ANOMALIES: "+hadi@long_name
  delete(hadiClm)
;===================================================================================================================================
  nyears = dimsizes(hadi&time)/nmos
  ;print(nyears)

  hadi_s     = onedtond(ndtooned(hadi),(/nyears,nmos,dimsizes(hadi&lat),dimsizes(hadi&lon)/))
  copy_VarCoords(hadi(0,:,:),hadi_s(0,0,:,:))                             ; copy coordinate arrays
  hadi_s!0 = "time"
  hadi_s&time = fspan(yrStrt,yrLast,nyears)
  hadi_s!1 = "month"
  ;printVarSummary(hadi_s)

  hadi_s = dtrend_msg_n(hadi_s&time,hadi_s,True,False,0)        ;removes linear trends
  ;printVarSummary(hadi_s)
  ;printMinMax(hadi_s,True)
;---------------------------------------------------------------------------------------------------------------------------------------------------
  ahadi_s = wgt_areaave_Wrap(hadi_s(time|:,month|:,{lon|-20:0}, {lat|-3:3}),1.0, 1.0, 0)        ;atl3 index
  printVarSummary(ahadi_s)
  ;printMinMax(ahadi_s,True)
;------------------------------------------------------------------------------------------------------------------------------------
  delete(hadi)
 ; delete(hadi_s)
;------------------------------------------------------------------------------------------------------------------------------------
  jja_ahadi = dim_avg_n_Wrap(ahadi_s(time|:,month|5:7),1)		;JJA, ncl counts from zero
 ; print(jja_ahadi)
;------------------------------------------------------------------------------------------------------------------------------------
  jja_ahadi = dim_standardize_n(jja_ahadi, 0, 0)  ;Normalise deviations
 ; print(jja_ahadi)
;----------------------------------------------------------------------------
; Positive cases
;============================================================================
  n_atlP = dim_num_n(jja_ahadi.ge.thres,0)        ;Number of cases
 ; print(n_atlP)

  valsP = ind(jja_ahadi.ge.thres)
 ; print(valsP)
 ; print(valsP+yrStrt)		;the years
;------------------------------------------------------------------------------------------------------------------------------------
  atlP_sst = ahadi_s(valsP,:)
 ; print(atlP_sst)
;----------------------------------------------------------------------------
; Negative cases,
; note the "-" sign attached to the threshold here
;============================================================================
  n_atlN = dim_num_n(jja_ahadi.le.-thres,0)        ;Number of cases
 ; print(n_atlN)
;------------------------------------------------------------------------------------------------------------------------------------
  valsN = ind(jja_ahadi.le.-thres)
  ;print(valsN)
  ;print(valsN+yrStrt)             ;the years
;------------------------------------------------------------------------------------------------------------------------------------
  atlN_sst = ahadi_s(valsN,:)
 ; print(atlN_sst)
 delete(jja_ahadi)
;====================================================================================================================================
  f      = addfile("/home/gyuk/Downloads/GYUK/GEOMAR/Final_work/precip.mon.mean.nc", "r")
  TIME   = f->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  rain   = f->precip(iYYYY,:,:)
  printVarSummary(rain)       ; variable overview
  delete(f)             ;       Not needed anymore
  delete(TIME)
  delete(YYYY)
  delete(iYYYY)

  rain!1 = "lat"
  rain!2 = "lon"
;======================================================================================================================================
  rainClm      = clmMonTLL(rain)
  printVarSummary(rainClm)                   ; (12,nlat,nlon)   
  rain         = calcMonAnomTLL(rain,rainClm)    ; replace with anonamlies
  rain@long_name = "ANOMALIES: "+rain@long_name
  delete(rainClm)
;----------------------------------------------------------------------------------------------------------------------------------------
  nyears = dimsizes(rain&time)/nmos
  ;print(nyears)
  rain_s     = onedtond(ndtooned(rain),(/nyears,nmos,dimsizes(rain&lat),dimsizes(rain&lon)/))
  copy_VarCoords(rain(0,:,:),rain_s(0,0,:,:))                             ; copy coordinate arrays
  rain_s!0 = "time"
  rain_s&time = fspan(yrStrt,yrLast,nyears)
  rain_s!1 = "month"
  ;printVarSummary(rain_s)

  rain_s = dtrend_msg_n(rain_s&time,rain_s,True,False,0)        ;removes linear trends
  ;printVarSummary(rain_s)
  ;printMinMax(rain_s,True)
;--------------------------------------------------------------------------------------------------------------------------------------- (3°N–3°S, 5°E–40°W)
  arain_s = wgt_areaave_Wrap(rain_s(time|:,month|:,{lon|-40:5}, {lat|-3:3}),1.0, 1.0, 0)        ;Nnamchi et al., 2021 definition
  ;printVarSummary(arain_s)
  ;printMinMax(arain_s,True)
  
  arain_sA = wgt_areaave_Wrap(rain_s(time|:,month|:,{lon|-15:10}, {lat|4:10}),1.0, 1.0, 0)        ;Nnamchi et al., 2021 definition
  ;printVarSummary(arain_sA)
  ;printMinMax(arain_sA,True)
;---------------------------------------------------------------------------------------------------------------------------------------4-10N 15W-10E
  delete(rain)
 ; delete(rain_s)
;-------------------------------------------------------------------------------
; Take out the time series
;===============================================================================
  atlP_rain = arain_s(valsP,:)
 ; print(atlP_rain)
;-------------------------------------------------------------------------------
; Take out the time series
;===============================================================================
  atlN_rain = arain_s(valsN,:)
 ; print(atlN_rain)
 
; Take out the time series
;===============================================================================
  atlP_rainA = arain_sA(valsP,:)
 ; print(atlP_rainA)
;-------------------------------------------------------------------------------
; Take out the time series
;===============================================================================
  atlN_rainA = arain_sA(valsN,:)
 ; print(atlN_rainA) 
;***************************************************************************************************************
;***************************************************************************************************************
;  Positive  Canonical events SSTA 
;***************************************************************************************************************  
  atlP_sst1 = (atlP_sst(time|2,month|:))
  printVarSummary(atlP_sst1)
  atlP_sst2 = (atlP_sst(time|3,month|:))
  printVarSummary(atlP_sst2)
  atlP_sst3 = (atlP_sst(time|4,month|:))
  printVarSummary(atlP_sst3)
  atlP_sst4 = (atlP_sst(time|6,month|:))
  printVarSummary(atlP_sst4)
  atlP_sst5 = (atlP_sst(time|7,month|:))
  printVarSummary(atlP_sst5)

  canonical_eventsP = (/atlP_sst1,atlP_sst2,atlP_sst3,atlP_sst4,atlP_sst5/)       
  copy_VarCoords_n(atlP_sst,canonical_eventsP,(/0,1/))
  canonical_eventsP!0 = "nyears"    
  printVarSummary(canonical_eventsP) 
 
  delete([/atlP_sst1, atlP_sst2, atlP_sst3, atlP_sst4, atlP_sst5/])
;************************************************************************************************************
;  Positive  Canonical events precipitation
;***************************************************************************************************************   
  atlP_rain1 = (atlP_rain(time|2,month|:))
  atlP_rain2 = (atlP_rain(time|3,month|:))
  atlP_rain3 = (atlP_rain(time|4,month|:))
  atlP_rain4 = (atlP_rain(time|6,month|:))
  atlP_rain5 = (atlP_rain(time|7,month|:))
  
  canonical_events_rainP = (/atlP_rain1,atlP_rain2,atlP_rain3,atlP_rain4,atlP_rain5/)       
  copy_VarCoords_n(atlP_rain,canonical_events_rainP,(/0,1/))
  canonical_events_rainP!0 = "nyears"
  printVarSummary(canonical_events_rainP) 
  
  delete([/atlP_rain1, atlP_rain2, atlP_rain3, atlP_rain4, atlP_rain5/])   
;***************************************************************************************************************
;  Negative  Canonical events SSTA 
;***************************************************************************************************************  
  atlN_sst1 = (atlN_sst(time|0,month|:))
  atlN_sst2 = (atlN_sst(time|1,month|:))
  atlN_sst3 = (atlN_sst(time|2,month|:))
  atlN_sst4 = (atlN_sst(time|5,month|:))

  canonical_eventsN = (/atlN_sst1,atlN_sst2,atlN_sst3,atlN_sst4/)       
  copy_VarCoords_n(atlP_sst,canonical_eventsN,(/0,1/))
  canonical_eventsN!0 = "nyears"    
  printVarSummary(canonical_eventsN) 
;===================================== 
  canonical_eventsN_s = dim_avg_n_Wrap(canonical_eventsN,0)
  printVarSummary(canonical_eventsN_s) 
  
  delete([/atlN_sst1, atlN_sst2, atlN_sst3, atlN_sst4/])
;-------------------------------------------------------------------------------
;************************************************************************************************************
;  Negative  Canonical events precipitation
;***************************************************************************************************************   
  atlN_rain1 = (atlN_rain(time|0,month|:))
  atlN_rain2 = (atlN_rain(time|1,month|:))
  atlN_rain3 = (atlN_rain(time|2,month|:))
  atlN_rain4 = (atlN_rain(time|5,month|:))

  canonical_events_rainN = (/atlN_rain1,atlN_rain2,atlN_rain3,atlN_rain4/)       
  copy_VarCoords_n(atlP_rain,canonical_events_rainN,(/0,1/))
  canonical_events_rainN!0 = "nyears"    
  printVarSummary(canonical_events_rainN) 
  
   delete([/atlN_rain1, atlN_rain2, atlN_rain3, atlN_rain4/])
;***************************************************************************************************************
;;;; Difference in Canonical events
;*************************************************************************************************************** 
  canonical_events_difference = (dim_avg_n_Wrap(canonical_eventsP,0) - dim_avg_n_Wrap(canonical_eventsN,0))/(2.)
  copy_VarCoords(canonical_eventsN_s,canonical_events_difference)                             ; copy coordinate arrays 
  printVarSummary(canonical_events_difference)
  printMinMax(canonical_events_difference,True)
;*************************************************************************************************************** 
;*************************************************************************************************************** 
  pc1SSTd = canonical_events_difference  ; Transfer metdata for plotting sig. values        

  pc1SSTd = ttest((dim_avg_n_Wrap(canonical_eventsP,0)),(dim_stddev_n_Wrap(canonical_eventsP,0))^2,5\
               ,(dim_avg_n_Wrap(canonical_eventsN,0)),(dim_stddev_n_Wrap(canonical_eventsN,0))^2,4,True, False)
  pc1SSTd@long_name = "Probability: difference between means"
  printVarSummary(pc1SSTd)

  pc1SSTd = where(abs(pc1SSTd).le.p90,canonical_events_difference,pc1SSTd@_FillValue) ; Sig.u winds, where UWND or VWND is significanlt
  printVarSummary(pc1SSTd)
;*************************************************************************************************************** 
  canonical_events__rain_difference = (dim_avg_n_Wrap(canonical_events_rainP,0) - dim_avg_n_Wrap(canonical_events_rainN,0))/(2.)
  copy_VarCoords(canonical_eventsN_s,canonical_events__rain_difference)                             ; copy coordinate arrays 
  printVarSummary(canonical_events__rain_difference)
  printMinMax(canonical_events__rain_difference,True)
  
  canonical_events__rain_difference2 = dim_avg_n_Wrap(canonical_events__rain_difference,0)
;*************************************************************************************************************** 
  pc1raind = canonical_events__rain_difference  ; Transfer metdata for plotting sig. values        
  pc1raind = ttest((dim_avg_n_Wrap(canonical_events_rainP,0)),(dim_stddev_n_Wrap(canonical_events_rainP,0))^2,5\
               ,(dim_avg_n_Wrap(canonical_events_rainN,0)),(dim_stddev_n_Wrap(canonical_events_rainN,0))^2,4,True, False)
  pc1raind@long_name = "Probability: difference between means"
  printVarSummary(pc1raind)

  pc1raind = where(abs(pc1raind).le.p90,canonical_events__rain_difference,pc1raind@_FillValue)
  printVarSummary(pc1raind)
  
  delete([/canonical_events_rainP, canonical_events_rainN, canonical_eventsP, canonical_eventsN/])
;***************************************************************************************************************
;  Positive  Non-Canonical events SSTA 
;***************************************************************************************************************  
  atlP_sst6 = (atlP_sst(time|0,month|:))
  atlP_sst7 = (atlP_sst(time|1,month|:))
  atlP_sst8 = (atlP_sst(time|5,month|:))
  atlP_sst9 = (atlP_sst(time|8,month|:))
  
  canonical_eventsP2 = (/atlP_sst6,atlP_sst7,atlP_sst8,atlP_sst9/)       
  copy_VarCoords_n(atlP_sst,canonical_eventsP2,(/0,1/))
  canonical_eventsP2!0 = "nyears"   
 ; printVarSummary(canonical_eventsP2)  
  
  delete([/atlP_sst6, atlP_sst7, atlP_sst8, atlP_sst9/])
;************************************************************************************************************
;  Positive  Non-Canonical events precipitation
;***************************************************************************************************************   
  atlP_rain6 = (atlP_rain(time|0,month|:))
  atlP_rain7 = (atlP_rain(time|1,month|:))
  atlP_rain8 = (atlP_rain(time|5,month|:))
  atlP_rain9 = (atlP_rain(time|8,month|:))
 
  canonical_events_rainH = (/atlP_rain6,atlP_rain7,atlP_rain8,atlP_rain9/)       
  copy_VarCoords_n(atlP_rain,canonical_events_rainH,(/0,1/))
  canonical_events_rainH!0 = "nyears"    
 ; printVarSummary(canonical_events_rainH) 
  
  delete([/atlP_rain6, atlP_rain7, atlP_rain8, atlP_rain9/])
;***************************************************************************************************************
;  Negative  Non-Canonical events SSTA 
;***************************************************************************************************************  
  atlN_sst6 = (atlN_sst(time|3,month|:))
  atlN_sst7 = (atlN_sst(time|4,month|:))
  atlN_sst8 = (atlN_sst(time|6,month|:))
  atlN_sst9 = (atlN_sst(time|7,month|:)) 
 
  canonical_eventsN2 = (/atlN_sst6,atlN_sst7,atlN_sst8,atlN_sst9/)       
  copy_VarCoords_n(atlP_sst,canonical_eventsN2,(/0,1/))
  canonical_eventsN2!0 = "lev"    
 ; printVarSummary(canonical_eventsN2) 
  
  delete([/atlN_sst6, atlN_sst7, atlN_sst8, atlN_sst9/]) 
;************************************************************************************************************
;  Negative  Non-Canonical events precipitation
;***************************************************************************************************************   
  atlN_rain6 = (atlN_rain(time|3,month|:))
  atlN_rain7 = (atlN_rain(time|4,month|:))
  atlN_rain8 = (atlN_rain(time|6,month|:))
  atlN_rain9 = (atlN_rain(time|7,month|:))
  
  canonical_events_rainK = (/atlN_rain6,atlN_rain7,atlN_rain8,atlN_rain9/)       
  copy_VarCoords_n(atlP_rain,canonical_events_rainK,(/0,1/))
  canonical_events_rainK!0 = "lev"    
 ; printVarSummary(canonical_events_rainK)  
  
  delete([/atlN_rain6, atlN_rain7, atlN_rain8, atlN_rain9/])
;-------------------------------------------------------------------------------
;;;; Difference in Non-Canonical events
;*************************************************************************************************************** 
;=======================SSTA==============================================================
  canonical_events_differenceNON = (dim_avg_n_Wrap(canonical_eventsP2,0) - dim_avg_n_Wrap(canonical_eventsN2,0))/(2.)
  copy_VarCoords(canonical_eventsN_s,canonical_events_differenceNON)                             ; copy coordinate arrays 
 ; printVarSummary(canonical_events_differenceNON)
  printMinMax(canonical_events_differenceNON,True)  
  canonical_events_differenceNON2 = dim_avg_n_Wrap(canonical_events_differenceNON,0)
;*************************************************************************************************************** 
  pc1SSTdNON = canonical_events_differenceNON  ; Transfer metdata for plotting sig. values        

  pc1SSTdNON = ttest((dim_avg_n_Wrap(canonical_eventsP2,0)),(dim_stddev_n_Wrap(canonical_eventsP2,0))^2,4\
               ,(dim_avg_n_Wrap(canonical_eventsN2,0)),(dim_stddev_n_Wrap(canonical_eventsN2,0))^2,4,True, False)
  pc1SSTdNON@long_name = "Probability: difference between means"
  printVarSummary(pc1SSTdNON)

  pc1SSTdNON = where(abs(pc1SSTdNON).le.p90,canonical_events_differenceNON,pc1SSTdNON@_FillValue) 
 ; printVarSummary(pc1SSTdNON)
;=========================================================================================================================
;=======================Precipitation==============================================================
  canonical_events__rain_differenceJ = (dim_avg_n_Wrap(canonical_events_rainH,0) - dim_avg_n_Wrap(canonical_events_rainK,0))/(2.)
  copy_VarCoords(canonical_eventsN_s,canonical_events__rain_differenceJ)                             ; copy coordinate arrays 
  printVarSummary(canonical_events__rain_differenceJ)
  
  canonical_events__rain_differenceJ1 = dim_avg_n_Wrap(canonical_events__rain_differenceJ,0)
;*************************************************************************************************************** 
;===============================================================================================================
  pc1raind3 = canonical_events__rain_differenceJ  ; Transfer metdata for plotting sig. values        

  pc1raind3 = ttest((dim_avg_n_Wrap(canonical_events_rainH,0)),(dim_stddev_n_Wrap(canonical_events_rainH,0))^2,4\
               ,(dim_avg_n_Wrap(canonical_events_rainK,0)),(dim_stddev_n_Wrap(canonical_events_rainK,0))^2,4,True, False)
  pc1raind3@long_name = "Probability: difference between means"
  printVarSummary(pc1raind3)

  pc1raind3 = where(abs(pc1raind3).le.p90,canonical_events__rain_differenceJ,pc1raind3@_FillValue) 
  printVarSummary(pc1raind3)
  
  delete([/canonical_events_rainH, canonical_events_rainK, canonical_eventsP2, canonical_eventsN2/]) 
;===============================================================
;==============================================================
;=================== Canonical SSTA =================================
;==============================================================
  dimZ = dimsizes(hadi_s)
  printVarSummary(dimZ)
      
        saodP_cases_maps  = new((/n_atlP,dimZ(1),dimZ(2),dimZ(3)/),float,hadi_s@_FillValue)
        printVarSummary(saodP_cases_maps)
    
        saodP_cases_maps = hadi_s(valsP,:,:,:)
        printVarSummary(saodP_cases_maps)
       
        sstb1p = (saodP_cases_maps(time|2,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(sstb1p)    
  
        sstb2p = (saodP_cases_maps(time|3,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(sstb2p)    
  
        sstb3p = (saodP_cases_maps(time|4,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(sstb3p)    
  
        sstb4p = (saodP_cases_maps(time|6,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(sstb4p) 
  
        sstb5p = (saodP_cases_maps(time|7,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(sstb5p) 
     
       sst_fp = (/sstb1p,sstb2p,sstb3p,sstb4p,sstb5p/)       
       copy_VarCoords_n(hadi_s,sst_fp,(/1,2,3/))
       sst_fp!0 = "time"    
       printVarSummary(sst_fp)   
;==============================================================
;==============================================================    
        saodN_cases_maps  = new((/n_atlN,dimZ(1),dimZ(2),dimZ(3)/),float,hadi_s@_FillValue)
        printVarSummary(saodN_cases_maps)
    
        saodN_cases_maps = hadi_s(valsN,:,:,:)
        printVarSummary(saodN_cases_maps)
        
        sstb1n = (saodN_cases_maps(time|0,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(sstb1n)    
  
        sstb2n = (saodN_cases_maps(time|1,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(sstb2n)    
  
        sstb3n = (saodN_cases_maps(time|2,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(sstb3n)    
  
        sstb4n = (saodN_cases_maps(time|5,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(sstb4n) 
     
       sst_fn = (/sstb1n,sstb2n,sstb3n,sstb4n/)        
       copy_VarCoords_n(hadi_s,sst_fn,(/1,2,3/))
       sst_fn!0 = "time"   
       printVarSummary(sst_fn) 
;==============================================================
  pc1SST_C = (dim_avg_n_Wrap(sst_fp,0) - dim_avg_n_Wrap(sst_fn,0))/(2.)
  copy_VarCoords(hadi_s(0,:,:,:),pc1SST_C)                             ; copy coordinate arrays 
  printVarSummary(pc1SST_C)
;=====================================================================================================
  pc1SSTd_C = pc1SST_C  ; Transfer metdata for plotting sig. values        

  pc1SSTd_C = ttest((dim_avg_n_Wrap(sst_fp,0)),(dim_stddev_n_Wrap(sst_fp,0))^2,5\
               ,(dim_avg_n_Wrap(sst_fn,0)),(dim_stddev_n_Wrap(sst_fn,0))^2,4,True, False)
  pc1SSTd_C@long_name = "Probability: difference between means"
  printVarSummary(pc1SSTd_C)

  pc1SSTd_C = where(abs(pc1SSTd_C).le.p90,pc1SSTd_C,pc1SSTd_C@_FillValue)
  printVarSummary(pc1SSTd_C)
 
  SST_C = dim_avg_n_Wrap(pc1SST_C(month|5:7,lat|:, lon|:),0)            ; SSTA difference in events
  
  SSTdP_C = dim_avg_n_Wrap(pc1SSTd_C(month|5:7,lat|:, lon|:),0) 
;*************************************************************************************************************** 
  delete([/sst_fp, sst_fn/]) 
;*************************************************************************************************************** 
  sst_new = (/sstb1n,sstb2n,sstb3n,sstb4n,sstb1p,sstb2p,sstb3p,sstb4p,sstb5p/) 
  copy_VarCoords_n(saodP_cases_maps,sst_new,(/1,2,3/)) 
  sst_new!0 = "time"         
  printVarSummary(sst_new)
 
  delete([/sstb1n, sstb2n, sstb3n, sstb4n/])
  delete([/sstb1p, sstb2p, sstb3p, sstb4p, sstb5p/]) 
;*************************************************************************************************************** 
;*************************************************************************************************************** 
;==============================================================================================================
  delete(dimZ)
;=======================================================================================================================
;=================== Non-Canonical SSTA ======================================================  
;===============================================================================================================       
        sstb1p1 = (saodP_cases_maps(time|0,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(sstb1p1)    
  
        sstb2p1 = (saodP_cases_maps(time|1,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(sstb2p1)    
  
        sstb3p1 = (saodP_cases_maps(time|5,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(sstb3p1)    
  
        sstb4p1 = (saodP_cases_maps(time|8,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(sstb4p1) 
        
       sst_fp1 = (/sstb1p1,sstb2p1,sstb3p1,sstb4p1/)        
       copy_VarCoords_n(hadi_s,sst_fp1,(/1,2,3/))
       sst_fp1!0 = "time"  
       printVarSummary(sst_fp1) 
;================================================================================================================
        sstb1n1 = (saodN_cases_maps(time|3,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(sstb1n1)    
  
        sstb2n1 = (saodN_cases_maps(time|4,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(sstb2n1)    
  
        sstb3n1 = (saodN_cases_maps(time|6,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(sstb3n1)    
  
        sstb4n1 = (saodN_cases_maps(time|7,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(sstb4n1) 
       
       sst_fn1 = (/sstb1n1,sstb2n1,sstb3n1,sstb4n1/)       
       copy_VarCoords_n(hadi_s,sst_fn1,(/1,2,3/))
       sst_fn1!0 = "time"    
       printVarSummary(sst_fn1) 
;================================================================================================================
  pc1SST_N = (dim_avg_n_Wrap(sst_fp1,0) - dim_avg_n_Wrap(sst_fn1,0))/(2.)
  copy_VarCoords(hadi_s(0,:,:,:),pc1SST_N)                             ; copy coordinate arrays 
  printVarSummary(pc1SST_N)
;================================================================================================================ 
  pc1SSTd_N = pc1SST_N  ; Transfer metdata for plotting sig. values        

  pc1SSTd_N = ttest((dim_avg_n_Wrap(sst_fp1,0)),(dim_stddev_n_Wrap(sst_fp1,0))^2,4\
               ,(dim_avg_n_Wrap(sst_fn1,0)),(dim_stddev_n_Wrap(sst_fn1,0))^2,4,True, False)
  pc1SSTd_N@long_name = "Probability: difference between means"
  printVarSummary(pc1SSTd_N)

  pc1SSTd_N = where(abs(pc1SSTd_N).le.p90,pc1SSTd_N,pc1SSTd_N@_FillValue)
  ;print("==============================================================")  
  printVarSummary(pc1SSTd_N)

  SST_N = dim_avg_n_Wrap(pc1SST_N(month|5:7,lat|:, lon|:),0)

  SSTdp_N = dim_avg_n_Wrap(pc1SSTd_N(month|5:7,lat|:, lon|:),0)
;***************************************************************************************************************   
;*************************************************************************************************************** 
  sst_new_non = (/sstb1p1,sstb2p1,sstb3p1,sstb4p1,sstb1n1,sstb2n1,sstb3n1,sstb4n1/) 
  copy_VarCoords_n(saodP_cases_maps,sst_new_non,(/1,2,3/)) 
  sst_new_non!0 = "time"         
  printVarSummary(sst_new_non)
  delete([/sstb1p1, sstb2p1, sstb3p1, sstb4p1/])
  delete([/sstb1n1, sstb2n1, sstb3n1, sstb4n1/])

  ahadi_non = wgt_areaave_Wrap(sst_new_non(time|:,month|5,{lon|-20:0}, {lat|-3:3}),1.0, 1.0, 0)        ;atl3 index
  printVarSummary(ahadi_non)
  printMinMax(ahadi_non,True)
;*************************************************************************************************************** 
;*************************************************************************************************************** 
  x_can_non = new((/2,nmos/),float,canonical_events_difference@_FillValue)
  x_can_non(0:1,:) = (/canonical_events_difference,canonical_events_differenceNON/)
  printVarSummary(x_can_non)

  x_can_non!0 = "lev"
  x_can_non!1 = "month"
  printVarSummary(x_can_non)
;***************************************************************************************************************
  x_can_non_sta = new((/2,nmos/),float,canonical_events_difference@_FillValue)
  x_can_non_sta(0:1,:) = (/pc1SSTd,pc1SSTdNON/)
  printVarSummary(x_can_non_sta)
 
  x_can_non_sta!0 = "lev"
  x_can_non_sta!1 = "month"
  printVarSummary(x_can_non_sta)
;==========================================================================================================================
;***************************************************************************************************************
  rain_can_non = new((/2,nmos/),float,canonical_events__rain_difference@_FillValue)
  rain_can_non(0:1,:) = (/canonical_events__rain_difference,canonical_events__rain_differenceJ/)
  printVarSummary(rain_can_non)

  rain_can_non!0 = "lev"
  rain_can_non!1 = "month"
  printVarSummary(rain_can_non)
;***************************************************************************************************************
  rain_can_non_sta = new((/2,nmos/),float,pc1raind3@_FillValue)
  rain_can_non_sta(0:1,:) = (/pc1raind,pc1raind3/)
  printVarSummary(rain_can_non_sta)
 
  rain_can_non_sta!0 = "lev"
  rain_can_non_sta!1 = "month"
  printVarSummary(rain_can_non_sta) 
;=================================================================================================================================
  delete(nyears)
  delete([/pc1raind,pc1raind3, canonical_events__rain_difference,canonical_events__rain_differenceJ,pc1SSTd,pc1SSTdNON,canonical_events_difference,canonical_events_differenceNON,hadi_s/])
;==============================================================
;=================================================================================================================================
;  Positive  Canonical events precipitation
;***************************************************************************************************************   
  atlP_rain1A = (atlP_rainA(time|2,month|:))
  atlP_rain2A = (atlP_rainA(time|3,month|:))
  atlP_rain3A = (atlP_rainA(time|4,month|:))
  atlP_rain4A = (atlP_rainA(time|6,month|:))
  atlP_rain5A = (atlP_rainA(time|7,month|:))
  
  canonical_events_rainPA = (/atlP_rain1A,atlP_rain2A,atlP_rain3A,atlP_rain4A,atlP_rain5A/)       
  copy_VarCoords_n(atlP_rainA,canonical_events_rainPA,(/0,1/))
  canonical_events_rainPA!0 = "nyears"
  printVarSummary(canonical_events_rainPA) 
  
  delete([/atlP_rain1A, atlP_rain2A, atlP_rain3A, atlP_rain4A, atlP_rain5A/])   
;***************************************************************************************************************
;  Negative  Canonical events precipitation
;***************************************************************************************************************   
  atlN_rain1A = (atlN_rainA(time|0,month|:))
  atlN_rain2A = (atlN_rainA(time|1,month|:))
  atlN_rain3A = (atlN_rainA(time|2,month|:))
  atlN_rain4A = (atlN_rainA(time|5,month|:))

  canonical_events_rainNA = (/atlN_rain1A,atlN_rain2A,atlN_rain3A,atlN_rain4A/)       
  copy_VarCoords_n(atlP_rainA,canonical_events_rainNA,(/0,1/))
  canonical_events_rainNA!0 = "nyears"    
  printVarSummary(canonical_events_rainNA) 
  
   delete([/atlN_rain1A, atlN_rain2A, atlN_rain3A, atlN_rain4A/])
;***************************************************************************************************************
  canonical_events__rain_differenceA = (dim_avg_n_Wrap(canonical_events_rainPA,0) - dim_avg_n_Wrap(canonical_events_rainNA,0))/(2.)
  copy_VarCoords(canonical_eventsN_s,canonical_events__rain_differenceA)                             ; copy coordinate arrays 
  printVarSummary(canonical_events__rain_differenceA)
  printMinMax(canonical_events__rain_differenceA,True)
  
  canonical_events__rain_difference2A = dim_avg_n_Wrap(canonical_events__rain_differenceA,0)
;*************************************************************************************************************** 
  pc1raindA = canonical_events__rain_differenceA  ; Transfer metdata for plotting sig. values        
  pc1raindA = ttest((dim_avg_n_Wrap(canonical_events_rainPA,0)),(dim_stddev_n_Wrap(canonical_events_rainPA,0))^2,5\
               ,(dim_avg_n_Wrap(canonical_events_rainNA,0)),(dim_stddev_n_Wrap(canonical_events_rainNA,0))^2,4,True, False)
  pc1raindA@long_name = "Probability: difference between means"
  printVarSummary(pc1raindA)

  pc1raindA = where(abs(pc1raindA).le.p90,canonical_events__rain_differenceA,pc1raindA@_FillValue)
  printVarSummary(pc1raindA)
  
  delete([/canonical_events_rainPA, canonical_events_rainNA/])
;***************************************************************************************************************
;  Positive  Non-Canonical events precipitation
;***************************************************************************************************************   
  atlP_rain6A = (atlP_rainA(time|0,month|:))
  atlP_rain7A = (atlP_rainA(time|1,month|:))
  atlP_rain8A = (atlP_rainA(time|5,month|:))
  atlP_rain9A = (atlP_rainA(time|8,month|:))
 
  canonical_events_rainHA = (/atlP_rain6A,atlP_rain7A,atlP_rain8A,atlP_rain9A/)       
  copy_VarCoords_n(atlP_rainA,canonical_events_rainHA,(/0,1/))
  canonical_events_rainHA!0 = "nyears"    
 ; printVarSummary(canonical_events_rainHA) 
  
  delete([/atlP_rain6A, atlP_rain7A, atlP_rain8A, atlP_rain9A/])
;***************************************************************************************************************
;  Negative  Non-Canonical events precipitation
;***************************************************************************************************************   
  atlN_rain6A = (atlN_rainA(time|3,month|:))
  atlN_rain7A = (atlN_rainA(time|4,month|:))
  atlN_rain8A = (atlN_rainA(time|6,month|:))
  atlN_rain9A = (atlN_rainA(time|7,month|:))
  
  canonical_events_rainKA = (/atlN_rain6A,atlN_rain7A,atlN_rain8A,atlN_rain9A/)       
  copy_VarCoords_n(atlP_rainA,canonical_events_rainKA,(/0,1/))
  canonical_events_rainKA!0 = "lev"    
 ; printVarSummary(canonical_events_rainKA)  
  
  delete([/atlN_rain6A, atlN_rain7A, atlN_rain8A, atlN_rain9A/])
;-------------------------------------------------------------------------------
;=======================Precipitation==============================================================
  canonical_events__rain_differenceJA = (dim_avg_n_Wrap(canonical_events_rainHA,0) - dim_avg_n_Wrap(canonical_events_rainKA,0))/(2.)
  copy_VarCoords(canonical_eventsN_s,canonical_events__rain_differenceJA)                             ; copy coordinate arrays 
  printVarSummary(canonical_events__rain_differenceJA)
  printMinMax(canonical_events__rain_differenceJA,True)
  
  canonical_events__rain_differenceJ1A = dim_avg_n_Wrap(canonical_events__rain_differenceJA,0)
;*************************************************************************************************************** 
;===============================================================================================================
  pc1raind3A = canonical_events__rain_differenceJA  ; Transfer metdata for plotting sig. values        

  pc1raind3A = ttest((dim_avg_n_Wrap(canonical_events_rainHA,0)),(dim_stddev_n_Wrap(canonical_events_rainHA,0))^2,4\
               ,(dim_avg_n_Wrap(canonical_events_rainKA,0)),(dim_stddev_n_Wrap(canonical_events_rainKA,0))^2,4,True, False)
  pc1raind3A@long_name = "Probability: difference between means"
  printVarSummary(pc1raind3A)

  pc1raind3A = where(abs(pc1raind3A).le.p90,canonical_events__rain_differenceJA,pc1raind3A@_FillValue) 
  printVarSummary(pc1raind3A)
  
  delete([/canonical_events_rainHA, canonical_events_rainKA/]) 
;===============================================================
;***************************************************************************************************************
  rain_can_nonA = new((/2,nmos/),float,canonical_events__rain_differenceA@_FillValue)
  rain_can_nonA(0:1,:) = (/canonical_events__rain_differenceA,canonical_events__rain_differenceJA/)
  printVarSummary(rain_can_nonA)
 
  rain_can_nonA!0 = "lev"
  rain_can_nonA!1 = "month"
  printVarSummary(rain_can_nonA)
 
;***************************************************************************************************************
  rain_can_non_staA = new((/2,nmos/),float,pc1raind3A@_FillValue)
  rain_can_non_staA(0:1,:) = (/pc1raindA,pc1raind3A/)
  printVarSummary(rain_can_non_staA)

  rain_can_non_staA!0 = "lev"
  rain_can_non_staA!1 = "month"
  printVarSummary(rain_can_non_staA) 
  
  printMinMax(canonical_events__rain_differenceA,True)
  printMinMax(canonical_events__rain_differenceJA,True)

;=================================================================================================================================
;=================================================================================================================================
   wks = gsn_open_wks("eps","Figure_1")         ; send graphics to PNG file
;***************************************************************************************************************************************************
   plot = new(4,graphic)                ; create graphic array
;***************************************************************************************************************************************************
;===============================================================================================================
  res                      = True         
  res@gsnDraw              = False        ; don't draw yet
  res@gsnFrame             = False        ; don't advance frame yet
  res@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels
  res@cnMinLevelValF       = -3.            ; set min contour level
  res@cnMaxLevelValF       =  3.             ; set max contour level
  res@cnLevelSpacingF      =  .3            ; set contour spacing   
  res@cnFillOn             = True         ; turn on color fill
  res@cnFillPalette        = "cmp_b2r"  
 
 
  res@cnLinesOn            = False        ; True is default
  res@cnLineLabelsOn       = False        ; True is default
  res@cnInfoLabelOn        = False
  res@lbLabelBarOn         = True       ; turn off individual lb's 
  
  res@lbOrientation        = "vertical"
  res@lbTitleFontHeightF  = 0.020 
  res@lbLabelFontHeightF  = 0.020
  res@lbLeftMarginF       = 0.10
  res@lbRightMarginF      = 0.55
  ;res@lbLabelOffsetF	  = 0.20
  res@lbLabelFont   	  = "Helvetica-Bold"     ; label font
  res@lbTitleString       = "Precipitation (mm day~S~-1~N~)"
  ;res@lbLabelPosition	  = "Top"
  res@lbTitleAngleF       = 270.
  res@lbTitlePosition      = "Right"
  res@lbTitleDirection     = "Across"
  
  
  res@lbLabelStride        = 10 
  
  res@mpFillOn             = True        ; turn off map fill
  res@mpMinLatF            = latS         ; zoom in on map
  res@mpMaxLatF            = latN
  res@mpMinLonF            = lonL
  res@mpMaxLonF            = lonR
  res@pmTickMarkDisplayMode = "Always"
  ;res@tmXBLabelsOn         = False      ;no lon ticks
  ;res@tmYLLabelsOn         = False      ;no lon ticks
  
  ;res@tmYLMode   = "Explicit"
  ;res@tmYLValues = ispan(-15,15,15)
  ;res@tmYLLabels = (/"15~S~o~N~S", "0~S~o~N","15~S~o~N~N"/)
  ;res@tmYLMinorValues  = ispan(-15,15,5) 
  ;res@tmYLLabelFontHeightF = 0.017

  ;res@tmXBMode   = "Explicit"
  ;res@tmXBValues = ispan(-60,20,20)
  ;res@tmXBLabels = (/"60~S~o~N~W","40~S~o~N~W","20~S~o~N~W","0~S~o~N~","20~S~o~N~E"/)
  ;res@tmXBMinorValues  = ispan(-50,20,5)
  ;res@tmXBLabelFontHeightF = 0.017
  ;res@mpLandFillColor	   =  "grey"
  ;res@tmXTOn = True ; Turn off tickmarks and labels
  ;res@tmYROn = False 	
  ;res@tmXBOn = True ; Turn off tickmarks and labels
  ;res@tmYLOn = False 	    
 
  res@mpGeophysicalLineColor = "grey" 
  res@mpFillDrawOrder      = "postDraw"

  res@vpHeightF 	  = 0.80        ; Changes the aspect ratio
  res@vpWidthF 		  = 0.80
  res@gsnCenterStringFontHeightF   = 0.030
    

  res@gsnCenterString  = ""   
  res@gsnLeftString  = ""  
  res@gsnRightString  = ""  
    
    
 ; res@vcRefAnnoOrthogonalPosF = -1.0             ; move ref vector up
 ; res@vcRefMagnitudeF         = 0.02             ; define vector ref mag
 ; res@vcRefAnnoString1          = "0.02 N m~S~-2~N~ "
 ; res@vcRefLengthF            = 0.02            ; define length of vec ref
 ; res@vcRefLengthF             = 0.035            ; define length of vec ref
 ; vcres@vcGlyphStyle            = "CurlyVector"    ; turn on curly vectors
 ; res@vcMinDistanceF          = 0.030
 ; res@vcLineArrowThicknessF   = 1.
 ; res@vcLineArrowColor        = "black"           ; change vector color
 ; res@gsnScalarContour	      = True
 ; res@vcLineArrowHeadMaxSizeF = 0.05
 ; res@vcLineArrowHeadMinSizeF = 0.01
 ; res@vcMaxLevelCount         = 100
  
 
 ;---For higher quality large images
  res@tiYAxisFontThicknessF  = 1.
  res@tmXBMajorThicknessF    = 1.
  res@tmXBMinorThicknessF    = 1.
  res@tmYLMajorThicknessF    = 1.
  res@tmYLMinorThicknessF    = 1.
  res@tmYLLabelFontThicknessF= 1.
  res@tmBorderThicknessF     = 1. 
;=============================================================================
;=============================================================================

  res1                      =   True
  res1@gsnDraw              = False        ; don't draw yet
  res1@gsnFrame             = False        ; don't advance frame yet
  res1@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels
  res1@cnMinLevelValF       = res@cnMinLevelValF            ; set min contour level
  res1@cnMaxLevelValF       = res@cnMaxLevelValF             ; set max contour level
  res1@cnLevelSpacingF      = res@cnLevelSpacingF           ; set contour spacing           ; set contour spacing  
  res1@gsnAddCyclic         = False   
  res1@cnLinesOn            = False       ; True is default
  res1@cnLineLabelsOn       = False        ; True is default
  res1@cnInfoLabelOn        = False
  res1@gsnRightString       = ""
  res1@gsnLeftString        = ""
  res1@cnFillScaleF         =  1.5         ; add extra density
  ;res1@cnLineLabelDensityF =  -10.
  res1@cnFillDotSizeF       = 0.0001            ; increase dot size of contours
;-----------------------------------------------------------------------------
  opt = True
  opt@gsnShadeFillType = "pattern"        ; color is the default
  opt@gsnShadeLow = 17                     ; use fill pattern #17
  opt@gsnShadeHigh = 17                     ; use fill pattern #17
  opt@gsnShadeFillScaleF = 1.0
 ; opt@gsnShadeFillDotSizeF = 0.001
  
   delete(res@cnFillPalette)
   res@cnFillPalette        	= "amwg_blueyellowred"       ;different colour bar
   res@gsnStringFontHeightF = 0.013
   res@gsnLeftString  		= "Canonical events"
;============================================================================= 
   plotA07 = gsn_csm_contour_map_ce(wks,SST_C(:,:),res)  ; 
   plotS07 = gsn_csm_contour(wks,SSTdP_C(:,:),res1)       ;
   plotS07 = gsn_contour_shade(plotS07,p90,SSTdP_C@_FillValue,opt)
   overlay(plotA07,plotS07)
   plot(2) = plotA07
   
   res@gsnLeftString  		= "Non-canonical events"
   plotA17 = gsn_csm_contour_map_ce(wks,SST_N(:,:),res)  ;           ; 
   plotS17 = gsn_csm_contour(wks,SSTdP_N(:,:),res1)       ;
   plotS17 = gsn_contour_shade(plotS17,p90,SSTdP_N@_FillValue,opt)
   overlay(plotA17,plotS17)
   plot(3) = plotA17
   
;***************************************************************************************************************************************************  
   r2               = True                        ; plot mods desired
   r2@gsnDraw       = False                       ; do not draw
   r2@gsnFrame      = False                       ; do not advance frame
   r2@tiMainString  = ""                       ; title
   r2@tiXAxisFontHeightF = 0.015
   r2@tiXAxisString = "Calendar month"  ; xaxis
   r2@tiYAxisFontHeightF = 0.015
   r2@tiYAxisString = "Atl3 sea surface temperature (~S~o~N~C)"                  ; yaxis
  
   r2@tmXBLabelFontHeightF = 0.012
   r2@tmYLLabelFontHeightF = 0.012
   r2@trYLog              = False                ; log scaling
   r2@trXLog              = False               ; log scaling
   
    
   r2@trXMinF                     = min(xlabel)        ; set minimum X-axis value
   r2@trXMaxF                     = max(xlabel)        ; set maximum X-axis value
   
   r2@trYMinF                     =  -0.2                 ; manually set lower limit
   r2@trYMaxF                     =   1.0           

   
   r2@gsnXRefLineThicknesses = (/0.5/)
   r2@gsnXRefLineDashPatterns = (/2./)
   r2@xyLineColors        = (/"red","transparent"/)
   r2@tmXTBorderOn            =  True
   r2@tmXTOn                  =  True
      
   r2@tmYLFormat = "@*+^sg"	; use "e" symbol, don't force 0's (see documentation)
   ;r2@gsnXRefLine            =  (/1988./) 
   
   r2@vpWidthF                = 0.60
   r2@vpHeightF               = 0.40

   r2@xyLineThicknesses   = (/1.,1.,1./)       ; Define line thicknesses 
   r2@xyDashPatterns      = (/0.,2.,2./)           ; Dash patterns 
   r2@tmYLLabelFontColor = "red"   ; and right axes are same.
   r2@tiYAxisFontColor  = "red"
   
      r2@tmXBMode                    = "Explicit"   ; Define your own tick mark labels.
 ; r2@tmXBLabelFont               = 25           ; Change font of labels.
   r2@tmXBLabelFontHeightF        = 0.015 ; Change font height of labels.
   r2@tmXBMinorOn                 = False        ; No minor tick marks.
   r2@tmXBValues                  = ispan(0,12,1); Values from 0 to 12.
   r2@tmXBLabels                  = (/"Jan","Feb","Mar","Apr","May","Jun",\
                               "Jul","Aug","Sep","Oct","Nov","Dec","Jan"/)
   
   ;===legends
   r2@pmLegendDisplayMode    = "Always"            ; turn on legend
   r2@pmLegendSide           = "Top"               ; Change location of 
   r2@pmLegendParallelPosF   = 0.78                 ; move units right
   r2@pmLegendOrthogonalPosF = -0.33                ; move units down
   r2@pmLegendWidthF         = 0.12                ; Change width and
   r2@pmLegendHeightF        = 0.10                ; height of legend.
   r2@lgLabelFontHeightF     = 0.010                 ; change font height
   r2@lgTitleOn              = False               ; turn on legend title
   r2@lgTitleString          = "SST indices"           ; create legend title
   r2@lgTitleFontHeightF     = .003                ; font of legend title
   r2@lgPerimOn              = False               ; turn off box around
   r2@xyExplicitLegendLabels = (/"Atl3-SSTA"/)         ; explicit labels
   		
   				
;----------------------------------------------------------------------------------------------------------------   
;*************************************************************************************************************** 
   r1               = True                        ; plot mods desired
   r1@gsnDraw       = False                       ; do not draw
   r1@gsnFrame      = False                       ; do not advance frame
   r1@tiXAxisString = ""  ; xaxis
   r1@tiYAxisString = "Precipitation anomalies (mm day~S~-1~N~)"      
   r1@tiYAxisAngleF  =  270. 
   r1@xyLineColors   =  (/"blue","deepskyblue","transparent"/)
   r1@xyLineThicknesses  =  (/1.,1.,1./)
   r1@xyDashPatterns     =  (/0.,0.,0./)
   r1@tiXAxisFontHeightF = 0.015
   r1@tiYAxisFontHeightF = 0.015
   r1@tmXBLabelFontHeightF = 0.012
   r1@tmYLLabelFontHeightF = 0.012
       
   r1@tmYUseLeft          =  True
   r1@trYMinF                     =  -.5                 ; manually set lower limit
   r1@trYMaxF                     =   2.5           
   
   r1@trXMinF                     = min(xlabel)        ; set minimum X-axis value
   r1@trXMaxF                     = max(xlabel)        ; set maximum X-axis value
     
   r1@vpWidthF                = 0.60
   r1@vpHeightF               = 0.40
   r1@tmXTBorderOn            =  True
   r1@tmXTOn                  =  False
   ;r1@tmXTMode = "Explicit"
   ;r1@tmXTValues = (/"0", "1","2","3","4","5","6","7","8","9","10","11","12"/)
   ;r1@tmXTLabels = (/"1", "2","3","4","5","6","7","8","10","10","11","12"/)
   
  r1@tmYRLabelFontColor = "blue"   ; and right axes are same.
  r1@tiYAxisFontColor  = "blue"
  ;r1@tiYAxisString   = "SSH anomaly (cm)"
  
  r1@tmXTLabelsOn        = False    ; have tick mark labels
  r1@tmXTOn              = True      ; have tick marks
  r1@tmYRLabelFontColor = "blue"   ; and right axes are same.
  r1@tiYAxisFontColor  = "blue"

  r1@tmYLBorderOn 	 = False
  r1@tmYRBorderOn	 = True
  r1@tmXTBorderOn	 = False
  r1@tmXBBorderOn	 = False
 
  r1@gsnYRefLineColor     = "blue"
  r1@tmYRFormat   	=  "@*+^sg"          ; trim the uneccecary "0" in labels. 
  
    ;===legends
   r1@pmLegendDisplayMode    = "Always"            ; turn on legend
   r1@pmLegendSide           = "Top"               ; Change location of 
   r1@pmLegendParallelPosF   = 0.82                 ; move units right
   r1@pmLegendOrthogonalPosF = -0.35                ; move units down
   r1@pmLegendWidthF         = 0.12                ; Change width and
   r1@pmLegendHeightF        = 0.10                ; height of legend.
   r1@lgLabelFontHeightF     = 0.010                 ; change font height
   r1@lgTitleOn              = False               ; turn on legend title
   r1@lgTitleString          = "SST indices"           ; create legend title
   r1@lgTitleFontHeightF     = .003                ; font of legend title
   r1@lgPerimOn              = False               ; turn off box around
   r1@xyExplicitLegendLabels = (/"Eq. Atl. Precipitation"/)         ; explicit labels
   		
      
;*************************************************************************************************************************************************** 
;*************************************************************************************************************** 
;========================================================
  r3             = True
  r3@trYMinF             =  0.0                 ;   "          upper
  r3@trYMaxF             =  1.0                 ;   "          upper
  r3@xyDashPatterns     =  (/0.,0.,0./)
  r3@tmXBLabelFontHeightF = 0.012
  r3@tmYLLabelFontHeightF = 0.012
  r3@tmYLLabelFontColor = "deepskyblue"   ; and right axes are same.
  r3@tiYAxisFontColor  = "deepskyblue"
  r3@xyLineColors        = (/"deepskyblue","red","transparent"/)
  r3@tmYRLabelFontColor = "deepskyblue"   ; and right axes are same.
  r3@tiYAxisFontColor  = "deepskyblue"
  r3@tiYAxisString = "Precipitation anomalies (mm day~S~-1~N~)"                  ; yaxis
  r3@tmYRFormat = "@*+^sg"
  r3@trYMinF             =  -0.60                 ; manually set lower limit
  r3@trYMaxF             =  2.0                 ;   "          upper
  
  r3@tiYAxisAngleF  =  270. 
  
 
     ;===legends
   r3@pmLegendDisplayMode    = "Always"            ; turn on legend
   r3@pmLegendSide           = "Top"               ; Change location of 
   r3@pmLegendParallelPosF   = 0.82                 ; move units right
   r3@pmLegendOrthogonalPosF = -0.18                ; move units down
   r3@pmLegendWidthF         = 0.12                ; Change width and
   r3@pmLegendHeightF        = 0.10                ; height of legend.
   r3@lgLabelFontHeightF     = 0.010                 ; change font height
   r3@lgTitleOn              = False               ; turn on legend title
   r3@lgTitleString          = "SST indices"           ; create legend title
   r3@lgTitleFontHeightF     = .003                ; font of legend title
   r3@lgPerimOn              = False               ; turn off box around
   r3@xyExplicitLegendLabels = (/"G. Coast precipitation"/)         ; explicit labels
   		
;======================================================== 
;===============================================================================================================   
     
   plotAS   =  gsn_csm_xy3(wks,xlabel,x_can_non(lev|0,month|:),rain_can_non(lev|0,month|:),rain_can_nonA(lev|0,month|:),r2,r1,r3)   
   plotAS1  =  gsn_csm_xy3(wks,xlabel,x_can_non(lev|1,month|:),rain_can_non(lev|1,month|:),rain_can_nonA(lev|1,month|:),r2,r1,r3)   
;===============================================================================================================
;*************************************************************************************************************** 
;===============================================================================================================
   ; r2@trYMinF             =  -0.2                 ; manually set lower limit
  ; r2@trYMaxF             =   0.2           

   r2@pmLegendDisplayMode    = "Always"            ; turn on legend
  ; r2@xyDashPatterns      = (/0.,0.,0.,0.,0.,0.,0.,0.,0./)           ; Dash patterns 
 
   r2@pmLegendDisplayMode    = "Never"            ; turn on legend
 

 ;  delete(r2@xyLineColors)
   r2@xyMarkLineMode	 =   "Markers"
   r2@xyMarkers     =  (/16.,16.,16.,16.,16.,16.,16.,16.,16./)
   r2@xyMarkerColors =  xy_colors_dheating
   r2@xyMarkerSizeF   = 0.010
   r2@pmLegendDisplayMode    = "Never"            ; turn on legend
   r2@gsnStringFontHeightF = 0.012
 ;  r2@gsnLeftString  		= "Difference in events"
;===============================================================================================================
   ; r1@trYMinF             =  -0.2                 ; manually set lower limit
  ; r1@trYMaxF             =   0.2           

   r1@pmLegendDisplayMode    = "Always"            ; turn on legend
 ;  r1@xyDashPatterns      = (/0.,0.,0.,0.,0.,0.,0.,0.,0./)           ; Dash patterns 
 
   r1@pmLegendDisplayMode    = "Never"            ; turn on legend
 
   r1@xyMarkLineMode	 =   "Markers"
   r1@xyMarkers     =  (/16.,16.,16.,16.,16.,16.,16.,16.,16./)
   r1@xyMarkerColors =  xy_colors_dheating2
   r1@xyMarkerSizeF   = 0.010
   r1@pmLegendDisplayMode    = "Never"            ; turn on legend
   r1@gsnStringFontHeightF = 0.012
 ;  r1@gsnLeftString  		= "Difference in events" 
   
;===============================================================================================================
   ; r3@trYMinF             =  -0.2                 ; manually set lower limit
  ; r3@trYMaxF             =   0.2           

   r3@pmLegendDisplayMode    = "Always"            ; turn on legend
  ; r3@xyDashPatterns      = (/0.,0.,0.,0.,0.,0.,0.,0.,0./)           ; Dash patterns 
 
   r3@pmLegendDisplayMode    = "Never"            ; turn on legend
 
   r3@xyMarkLineMode	 =   "Markers"
   r3@xyMarkers     =  (/16.,16.,16.,16.,16.,16.,16.,16.,16./)
   r3@xyMarkerColors =  xy_colors_dheating1
   r3@xyMarkerSizeF   = 0.010
   r3@pmLegendDisplayMode    = "Never"            ; turn on legend
   r3@gsnStringFontHeightF = 0.012
 ;  r3@gsnLeftString  		= "Difference in events" 
      
;***************************************************************************************************************************************************  
   plotDS   =  gsn_csm_xy3(wks,xlabel,x_can_non_sta(lev|0,month|:),rain_can_non_sta(lev|0,month|:),rain_can_non_staA(lev|0,month|:),r2,r1,r3)
   
   plotDS1  =  gsn_csm_xy3(wks,xlabel,x_can_non_sta(lev|1,month|:),rain_can_non_sta(lev|1,month|:),rain_can_non_staA(lev|1,month|:),r2,r1,r3)
;*************************************************************************************************************** 
;===============================================================================================================  
  overlay(plotAS,plotDS)
  plot(0) = plotAS
  
  overlay(plotAS1,plotDS1)
  plot(1) = plotAS1
;*************************************************************************************************************** 
;=============================================================================================================== 
  x_0 = 5.
  x_1 = 5.
  y_0 = -2.5 
  y_1 = 2.5 
 
 ;---Resources for the polyline
  pres1                  = True           ; polyline mods desired
  pres1@gsLineThicknessF = 0.8            ; line thickness   
  pres1@gsLineColor      = "black"         ; color of lines
  pres1@gsLineDashPattern = 2
  pres1@tfPolyDrawOrder   = "postDraw"
  pres1@gsLineColor      = "black"         ; color of lines
  dumL0 = gsn_add_polyline(wks,plot(0), (/x_0,x_1/),(/y_0,y_1/) ,pres1)
  dumL2 = gsn_add_polyline(wks,plot(1), (/x_0,x_1/),(/y_0,y_1/) ,pres1)
  x_2 = 0.
  x_3 = 11.
  y_2 = 0. 
  y_3 = 0.  
  
  dumL1 = gsn_add_polyline(wks,plot(0), (/x_2,x_3/),(/y_2,y_3/) ,pres1)
  dumL3 = gsn_add_polyline(wks,plot(1), (/x_2,x_3/),(/y_2,y_3/) ,pres1)
;***************************************************************************************************************************************************
;=====================================================================
  txres5               = True
  txres5@txAngleF      = 0.          ; Rotate string clockwise 
  txres5@txFontHeightF = 0.025   ; Use same font height as left axis
  txres5@txFont                = "Helvetica-Bold"
 ; txid = gsn_create_text(wks, "DT anomaly (~S~o~N~C)", txres)
   txid5 = gsn_create_text(wks, "a", txres5);
; Move text string to center/right edge of plot.;
  amres5                  = True
  amres5@amParallelPosF   =  -0.5 ; 0.5 is the right edge of the plot, so
                                 ; 0.6 is a little further to the right.
  amres5@amOrthogonalPosF =  -0.5   ; This is the center of the plot.
 ; tt_p0X = gsn_add_annotation(plot(0),txid5,amres5)
  
  txres5@txFontHeightF = 0.030   ; Use same font height as left axis
  amres5@amParallelPosF   =  -0.63 ; 0.5 is the right edge of the plot, so
  amres5@amOrthogonalPosF =  -0.55    ; This is the center of the plot.
  amres5@amJust           = "TopLeft"

 tt = (/"a","b","c","d","e","f","g","g"/)
create_tt = gsn_create_text(wks,tt,txres5);
; tt_p0 = new(nplots,graphic)                ; create graphic array ; only needed if paneling
tt_p0 = new(4,graphic)                ; create graphic array ; only needed if paneling
;---Attach the polyline
     do i = 0,3
	print(i)
        tt_p0(i) = gsn_add_annotation(plot(i),create_tt(i),amres5)
     end do
delete(i)
;=============================================================================  
;===============================================================
   resP             			  = True                        ; panel mods desired
   resP@gsnDraw              		  = True       ; don't draw yet
   resP@gsnFrame         		  = False        ; don't advance frame yet
   resP@gsnPanelYWhiteSpacePercent       = 0.0    ; default is 1.0
   resP@gsnPanelXWhiteSpacePercent       = 0.0    ; default is 1.0
   resP@gsnPanelTop                      = 0.65
   resP@gsnPanelBottom                   = 0.28
   resP@gsnPanelLeft                     = 0.10
   resP@gsnPanelRight                    = 0.95
   ;resP@gsnPanelTop                      = 0.99
   ;resP@gsnPanelBottom                   = 0.60
   ;resP@gsnPanelLeft                     = 0.10
   ;resP@gsnPanelRight                    = 0.90
   resP@gsnPanelFigureStringsPerimOn 	  = False 
   resP@gsnPaperOrientation 		  = "portrait"
   resP@gsnPanelDebug      		  = True
   resP@amJust   			  = "BottomLeft"
  ; resP@gsnPanelFigureStringsPerimOn     = False
  ; resP@gsnPanelFigureStringsFontHeightF = 0.010
  ; resP@gsnPanelFigureStringsBackgroundFillColor = "transparent"
  ; resP@gsnPanelFigureStrings= (/"lag = 0","lag = -1","c","d","e","f","g","h","i","j","k","l"/) ; add strings to panel
   resP@gsnMaximize = False                        ; blow up plot
;  gsn_panel(wks,plot(2:),(/1,2/),resP)               ; create panel plots

;===============================================================
;***********************************************
  resP                                  = True         ; modify the panel plot
  resP@gsnDraw                          = True      ; don't draw yet
  resP@gsnFrame                         = False        ; don't advance frame yet
  resP@gsnPanelLabelBar                 = False        ; add common colorbar
  resP@gsnPanelYWhiteSpacePercent       = 1.0    ; default is 1.0
  resP@gsnPanelXWhiteSpacePercent       = 1.0    ; default is 1.0
  resP@gsnPanelTop                      = 0.89
  resP@gsnPanelBottom                   = 0.59
  resP@gsnPanelLeft                     = 0.05
  resP@gsnPanelRight                    = 0.99
 ; resP@gsnMaximize                      = True         ; large format
  resP@gsnPaperOrientation              = "portrait"   ; force portrait
  resP@gsnLeftString                    = ""
  resP@lbLabelStride                    = 5
  ;resP@pmLabelBarParallelPosF          = 0.50
  ;resP@pmLabelBarOrthogonalPosF        = -0.11
  resP@lbOrientation                    = "vertical"
  resP@lbTitleFontHeightF               = 0.008
  resP@lbLabelFontHeightF               = 0.008
  resP@lbLeftMarginF                    = 0.10
  resP@lbRightMarginF                   = 0.55
  resP@lbTopMarginF                     = -0.65
  resP@lbBottomMarginF                  = -0.55
 ;resP@lbLabelOffsetF                   = 0.20
  resP@lbLabelFont                      = "Helvetica-Bold"     ; label font
  resP@lbTitleFont                      = "Helvetica-Bold"     ; label font
  ;resP@lbLabelPosition                 = "Top"
  resP@lbTitleAngleF                    = 270.
  resP@lbTitlePosition                  = "Right"
  resP@lbTitleDirection                 = "Across"
  resP@lbTitleString                    = "Sea surface temperature (~S~o~N~C)"
  resP@gsnPanelFigureStrings            = (/"Canonical events","Non-canonical events"/) ; add strings to panel
  resP@amJust                           = "TopLeft"
  resP@gsnPanelFigureStringsFontHeightF =  0.008
  resP@gsnPanelFigureStringsFont        =  "Helvetica-Bold"

  resP@txFont                           = "Helvetica-Bold"
 ;resP@txString                        = "SST trends: 1984-2008"
  ;resP@txString                        = "JAS SST trends: 1979-2018"
;**************************************************************
  gsn_panel(wks,plot(0:1),(/1,2/),resP)     ; now draw as one plot
;======================================================
;======================================================
  frame(wks)
end
