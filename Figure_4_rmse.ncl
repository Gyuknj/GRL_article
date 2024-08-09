begin
;======================================
  yrStrt = 1982
  yrLast = 2022

  nmos = 12	;months in a year
  thres = 1.0	;set a threshold
 
  xlabel = ispan(0,nmos-1,1)
 ; print(xlabel)
  
  p95 = 0.1
  diri = "/home/gyuk/Downloads/GYUK/GEOMAR/Final_work/"
  
  sst_data = "HadISST_sst.nc"
;==============================================================
;==============================================================
  f  =   addfile(diri+sst_data,"r")
  TIME   = f->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  hadi   = f->sst(iYYYY,:,:)
  hadi!1 = "lat"
  hadi!2 = "lon"
  printVarSummary(hadi)       ; variable overview
		
  delete(f)		;	Not needed anymore
  delete(TIME)
  delete(YYYY)
  delete(iYYYY)

  ;printMinMax(hadi,True)
  hadi = where(hadi.gt.-30.,hadi,hadi@_FillValue)	;removes padded values
  ;printVarSummary(hadi)       ; variable overview
  ;printMinMax(hadi,True)
;==========================
  hadiClm      = clmMonTLL(hadi)
  ;printVarSummary(hadiClm)                   ; (12,nlat,nlon)   
  
  hadi         = calcMonAnomTLL(hadi,hadiClm)    ; replace with anonamlies
  hadi@long_name = "ANOMALIES: "+hadi@long_name
  ;delete(hadiClm)
  ;printMinMax(hadi,True)
;=========================
;----------------------------
  nyears = dimsizes(hadi&time)/nmos
  ;print(nyears)

  hadi_s     = onedtond(ndtooned(hadi),(/nyears,nmos,dimsizes(hadi&lat),dimsizes(hadi&lon)/))
  copy_VarCoords(hadi(0,:,:),hadi_s(0,0,:,:))                             ; copy coordinate arrays
  hadi_s!0 = "time"
  hadi_s&time = fspan(yrStrt,yrLast,nyears)
  hadi_s!1 = "month"
  printVarSummary(hadi_s)

  hadi_s = dtrend_msg_n(hadi_s&time,hadi_s,True,False,0)        ;removes linear trends
  ;printVarSummary(hadi_s)
  ;printMinMax(hadi_s,True)

;---------------------------------------------------------------------------
  ahadi_s = wgt_areaave_Wrap(hadi_s(time|:,month|:,{lon|-20:0}, {lat|-3:3}),1.0, 1.0, 0)        ;atl3 index
 ; printVarSummary(ahadi_s)
  ;printMinMax(ahadi_s,True)
;----------------------------
  delete(hadi)

  jja_ahadi = dim_avg_n_Wrap(ahadi_s(time|:,month|5:7),1)		;JJA, ncl counts from zero
 ; print(jja_ahadi)

   jja_ahadi = dim_standardize_n(jja_ahadi, 0, 0)  ;Normalise deviations
 ;  print(jja_ahadi)
;----------------------------------------------------------------------------
; Positive cases
;============================================================================
  n_atlP = dim_num_n(jja_ahadi.ge.thres,0)        ;Number of cases
 ; print(n_atlP)

  valsP = ind(jja_ahadi.ge.thres)
 ; print(valsP)
 ; print(valsP+yrStrt)		;the years
;----------------------------------------------------------------------------
; Negative cases,
; note the "-" sign attached to the threshold here
;============================================================================
  n_atlN = dim_num_n(jja_ahadi.le.-thres,0)        ;Number of cases
 ; print(n_atlN)

  valsN = ind(jja_ahadi.le.-thres)
 ; print(valsN)
 ; print(valsN+yrStrt)             ;the years
;==============================================================
  dimZ = dimsizes(hadi_s)
  printVarSummary(dimZ)
        
        saodP_cases_maps  = new((/n_atlP,dimZ(1),dimZ(2),dimZ(3)/),float,hadi_s@_FillValue)
      ;  printVarSummary(saodP_cases_maps)
          
        saodP_cases_maps = hadi_s(valsP,:,:,:)
      ;  printVarSummary(saodP_cases_maps)        
;==============================================================
;   Positive Canonical event years
;==============================================================        
        sstb1p = (saodP_cases_maps(time|2,month|:,lat|:,lon|:))	;now based on the basin average
    ;    printVarSummary(sstb1p)    
  
        sstb2p = (saodP_cases_maps(time|3,month|:,lat|:,lon|:))	;now based on the basin average
     ;   printVarSummary(sstb2p)    
  
        sstb3p = (saodP_cases_maps(time|4,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(sstb3p)    
  
        sstb4p = (saodP_cases_maps(time|6,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(sstb4p) 
  
        sstb5p = (saodP_cases_maps(time|7,month|:,lat|:,lon|:))	;now based on the basin average
     ;   printVarSummary(sstb5p) 
         
       sst_fp = (/sstb1p,sstb2p,sstb3p,sstb4p,sstb5p/)       
 
       copy_VarCoords_n(hadi_s,sst_fp,(/1,2,3/))
       sst_fp!0 = "time"    
       printVarSummary(sst_fp)  
;==============================================================
;   Negative canonical event years 
;==============================================================     
        saodN_cases_maps  = new((/n_atlN,dimZ(1),dimZ(2),dimZ(3)/),float,hadi_s@_FillValue)
       ; printVarSummary(saodN_cases_maps)
    
        saodN_cases_maps = hadi_s(valsN,:,:,:)
      ;  printVarSummary(saodN_cases_maps)
        
        sstb1n = (saodN_cases_maps(time|0,month|:,lat|:,lon|:))	;now based on the basin average
       ; printVarSummary(sstb1n)    
  
        sstb2n = (saodN_cases_maps(time|1,month|:,lat|:,lon|:))	;now based on the basin average
       ; printVarSummary(sstb2n)    
  
        sstb3n = (saodN_cases_maps(time|2,month|:,lat|:,lon|:))	;now based on the basin average
       ; printVarSummary(sstb3n)    
  
        sstb4n = (saodN_cases_maps(time|5,month|:,lat|:,lon|:))	;now based on the basin average
       ; printVarSummary(sstb4n)   
      
       sst_fn = (/sstb1n,sstb2n,sstb3n,sstb4n/)        
       copy_VarCoords_n(hadi_s,sst_fn,(/1,2,3/))
       sst_fn!0 = "time"   
       printVarSummary(sst_fn)   
;==============================================================
;   Canonical difference
;============================================================================================
  canonical_diff = (dim_avg_n_Wrap(sst_fp,0) - dim_avg_n_Wrap(sst_fn,0))/(2.)
  copy_VarCoords(hadi_s(0,:,:,:),canonical_diff)                             ; copy coordinate arrays 
  printVarSummary(canonical_diff)
;---------------------------------------------------------------------------
 ; atl3_can =  dim_avg_n_Wrap(canonical_diff(month|5:7,{lon|-20:0}, {lat|-3:3}),0)        ;atl3 index
   wgt = dim_avg_n_Wrap(canonical_diff(month|5:7,{lon|-20:0}, {lat|-3:3}),0)       ;atl3 index
  atl3_can =  canonical_diff(month|5:7,{lon|-20:0}, {lat|-3:3})        ;atl3 index
 ; atl3_can =  canonical_diff(month|5:7,lon|:, lat|:)        ;atl3 index
  printVarSummary(atl3_can)
  printVarSummary(wgt)
;===================================== 
; Non-canonical difference
;==============================================================        
        sstb1p1 = (saodP_cases_maps(time|0,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(sstb1p1)    
  
        sstb2p1 = (saodP_cases_maps(time|1,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(sstb2p1)    
  
        sstb3p1 = (saodP_cases_maps(time|5,month|:,lat|:,lon|:))	;now based on the basin average
       ; printVarSummary(sstb3p1)    
  
        sstb4p1 = (saodP_cases_maps(time|8,month|:,lat|:,lon|:))	;now based on the basin average
       ; printVarSummary(sstb4p1) 
        
       sst_fp1 = (/sstb1p1,sstb2p1,sstb3p1,sstb4p1/)       
 
       copy_VarCoords_n(hadi_s,sst_fp1,(/1,2,3/))
       sst_fp1!0 = "time"    
       printVarSummary(sst_fp1) 
;===================================================================       
        sstb1n1 = (saodN_cases_maps(time|3,month|:,lat|:,lon|:))	;now based on the basin average
     ;   printVarSummary(sstb1n1)    
  
        sstb2n1 = (saodN_cases_maps(time|4,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(sstb2n1)    
  
        sstb3n1 = (saodN_cases_maps(time|6,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(sstb3n1)    
  
        sstb4n1 = (saodN_cases_maps(time|7,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(sstb4n1) 
       
       sst_fn1 = (/sstb1n1,sstb2n1,sstb3n1,sstb4n1/)        
       copy_VarCoords_n(hadi_s,sst_fn1,(/1,2,3/))
       sst_fn1!0 = "time"  
       printVarSummary(sst_fn1) 
;==============================================================
  non_canonical_diff = (dim_avg_n_Wrap(sst_fp1,0) - dim_avg_n_Wrap(sst_fn1,0))/(2.)
  copy_VarCoords(hadi_s(0,:,:,:),non_canonical_diff)                             ; copy coordinate arrays 
  printVarSummary(non_canonical_diff)
;----------------------------------------------------------------------------------------------------------------------
 ; atl3_non =  dim_avg_n_Wrap(non_canonical_diff(month|5:7,{lon|-20:0}, {lat|-3:3}),0)        ;atl3 index
  atl3_non =  non_canonical_diff(month|5:7,{lon|-20:0}, {lat|-3:3})        ;atl3 index
 ; atl3_non =  non_canonical_diff(month|5:7,lon|:, lat|:)        ;atl3 index
  printVarSummary(atl3_non)
;----------------------------------------------------------------------------------------------------------------------
  f      = addfile("/home/gyuk/Downloads/GYUK/GEOMAR/Final_work/sst_ensemble_data_results/regrid_sst1_mean.nc", "r")
  TIME   = f->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  sst_0   = short2flt(f->sst(iYYYY,:,:))
  printVarSummary(sst_0)       ; variable overview
  sst_1 = sst_0 - 273.15
  copy_VarCoords(sst_0,sst_1)                             ; copy coordinate arrays
  sst_1!0 = "time"
  sst_1!1 = "lat"
  sst_1!2 = "lon"
  printVarSummary(sst_1)       ; variable overview
  printMinMax(sst_1,True)
;====================================================================================================================================== 
  delete(f)		;	Not needed anymore
  delete(TIME)
  delete(YYYY)
  delete(iYYYY)
;==========================
  sst_1Clm      = clmMonTLL(sst_1)
  printVarSummary(sst_1Clm)                   ; (12,nlat,nlon)   
 
  sst_1         = calcMonAnomTLL(sst_1,sst_1Clm)    ; replace with anonamlies
  sst_1@long_name = "ANOMALIES: "+sst_1@long_name
  delete(sst_1Clm)
  ;printMinMax(sst_1,True)
  printVarSummary(sst_1)
;=========================
;----------------------------
  nyears = dimsizes(sst_1&time)/nmos
  ;print(nyears)

  sst_1_s     = onedtond(ndtooned(sst_1),(/nyears,nmos,dimsizes(sst_1&lat),dimsizes(sst_1&lon)/))
  copy_VarCoords(sst_1(0,:,:),sst_1_s(0,0,:,:))                             ; copy coordinate arrays
  sst_1_s!0 = "time"
  sst_1_s&time = fspan(yrStrt,yrLast,nyears)
  sst_1_s!1 = "month"
  printVarSummary(sst_1_s)

  sst_1_s = dtrend_msg_n(sst_1_s&time,sst_1_s,True,False,0)        ;removes linear trends
  ;printVarSummary(sst_1_s)
  ;printMinMax(sst_1_s,True)
;==============================================================
;=================== Canonical
;==============================================================
  dimZ = dimsizes(sst_1_s)
  printVarSummary(dimZ)
        raindP0_cases_maps  = new((/n_atlP,dimZ(1),dimZ(2),dimZ(3)/),float,sst_1_s@_FillValue)
     ;   printVarSummary(raindP0_cases_maps)
    
        raindP0_cases_maps = sst_1_s(valsP,:,:,:)
     ;   printVarSummary(raindP0_cases_maps)
;==============================================================
;   Positive Canonical event years
;==============================================================                    
        rainb1p = (raindP0_cases_maps(time|2,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb1p)    
  
        rainb2p = (raindP0_cases_maps(time|3,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb2p)    
  
        rainb3p = (raindP0_cases_maps(time|4,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb3p)    
  
       rainb4p = (raindP0_cases_maps(time|6,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb4p) 
  
        rainb5p = (raindP0_cases_maps(time|7,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb5p) 
         
       rain_fp = (/rainb1p,rainb2p,rainb3p,rainb4p,rainb5p/)       
 
       copy_VarCoords_n(sst_1_s,rain_fp,(/1,2,3/))
       rain_fp!0 = "time"   
       printVarSummary(rain_fp)  
;==============================================================
;   Negative Canonical event years
;==============================================================              
        raindN0_cases_maps  = new((/n_atlN,dimZ(1),dimZ(2),dimZ(3)/),float,sst_1_s@_FillValue)
        printVarSummary(raindN0_cases_maps)
    
        raindN0_cases_maps = sst_1_s(valsN,:,:,:)
      ;  printVarSummary(raindN0_cases_maps)
        
        rainb1n = (raindN0_cases_maps(time|0,month|:,lat|:,lon|:))	;now based on the basin average
     ;   printVarSummary(rainb1n)    
  
        rainb2n = (raindN0_cases_maps(time|1,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb2n)    
  
        rainb3n = (raindN0_cases_maps(time|2,month|:,lat|:,lon|:))	;now based on the basin average
       ; printVarSummary(rainb3n)    
  
        rainb4n = (raindN0_cases_maps(time|5,month|:,lat|:,lon|:))	;now based on the basin average
       ; printVarSummary(rainb4n) 
       
       rain_fn = (/rainb1n,rainb2n,rainb3n,rainb4n/)       
       copy_VarCoords_n(sst_1_s,rain_fn,(/1,2,3/))
       rain_fn!0 = "time"   
       printVarSummary(rain_fn) 
;==============================================================
  canonical_diff1 = (dim_avg_n_Wrap(rain_fp,0) - dim_avg_n_Wrap(rain_fn,0))/(2.)                       ; Difference
  copy_VarCoords(sst_1_s(0,:,:,:),canonical_diff1)                             ; copy coordinate arrays 
  printVarSummary(canonical_diff1)
  printMinMax(canonical_diff1,True)
;-----------------------------------------------------------------------------------------------------------------------
 ; atl3_lag1_can =  dim_avg_n_Wrap(canonical_diff1(month|5:7,{lon|-20:0}, {lat|-3:3}),0)        ;atl3 index
  atl3_lag1_can =  canonical_diff1(month|5:7,{lon|-20:0}, {lat|-3:3})        ;atl3 index
 ; atl3_lag1_can =  canonical_diff1(month|5:7,lon|:, lat|:)        ;atl3 index
  printVarSummary(atl3_lag1_can)
;----------------------------------------------------------------------------------------------------------------------
;==========Non-canonical
;==============================================================
        raindP1_cases_map1  = new((/n_atlP,dimZ(1),dimZ(2),dimZ(3)/),float,sst_1_s@_FillValue)
        printVarSummary(raindP1_cases_map1)
    
        raindP1_cases_maps1 = sst_1_s(valsP,:,:,:)
        printVarSummary(raindP1_cases_maps1)
        
        rainb1p1 = (raindP1_cases_maps1(time|0,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1p)    
  
        rainb2p1 = (raindP1_cases_maps1(time|1,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2p1)    
  
        rainb3p1 = (raindP1_cases_maps1(time|5,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3p1)    
  
        rainb4p1 = (raindP1_cases_maps1(time|8,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4p1) 
         
       rain_fp1 = (/rainb1p1,rainb2p1,rainb3p1,rainb4p1/)       
 
       copy_VarCoords_n(sst_1_s,rain_fp1,(/1,2,3/))
       rain_fp1!0 = "time"    
       printVarSummary(rain_fp1) 
;===================================================================================================================     
        raindN1_cases_maps1  = new((/n_atlN,dimZ(1),dimZ(2),dimZ(3)/),float,sst_1_s@_FillValue)
        printVarSummary(raindN1_cases_maps1)
    
        raindN1_cases_maps1 = sst_1_s(valsN,:,:,:)
        printVarSummary(raindN1_cases_maps1)
        
        rainb1n1 = (raindN1_cases_maps1(time|3,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1n1)    
  
        rainb2n1 = (raindN1_cases_maps1(time|4,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2n1)    
  
        rainb3n1 = (raindN1_cases_maps1(time|6,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3n1)    
  
        rainb4n1 = (raindN1_cases_maps1(time|7,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4n1)  
      
       rain_fn1 = (/rainb1n1,rainb2n1,rainb3n1,rainb4n1/)       
       copy_VarCoords_n(sst_1_s,rain_fn1,(/1,2,3/))
       rain_fn1!0 = "time"    
       printVarSummary(rain_fn1) 
;==============================================================
  non_canonical_diff1 = (dim_avg_n_Wrap(rain_fp1,0) - dim_avg_n_Wrap(rain_fn1,0))/(2.)
  copy_VarCoords(sst_1_s(0,:,:,:),non_canonical_diff1)                             ; copy coordinate arrays 
  printVarSummary(non_canonical_diff1)
;----------------------------------------------------------------------------------------------------------------------
 ; atl3_lag1_non =  dim_avg_n_Wrap(non_canonical_diff1(month|5:7,{lon|-20:0}, {lat|-3:3}),0)        ;atl3 index
  atl3_lag1_non =  non_canonical_diff1(month|5:7,{lon|-20:0}, {lat|-3:3})        ;atl3 index
 ; atl3_lag1_non =  non_canonical_diff1(month|5:7,lon|:, lat|:)        ;atl3 index
  printVarSummary(atl3_lag1_non)

;----------------------------------------------------------------------------------------------------------------------
  f      = addfile("/home/gyuk/Downloads/GYUK/GEOMAR/Final_work/sst_ensemble_data_results/regrid_sst2_mean.nc", "r")
  TIME   = f->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  sst_00   = short2flt(f->sst(iYYYY,:,:))
 
  sst_2 = sst_00 - 273.15
  copy_VarCoords(sst_00,sst_2)                             ; copy coordinate arrays
  
  sst_2!1 = "lat"
  sst_2!2 = "lon"
  
  delete(f)		;	Not needed anymore
  delete(TIME)
  delete(YYYY)
  delete(iYYYY)

  printMinMax(sst_2,True)
;==========================
  sst_2Clm      = clmMonTLL(sst_2)
  printVarSummary(sst_2Clm)                   ; (12,nlat,nlon)   
  
  sst_2         = calcMonAnomTLL(sst_2,sst_2Clm)    ; replace with anonamlies
  sst_2@long_name = "ANOMALIES: "+sst_2@long_name
  delete(sst_2Clm)
  ;printMinMax(sst_2,True)
  printVarSummary(sst_1)
;=========================
;----------------------------
  nyears = dimsizes(sst_2&time)/nmos
  ;print(nyears)

  sst_2_s     = onedtond(ndtooned(sst_2),(/nyears,nmos,dimsizes(sst_2&lat),dimsizes(sst_2&lon)/))
  copy_VarCoords(sst_2(0,:,:),sst_2_s(0,0,:,:))                             ; copy coordinate arrays
  sst_2_s!0 = "time"
  sst_2_s&time = fspan(yrStrt,yrLast,nyears)
  sst_2_s!1 = "month"
  ;printVarSummary(sst_2_s)


  sst_2_s = dtrend_msg_n(sst_2_s&time,sst_2_s,True,False,0)        ;removes linear trends
  ;printVarSummary(sst_2_s)
  ;printMinMax(sst_2_s,True)
;==============================================================
;=================== Canonical
;***************************************************************************************************************
  dimZ = dimsizes(sst_2_s)
  printVarSummary(dimZ)
        raindP2_cases_maps2  = new((/n_atlP,dimZ(1),dimZ(2),dimZ(3)/),float,sst_2_s@_FillValue)
     ;   printVarSummary(raindP2_cases_maps2)
    
        raindP2_cases_maps2 = sst_2_s(valsP,:,:,:)
     ;   printVarSummary(raindP2_cases_maps2)
;==============================================================
;   Positive Canonical event years
;==============================================================                    
        rainb1p2 = (raindP2_cases_maps2(time|2,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb1p2)    
  
        rainb2p2 = (raindP2_cases_maps2(time|3,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb2p2)    
  
        rainb3p2 = (raindP2_cases_maps2(time|4,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb3p2)    
  
       rainb4p2 = (raindP2_cases_maps2(time|6,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb4p2) 
  
        rainb5p2 = (raindP2_cases_maps2(time|7,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb5p2) 
         
       rain_fp2 = (/rainb1p2,rainb2p2,rainb3p2,rainb4p2,rainb5p2/)       
 
       copy_VarCoords_n(sst_2_s,rain_fp2,(/1,2,3/))
       rain_fp2!0 = "time"    
       printVarSummary(rain_fp2)   
;==============================================================
;   Negative Canonical event years
;==============================================================              
        raindN2_cases_maps2  = new((/n_atlN,dimZ(1),dimZ(2),dimZ(3)/),float,sst_2_s@_FillValue)
        printVarSummary(raindN2_cases_maps2)
    
        raindN2_cases_maps2 = sst_2_s(valsN,:,:,:)
      ;  printVarSummary(raindN2_cases_maps2)
        
        rainb1n2 = (raindN2_cases_maps2(time|0,month|:,lat|:,lon|:))	;now based on the basin average
     ;   printVarSummary(rainb1n2)    
  
        rainb2n2 = (raindN2_cases_maps2(time|1,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb2n2)    
  
        rainb3n2 = (raindN2_cases_maps2(time|2,month|:,lat|:,lon|:))	;now based on the basin average
       ; printVarSummary(rainb3n2)    
  
        rainb4n2 = (raindN2_cases_maps2(time|5,month|:,lat|:,lon|:))	;now based on the basin average
       ; printVarSummary(rainb4n2) 
       
       rain_fn2 = (/rainb1n2,rainb2n2,rainb3n2,rainb4n2/)       
       copy_VarCoords_n(sst_2_s,rain_fn2,(/1,2,3/))
       rain_fn2!0 = "time"   
       printVarSummary(rain_fn2) 
;==============================================================
  canonical_diff12 = (dim_avg_n_Wrap(rain_fp2,0) - dim_avg_n_Wrap(rain_fn2,0))/(2.)                       ; Difference
  copy_VarCoords(sst_2_s(0,:,:,:),canonical_diff12)                             ; copy coordinate arrays 
  printVarSummary(canonical_diff12)
;---------------------------------------------------------------------------
 ; atl3_lag2_can =  dim_avg_n_Wrap(canonical_diff12(month|5:7,{lon|-20:0}, {lat|-3:3}),0)        ;atl3 index
  atl3_lag2_can =  canonical_diff12(month|5:7,{lon|-20:0}, {lat|-3:3})        ;atl3 index
 ; atl3_lag2_can =  canonical_diff12(month|5:7,lon|:, lat|:)        ;atl3 index

;----------------------------------------------------------------------------------------------------------------------
;==============================================================
;==========Non-canonical
;==============================================================
  dimZ = dimsizes(sst_2_s)
  printVarSummary(dimZ)
      
        raindP12_cases_map12  = new((/n_atlP,dimZ(1),dimZ(2),dimZ(3)/),float,sst_2_s@_FillValue)
        printVarSummary(raindP12_cases_map12)
    
        raindP12_cases_map12 = sst_2_s(valsP,:,:,:)
        printVarSummary(raindP12_cases_map12)
        
        rainb1p12 = (raindP12_cases_map12(time|0,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1p2)    
  
        rainb2p12 = (raindP12_cases_map12(time|1,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2p12)    
  
        rainb3p12 = (raindP12_cases_map12(time|5,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3p12)    
  
        rainb4p12 = (raindP12_cases_map12(time|8,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4p12) 
         
       rain_fp12 = (/rainb1p12,rainb2p12,rainb3p12,rainb4p12/)        
       copy_VarCoords_n(sst_2_s,rain_fp12,(/1,2,3/))
       rain_fp12!0 = "time"    
       printVarSummary(rain_fp12) 
;===================================================================
        raindN12_cases_maps12  = new((/n_atlN,dimZ(1),dimZ(2),dimZ(3)/),float,sst_2_s@_FillValue)
        printVarSummary(raindN12_cases_maps12)
    
        raindN12_cases_maps12 = sst_2_s(valsN,:,:,:)
        printVarSummary(raindN12_cases_maps12)
        
        rainb1n12 = (raindN12_cases_maps12(time|3,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1n12)    
  
        rainb2n12 = (raindN12_cases_maps12(time|4,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2n12)    
  
        rainb3n12 = (raindN12_cases_maps12(time|6,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3n12)    
  
        rainb4n12 = (raindN12_cases_maps12(time|7,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4n12) 
       
       rain_fn12 = (/rainb1n12,rainb2n12,rainb3n12,rainb4n12/)       
       copy_VarCoords_n(sst_2_s,rain_fn12,(/1,2,3/))
       rain_fn12!0 = "time"   
       printVarSummary(rain_fn12) 
;==============================================================
  non_canonical_diff12 = (dim_avg_n_Wrap(rain_fp12,0) - dim_avg_n_Wrap(rain_fn12,0))/(2.)
  copy_VarCoords(sst_2_s(0,:,:,:),non_canonical_diff12)                             ; copy coordinate arrays 
  printVarSummary(non_canonical_diff12)
  printMinMax(non_canonical_diff12,True)
;----------------------------------------------------------------------------------------------------------------------
 ; atl3_lag2_non =  dim_avg_n_Wrap(non_canonical_diff12(month|5:7,{lon|-20:0}, {lat|-3:3}),0)        ;atl3 index
  atl3_lag2_non =  non_canonical_diff12(month|5:7,{lon|-20:0}, {lat|-3:3})        ;atl3 index
 ; atl3_lag2_non =  non_canonical_diff12(month|5:7,lon|:, lat|:)        ;atl3 index
  printVarSummary(atl3_lag2_non)

;----------------------------------------------------------------------------------------------------------------------
  f      = addfile("/home/gyuk/Downloads/GYUK/GEOMAR/Final_work/sst_ensemble_data_results/regrid_sst3_mean.nc", "r")
  TIME   = f->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  sst_000   = short2flt(f->sst(iYYYY,:,:))
  printVarSummary(sst_000)       ; variable overview
  
  sst_3 = sst_000 - 273.15
  copy_VarCoords(sst_000,sst_3)                             ; copy coordinate arrays
  sst_3!1 = "lat"
  sst_3!2 = "lon"
  printMinMax(sst_3,True)
  delete(f)		;	Not needed anymore
  delete(TIME)
  delete(YYYY)
  delete(iYYYY)
;==========================
  sst_3Clm      = clmMonTLL(sst_3)
  printVarSummary(sst_3Clm)                   ; (12,nlat,nlon)   
  
  sst_3         = calcMonAnomTLL(sst_3,sst_3Clm)    ; replace with anonamlies
  sst_3@long_name = "ANOMALIES: "+sst_3@long_name
  delete(sst_3Clm)
  ;printMinMax(sst_3,True)
  printVarSummary(sst_3)
;=========================
;----------------------------
  nyears = dimsizes(sst_3&time)/nmos
  ;print(nyears)

  sst_3_s     = onedtond(ndtooned(sst_3),(/nyears,nmos,dimsizes(sst_3&lat),dimsizes(sst_3&lon)/))
  copy_VarCoords(sst_3(0,:,:),sst_3_s(0,0,:,:))                             ; copy coordinate arrays
  sst_3_s!0 = "time"
  sst_3_s&time = fspan(yrStrt,yrLast,nyears)
  sst_3_s!1 = "month"
  ;printVarSummary(sst_3_s)

  sst_3_s = dtrend_msg_n(sst_3_s&time,sst_3_s,True,False,0)        ;removes linear trends
  ;printVarSummary(sst_3_s)
  ;printMinMax(sst_3_s,True)
;==============================================================
;=================== Canonical
  dimZ = dimsizes(sst_3_s)
  printVarSummary(dimZ)
        raindP13_cases_maps3  = new((/n_atlP,dimZ(1),dimZ(2),dimZ(3)/),float,sst_3_s@_FillValue)
     ;   printVarSummary(raindP13_cases_maps3)
    
        raindP13_cases_maps3 = sst_3_s(valsP,:,:,:)
     ;   printVarSummary(raindP13_cases_maps3)
;==============================================================
;   Positive Canonical event years
;==============================================================                    
        rainb1p3 = (raindP13_cases_maps3(time|2,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb1p3)    
  
        rainb2p3 = (raindP13_cases_maps3(time|3,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb2p3)    
  
        rainb3p3 = (raindP13_cases_maps3(time|4,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb3p3)    
  
       rainb4p3 = (raindP13_cases_maps3(time|6,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb4p3) 
  
        rainb5p3 = (raindP13_cases_maps3(time|7,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb5p3) 
         
       rain_fp3 = (/rainb1p3,rainb2p3,rainb3p3,rainb4p3,rainb5p3/)        
       copy_VarCoords_n(sst_3_s,rain_fp3,(/1,2,3/))
       rain_fp3!0 = "time"   
       printVarSummary(rain_fp3)  
;==============================================================
;   Negative Canonical event years
;==============================================================              
        raindN13_cases_maps3  = new((/n_atlN,dimZ(1),dimZ(2),dimZ(3)/),float,sst_3_s@_FillValue)
        printVarSummary(raindN13_cases_maps3)
    
        raindN13_cases_maps3 = sst_3_s(valsN,:,:,:)
      ;  printVarSummary(raindN13_cases_maps3)
        
        rainb1n3 = (raindN13_cases_maps3(time|0,month|:,lat|:,lon|:))	;now based on the basin average
     ;   printVarSummary(rainb1n3)     
        rainb2n3 = (raindN13_cases_maps3(time|1,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb2n3)     
        rainb3n3 = (raindN13_cases_maps3(time|2,month|:,lat|:,lon|:))	;now based on the basin average
       ; printVarSummary(rainb3n3)     
        rainb4n3 = (raindN13_cases_maps3(time|5,month|:,lat|:,lon|:))	;now based on the basin average
       ; printVarSummary(rainb4n3) 
            
       rain_fn3 = (/rainb1n3,rainb2n3,rainb3n3,rainb4n3/)       
       copy_VarCoords_n(sst_3_s,rain_fn3,(/1,2,3/))
       rain_fn3!0 = "time"    
       printVarSummary(rain_fn3) 
;==============================================================
  canonical_diff13 = (dim_avg_n_Wrap(rain_fp3,0) - dim_avg_n_Wrap(rain_fn3,0))/(2.)                       ; Difference
  copy_VarCoords(sst_3_s(0,:,:,:),canonical_diff13)                             ; copy coordinate arrays 
  printVarSummary(canonical_diff13)
;---------------------------------------------------------------------------
 ; atl3_lag3_can =  dim_avg_n_Wrap(canonical_diff13(month|5:7,{lon|-20:0}, {lat|-3:3}),0)        ;atl3 index
  atl3_lag3_can =  canonical_diff13(month|5:7,{lon|-20:0}, {lat|-3:3})       ;atl3 index
 ; atl3_lag3_can =  canonical_diff13(month|5:7,lon|:, lat|:)       ;atl3 index

;==============================================================
;==========Non-canonical
;==============================================================
        raindP31_cases_map13  = new((/n_atlP,dimZ(1),dimZ(2),dimZ(3)/),float,sst_3_s@_FillValue)
        printVarSummary(raindP31_cases_map13)
    
        raindP31_cases_map13 = sst_3_s(valsP,:,:,:)
        printVarSummary(raindP31_cases_map13)
   
       rainb1p13 = (raindP31_cases_map13(time|0,month|:,lat|:,lon|:))	;now based on the basin average
       printVarSummary(rainb1p13)    
  
       rainb2p13 = (raindP31_cases_map13(time|1,month|:,lat|:,lon|:))	;now based on the basin average
       printVarSummary(rainb2p13)    
  
       rainb3p13 = (raindP31_cases_map13(time|5,month|:,lat|:,lon|:))	;now based on the basin average
       printVarSummary(rainb3p13)    
  
       rainb4p13 = (raindP31_cases_map13(time|8,month|:,lat|:,lon|:))	;now based on the basin average
       printVarSummary(rainb4p13) 
       
       rain_fp13 = (/rainb1p13,rainb2p13,rainb3p13,rainb4p13/)        
       copy_VarCoords_n(sst_3_s,rain_fp13,(/1,2,3/))
       rain_fp13!0 = "time"   
       printVarSummary(rain_fp13) 
;===================================================================
 ; rain_fp13 = (/raindP31_cases_map13(0,:,:,:),raindP31_cases_map13(1,:,:,:),raindP31_cases_map13(5,:,:,:),raindP31_cases_map13(8,:,:,:)/) 
   
       ; copy_VarCoords_n(sst_3_s,rain_fp13,(/1,2,3/))
       ;rain_fp13!0 = "time"
      ;  printVarSummary(rain_fp13)      
;==============================================================
;==============================================================      
        raindN31_cases_maps13  = new((/n_atlN,dimZ(1),dimZ(2),dimZ(3)/),float,sst_3_s@_FillValue)
        printVarSummary(raindN31_cases_maps13)
    
        raindN31_cases_maps13 = sst_3_s(valsN,:,:,:)
        printVarSummary(raindN31_cases_maps13)
        
        rainb1n13 = (raindN31_cases_maps13(time|3,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1n13)    
  
        rainb2n13 = (raindN31_cases_maps13(time|4,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2n13)    
  
        rainb3n13 = (raindN31_cases_maps13(time|6,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3n13)    
  
        rainb4n13 = (raindN31_cases_maps13(time|7,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4n13)  
      
       rain_fn13 = (/rainb1n13,rainb2n13,rainb3n13,rainb4n13/)       
       copy_VarCoords_n(sst_3_s,rain_fn13,(/1,2,3/))
       rain_fn13!0 = "time"    
       printVarSummary(rain_fn13) 
;==============================================================
  non_canonical_diff13 = (dim_avg_n_Wrap(rain_fp13,0) - dim_avg_n_Wrap(rain_fn13,0))/(2.)
  copy_VarCoords(sst_3_s(0,:,:,:),non_canonical_diff13)                             ; copy coordinate arrays 
  printVarSummary(non_canonical_diff13)
;----------------------------------------------------------------------------------------------------------------------
 ; atl3_lag3_non =  dim_avg_n_Wrap(non_canonical_diff13(month|5:7,{lon|-20:0}, {lat|-3:3}),0)        ;atl3 index
  atl3_lag3_non =  non_canonical_diff13(month|5:7,{lon|-20:0}, {lat|-3:3})        ;atl3 index
 ; atl3_lag3_non =  non_canonical_diff13(month|5:7,lon|:, lat|:)        ;atl3 index
  printVarSummary(atl3_lag3_non)

;----------------------------------------------------------------------------------------------------------------------
  f      = addfile("/home/gyuk/Downloads/GYUK/GEOMAR/Final_work/sst_ensemble_data_results/regrid_sst4_mean.nc", "r")
  TIME   = f->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  sst_0000   = short2flt(f->sst(iYYYY,:,:))
  printVarSummary(sst_0000)       ; variable overview
  
  sst_4 = sst_0000 - 273.15
  copy_VarCoords(sst_0000,sst_4)                             ; copy coordinate arrays
  sst_4!1 = "lat"
  sst_4!2 = "lon"
  printMinMax(sst_4,True)
  delete(f)		;	Not needed anymore
  delete(TIME)
  delete(YYYY)
  delete(iYYYY)
;========================================================================================================================
  sst_4Clm      = clmMonTLL(sst_4)
  printVarSummary(sst_4Clm)                   ; (12,nlat,nlon)   
  
  sst_4         = calcMonAnomTLL(sst_4,sst_4Clm)    ; replace with anonamlies
  sst_4@long_name = "ANOMALIES: "+sst_4@long_name
  delete(sst_4Clm)
  ;printMinMax(sst_4,True)
  printVarSummary(sst_4)
;=========================
;----------------------------
  nyears = dimsizes(sst_4&time)/nmos
  ;print(nyears)

  sst_4_s     = onedtond(ndtooned(sst_4),(/nyears,nmos,dimsizes(sst_4&lat),dimsizes(sst_4&lon)/))
  copy_VarCoords(sst_4(0,:,:),sst_4_s(0,0,:,:))                             ; copy coordinate arrays
  sst_4_s!0 = "time"
  sst_4_s&time = fspan(yrStrt,yrLast,nyears)
  sst_4_s!1 = "month"
  ;printVarSummary(sst_4_s)

  sst_4_s = dtrend_msg_n(sst_4_s&time,sst_4_s,True,False,0)        ;removes linear trends
  ;printVarSummary(sst_4_s)
;==============================================================
;=================== Canonical
;==============================================================
  dimZ = dimsizes(sst_4_s)
  printVarSummary(dimZ)
        raindP14_cases_maps4  = new((/n_atlP,dimZ(1),dimZ(2),dimZ(3)/),float,sst_4_s@_FillValue)
     ;   printVarSummary(raindP_cases_maps4)
    
        raindP14_cases_maps4 = sst_4_s(valsP,:,:,:)
     ;   printVarSummary(raindP14_cases_maps4)
;==============================================================
;   Positive Canonical event years
;==============================================================                    
        rainb1p4 = (raindP14_cases_maps4(time|2,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb1p4)    
  
        rainb2p4 = (raindP14_cases_maps4(time|3,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb2p4)    
  
        rainb3p4 = (raindP14_cases_maps4(time|4,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb3p4)    
  
       rainb4p4 = (raindP14_cases_maps4(time|6,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb4p4) 
  
        rainb5p4 = (raindP14_cases_maps4(time|7,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb5p4) 
         
       rain_fp4 = (/rainb1p4,rainb2p4,rainb3p4,rainb4p4,rainb5p4/)        
       copy_VarCoords_n(sst_4_s,rain_fp4,(/1,2,3/))
       rain_fp4!0 = "time"   
       printVarSummary(rain_fp4)   
;==============================================================
;   Negative Canonical event years
;==============================================================             
        raindN14_cases_maps4  = new((/n_atlN,dimZ(1),dimZ(2),dimZ(3)/),float,sst_4_s@_FillValue)
        printVarSummary(raindN14_cases_maps4)
    
        raindN14_cases_maps4 = sst_4_s(valsN,:,:,:)
      ;  printVarSummary(raindN14_cases_maps4)
        
        rainb1n4 = (raindN14_cases_maps4(time|0,month|:,lat|:,lon|:))	;now based on the basin average
     ;   printVarSummary(rainb1n4)    
  
        rainb2n4 = (raindN14_cases_maps4(time|1,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb2n4)    
  
        rainb3n4 = (raindN14_cases_maps4(time|2,month|:,lat|:,lon|:))	;now based on the basin average
       ; printVarSummary(rainb3n4)    
  
        rainb4n4 = (raindN14_cases_maps4(time|5,month|:,lat|:,lon|:))	;now based on the basin average
       ; printVarSummary(rainb4n4) 
        
       rain_fn4 = (/rainb1n4,rainb2n4,rainb3n4,rainb4n4/)        
       copy_VarCoords_n(sst_4_s,rain_fn4,(/1,2,3/))
       rain_fn4!0 = "time"    
       printVarSummary(rain_fn4) 
;==============================================================
  canonical_diff14 = (dim_avg_n_Wrap(rain_fp4,0) - dim_avg_n_Wrap(rain_fn4,0))/(2.)                       ; Difference
  copy_VarCoords(sst_1_s(0,:,:,:),canonical_diff14)                             ; copy coordinate arrays 
  printVarSummary(canonical_diff14)
  printMinMax(canonical_diff14,True)
;---------------------------------------------------------------------------
 ; atl3_lag4_can =  dim_avg_n_Wrap(canonical_diff14(month|5:7,{lon|-20:0}, {lat|-3:3}),0)        ;atl3 index
  atl3_lag4_can =  canonical_diff14(month|5:7,{lon|-20:0}, {lat|-3:3})        ;atl3 index
 ; atl3_lag4_can =  canonical_diff14(month|5:7,lon|:, lat|:)        ;atl3 index

;----------------------------------------------------------------------------------------------------------------------
;==========Non-canonical
;==============================================================
        raindP41_cases_map14  = new((/n_atlP,dimZ(1),dimZ(2),dimZ(3)/),float,sst_4_s@_FillValue)
        printVarSummary(raindP41_cases_map14)
    
        raindP41_cases_map14 = sst_4_s(valsP,:,:,:)
        printVarSummary(raindP41_cases_map14)
        
        rainb1p14 = (raindP41_cases_map14(time|0,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1p14)    
  
        rainb2p14 = (raindP41_cases_map14(time|1,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2p14)    
  
        rainb3p14 = (raindP41_cases_map14(time|5,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3p14)    
  
        rainb4p14 = (raindP41_cases_map14(time|8,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4p14) 
        
       rain_fp14 = (/rainb1p14,rainb2p14,rainb3p14,rainb4p14/)        
       copy_VarCoords_n(sst_4_s,rain_fp14,(/1,2,3/))
       rain_fp14!0 = "time"
       printVarSummary(rain_fp14) 
;===================================================================
        raindN41_cases_maps14  = new((/n_atlN,dimZ(1),dimZ(2),dimZ(3)/),float,sst_4_s@_FillValue)
        printVarSummary(raindN41_cases_maps14)
    
        raindN41_cases_maps14 = sst_4_s(valsN,:,:,:)
        printVarSummary(raindN41_cases_maps14)
        
        rainb1n14 = (raindN41_cases_maps14(time|3,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1n14)    
  
        rainb2n14 = (raindN41_cases_maps14(time|4,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2n14)    
  
        rainb3n14 = (raindN41_cases_maps14(time|6,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3n14)    
  
        rainb4n14 = (raindN41_cases_maps14(time|7,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4n14) 
        
       rain_fn14 = (/rainb1n14,rainb2n14,rainb3n14,rainb4n14/)       
       copy_VarCoords_n(sst_4_s,rain_fn14,(/1,2,3/))
       rain_fn14!0 = "time"  
       printVarSummary(rain_fn14) 
;==============================================================
  non_canonical_diff14 = (dim_avg_n_Wrap(rain_fp14,0) - dim_avg_n_Wrap(rain_fn14,0))/(2.)
  copy_VarCoords(sst_4_s(0,:,:,:),non_canonical_diff14)                             ; copy coordinate arrays 
  printVarSummary(non_canonical_diff14)
;----------------------------------------------------------------------------------------------------------------------
 ; atl3_lag4_non =  dim_avg_n_Wrap(non_canonical_diff14(month|5:7,{lon|-20:0}, {lat|-3:3}),0)        ;atl3 index
  atl3_lag4_non =  non_canonical_diff14(month|5:7,{lon|-20:0}, {lat|-3:3})        ;atl3 index
 ; atl3_lag4_non =  non_canonical_diff14(month|5:7,lon|:, lat|:)        ;atl3 index
  printVarSummary(atl3_lag4_non)
 
;----------------------------------------------------------------------------------------------------------------------
  f           = addfile("/home/gyuk/Downloads/GYUK/GEOMAR/Final_work/sst_ensemble_data_results/regrid_sst5_mean.nc", "r")
  TIME        = f->time
  YYYY        = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY       = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  sst_00000   = short2flt(f->sst(iYYYY,:,:))
  printVarSummary(sst_00000)       ; variable overview
  
  sst_5 = sst_00000 - 273.15
  copy_VarCoords(sst_00000,sst_5)                             ; copy coordinate arrays
  sst_5!1 = "lat"
  sst_5!2 = "lon"
  printMinMax(sst_5,True)
  delete(f)		;	Not needed anymore
  delete(TIME)
  delete(YYYY)
  delete(iYYYY)

;==========================
  sst_5Clm      = clmMonTLL(sst_5)
  printVarSummary(sst_5Clm)                   ; (12,nlat,nlon)   
  
  sst_5         = calcMonAnomTLL(sst_5,sst_5Clm)    ; replace with anonamlies
  sst_5@long_name = "ANOMALIES: "+sst_5@long_name
  delete(sst_5Clm)
  ;printMinMax(sst_5,True)
  printVarSummary(sst_5)
;=========================
;----------------------------
  nyears = dimsizes(sst_5&time)/nmos
  ;print(nyears)

  sst_5_s     = onedtond(ndtooned(sst_5),(/nyears,nmos,dimsizes(sst_5&lat),dimsizes(sst_5&lon)/))
  copy_VarCoords(sst_5(0,:,:),sst_5_s(0,0,:,:))                             ; copy coordinate arrays
  sst_5_s!0 = "time"
  sst_5_s&time = fspan(yrStrt,yrLast,nyears)
  sst_5_s!1 = "month"
  ;printVarSummary(sst_5_s)

  sst_5_s = dtrend_msg_n(sst_5_s&time,sst_5_s,True,False,0)        ;removes linear trends
  ;printVarSummary(sst_5_s)
;==============================================================
;=================== Canonical
;***************************************************************************************************************
  dimZ = dimsizes(sst_5_s)
  printVarSummary(dimZ)
        raindP15_cases_maps5  = new((/n_atlP,dimZ(1),dimZ(2),dimZ(3)/),float,sst_5_s@_FillValue)
     ;   printVarSummary(raindP15_cases_maps5)
    
        raindP15_cases_maps5 = sst_5_s(valsP,:,:,:)
     ;   printVarSummary(raindP15_cases_maps5)     
;==============================================================
;   Positive Canonical event years
;==============================================================                    
        rainb1p5 = (raindP15_cases_maps5(time|2,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb1p5)    
  
        rainb2p5 = (raindP15_cases_maps5(time|3,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb2p5)    
  
        rainb3p5 = (raindP15_cases_maps5(time|4,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb3p5)    
  
        rainb4p5 = (raindP15_cases_maps5(time|6,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb4p5) 
  
        rainb5p5 = (raindP15_cases_maps5(time|7,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb5p5) 
         
       rain_fp5 = (/rainb1p5,rainb2p5,rainb3p5,rainb4p5,rainb5p5/)       
       copy_VarCoords_n(sst_5_s,rain_fp5,(/1,2,3/))
       rain_fp5!0 = "time"  
       printVarSummary(rain_fp5)  
;==============================================================
;   Negative Canonical event years
;==============================================================             
        raindN15_cases_maps5  = new((/n_atlN,dimZ(1),dimZ(2),dimZ(3)/),float,sst_5_s@_FillValue)
        printVarSummary(raindN15_cases_maps5)
    
        raindN15_cases_maps5 = sst_5_s(valsN,:,:,:)
      ;  printVarSummary(raindN15_cases_maps5)
        
        rainb1n5 = (raindN15_cases_maps5(time|0,month|:,lat|:,lon|:))	;now based on the basin average
     ;   printVarSummary(rainb1n5)    
  
        rainb2n5 = (raindN15_cases_maps5(time|1,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb2n5)    
  
        rainb3n5 = (raindN15_cases_maps5(time|2,month|:,lat|:,lon|:))	;now based on the basin average
       ; printVarSummary(rainb3n5)    
  
        rainb4n5 = (raindN15_cases_maps5(time|5,month|:,lat|:,lon|:))	;now based on the basin average
       ; printVarSummary(rainb4n5) 
        
       rain_fn5 = (/rainb1n5,rainb2n5,rainb3n5,rainb4n5/)       
       copy_VarCoords_n(sst_5_s,rain_fn5,(/1,2,3/))
       rain_fn5!0 = "time"   
       printVarSummary(rain_fn5) 
;==============================================================
  canonical_diff15 = (dim_avg_n_Wrap(rain_fp5,0) - dim_avg_n_Wrap(rain_fn5,0))/(2.)                       ; Difference
  copy_VarCoords(sst_5_s(0,:,:,:),canonical_diff15)                             ; copy coordinate arrays 
  printVarSummary(canonical_diff15)
;---------------------------------------------------------------------------
 ; atl3_lag5_can =  dim_avg_n_Wrap(canonical_diff15(month|5:7,{lon|-20:0}, {lat|-3:3}),0)        ;atl3 index
  atl3_lag5_can =  canonical_diff15(month|5:7,{lon|-20:0}, {lat|-3:3})        ;atl3 index
 ; atl3_lag5_can =  canonical_diff15(month|5:7,lon|:, lat|:)        ;atl3 index
  
;----------------------------------------------------------------------------------------------------------------------
;==========Non-canonical
;==============================================================  
        raindP51_cases_map15  = new((/n_atlP,dimZ(1),dimZ(2),dimZ(3)/),float,sst_5_s@_FillValue)
        printVarSummary(raindP51_cases_map15)
    
        raindP51_cases_map15 = sst_5_s(valsP,:,:,:)
        printVarSummary(raindP51_cases_map15)
        
        rainb1p15 = (raindP51_cases_map15(time|0,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1p15)    
  
        rainb2p15 = (raindP51_cases_map15(time|1,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2p15)    
  
        rainb3p15 = (raindP51_cases_map15(time|5,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3p15)    
  
        rainb4p15 = (raindP51_cases_map15(time|8,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4p15) 
        
       rain_fp15 = (/rainb1p15,rainb2p15,rainb3p15,rainb4p15/)       
       copy_VarCoords_n(sst_5_s,rain_fp15,(/1,2,3/))
       rain_fp15!0 = "time"   
       printVarSummary(rain_fp15) 
;===========================================================================================================================     
        raindN51_cases_maps15  = new((/n_atlN,dimZ(1),dimZ(2),dimZ(3)/),float,sst_5_s@_FillValue)
        printVarSummary(raindN51_cases_maps15)
    
        raindN51_cases_maps15 = sst_5_s(valsN,:,:,:)
        printVarSummary(raindN51_cases_maps15)
        
        rainb1n15 = (raindN51_cases_maps15(time|3,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1n15)    
  
        rainb2n15 = (raindN51_cases_maps15(time|4,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2n15)    
  
        rainb3n15 = (raindN51_cases_maps15(time|6,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3n15)    
  
        rainb4n15 = (raindN51_cases_maps15(time|7,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4n15) 
  
       rain_fn15 = (/rainb1n15,rainb2n15,rainb3n15,rainb4n15/)       
       copy_VarCoords_n(sst_5_s,rain_fn15,(/1,2,3/))
       rain_fn15!0 = "time"    
       printVarSummary(rain_fn15) 
;==============================================================
  non_canonical_diff15 = (dim_avg_n_Wrap(rain_fp15,0) - dim_avg_n_Wrap(rain_fn15,0))/(2.)
  copy_VarCoords(sst_5_s(0,:,:,:),non_canonical_diff15)                             ; copy coordinate arrays 
  printVarSummary(non_canonical_diff15)
;----------------------------------------------------------------------------------------------------------------------
 ; atl3_lag5_non =  dim_avg_n_Wrap(non_canonical_diff15(month|5:7,{lon|-20:0}, {lat|-3:3}),0)        ;atl3 index
  atl3_lag5_non =  non_canonical_diff15(month|5:7,{lon|-20:0}, {lat|-3:3})        ;atl3 index
 ; atl3_lag5_non =  non_canonical_diff15(month|5:7,lon|:, lat|:)        ;atl3 index
  printVarSummary(atl3_lag5_non)
 
;----------------------------------------------------------------------------------------------------------------------
  f      = addfile("/home/gyuk/Downloads/GYUK/GEOMAR/Final_work/sst_ensemble_data_results/regrid_sst6_mean.nc", "r")
  TIME   = f->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  sst_000000   = short2flt(f->sst(iYYYY,:,:))
  printVarSummary(sst_000000)       ; variable overview
  
  sst_6 = sst_000000 - 273.15
  copy_VarCoords(sst_000000,sst_6)                             ; copy coordinate arrays
  sst_6!1 = "lat"
  sst_6!2 = "lon"
  printMinMax(sst_6,True)
  delete(f)		;	Not needed anymore
  delete(TIME)
  delete(YYYY)
  delete(iYYYY)

;==================================================================================================================
  sst_6Clm      = clmMonTLL(sst_6)
  printVarSummary(sst_6Clm)                   ; (12,nlat,nlon)   
  
  sst_6         = calcMonAnomTLL(sst_6,sst_6Clm)    ; replace with anonamlies
  sst_6@long_name = "ANOMALIES: "+sst_6@long_name
  delete(sst_6Clm)
  ;printMinMax(sst_6,True)
  printVarSummary(sst_6)
;=========================
;----------------------------
  nyears = dimsizes(sst_6&time)/nmos
  ;print(nyears)
  
  sst_6_s     = onedtond(ndtooned(sst_6),(/nyears,nmos,dimsizes(sst_6&lat),dimsizes(sst_6&lon)/))
  copy_VarCoords(sst_6(0,:,:),sst_6_s(0,0,:,:))                             ; copy coordinate arrays
  sst_6_s!0 = "time"
  sst_6_s&time = fspan(yrStrt,yrLast,nyears)
  sst_6_s!1 = "month"
  ;printVarSummary(sst_6_s)

  sst_6_s = dtrend_msg_n(sst_6_s&time,sst_6_s,True,False,0)        ;removes linear trends
  ;printVarSummary(sst_6_s)
  ;printMinMax(sst_6_s,True)
;==============================================================
;=================== Canonical
;***************************************************************************************************************
  dimZ = dimsizes(sst_6_s)
  printVarSummary(dimZ)
        raindP16_cases_maps6  = new((/n_atlP,dimZ(1),dimZ(2),dimZ(3)/),float,sst_6_s@_FillValue)
     ;   printVarSummary(raindP16_cases_maps6)
    
        raindP16_cases_maps6 = sst_6_s(valsP,:,:,:)
     ;   printVarSummary(raindP16_cases_maps6)     
;==============================================================
;   Positive Canonical event years
;==============================================================                   
        rainb1p6 = (raindP16_cases_maps6(time|2,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb1p6)    
  
        rainb2p6 = (raindP16_cases_maps6(time|3,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb2p6)    
  
        rainb3p6 = (raindP16_cases_maps6(time|4,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb3p6)    
  
       rainb4p6 = (raindP16_cases_maps6(time|6,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb4p6) 
  
        rainb5p6 = (raindP16_cases_maps6(time|7,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb5p6) 
         
       rain_fp6 = (/rainb1p6,rainb2p6,rainb3p6,rainb4p6,rainb5p6/)        
       copy_VarCoords_n(sst_6_s,rain_fp6,(/1,2,3/))
       rain_fp6!0 = "time"    
       printVarSummary(rain_fp6)  
;==============================================================
;   Negative Canonical event years
;==============================================================             
        raindN16_cases_maps6  = new((/n_atlN,dimZ(1),dimZ(2),dimZ(3)/),float,sst_6_s@_FillValue)
        printVarSummary(raindN16_cases_maps6)
    
        raindN16_cases_maps6 = sst_6_s(valsN,:,:,:)
      ;  printVarSummary(raindN16_cases_maps6)
        
        rainb1n6 = (raindN16_cases_maps6(time|0,month|:,lat|:,lon|:))	;now based on the basin average
     ;   printVarSummary(rainb1n6)    
  
        rainb2n6 = (raindN16_cases_maps6(time|1,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(rainb2n6)    
  
        rainb3n6 = (raindN16_cases_maps6(time|2,month|:,lat|:,lon|:))	;now based on the basin average
       ; printVarSummary(rainb3n6)    
  
        rainb4n6 = (raindN16_cases_maps6(time|5,month|:,lat|:,lon|:))	;now based on the basin average
       ; printVarSummary(rainb4n6) 
        
       rain_fn6 = (/rainb1n6,rainb2n6,rainb3n6,rainb4n6/)        
       copy_VarCoords_n(sst_6_s,rain_fn6,(/1,2,3/))
       rain_fn6!0 = "time"    
       printVarSummary(rain_fn6) 
;==============================================================
  canonical_diff16 = (dim_avg_n_Wrap(rain_fp6,0) - dim_avg_n_Wrap(rain_fn6,0))/(2.)                       ; Difference
  copy_VarCoords(sst_6_s(0,:,:,:),canonical_diff16)                             ; copy coordinate arrays 
  printVarSummary(canonical_diff16)
;---------------------------------------------------------------------------
 ; atl3_lag6_can =  dim_avg_n_Wrap(canonical_diff16(month|5:7,{lon|-20:0}, {lat|-3:3}),0)        ;atl3 index
  atl3_lag6_can =  canonical_diff16(month|5:7,{lon|-20:0}, {lat|-3:3})        ;atl3 index
 ; atl3_lag6_can =  canonical_diff16(month|5:7,lon|:, lat|:)        ;atl3 index

;----------------------------------------------------------------------------------------------------------------------
;==========Non-canonical
;==============================================================  
        raindP61_cases_map16  = new((/n_atlP,dimZ(1),dimZ(2),dimZ(3)/),float,sst_6_s@_FillValue)
        printVarSummary(raindP61_cases_map16)
    
        raindP61_cases_map16 = sst_6_s(valsP,:,:,:)
        printVarSummary(raindP61_cases_map16)
        
        rainb1p16 = (raindP61_cases_map16(time|0,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1p6)    
  
        rainb2p16 = (raindP61_cases_map16(time|1,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2p16)    
  
        rainb3p16 = (raindP61_cases_map16(time|5,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3p1)    
  
        rainb4p16 = (raindP61_cases_map16(time|8,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4p16) 
  
       
       rain_fp16 = (/rainb1p16,rainb2p16,rainb3p16,rainb4p16/)        
       copy_VarCoords_n(sst_6_s,rain_fp16,(/1,2,3/))
       rain_fp16!0 = "time"    
       printVarSummary(rain_fp16) 
;==============================================================     
        raindN61_cases_maps16  = new((/n_atlN,dimZ(1),dimZ(2),dimZ(3)/),float,sst_6_s@_FillValue)
        printVarSummary(raindN61_cases_maps16)
    
        raindN61_cases_maps16 = sst_6_s(valsN,:,:,:)
        printVarSummary(raindN61_cases_maps16)
        
        rainb1n16 = (raindN61_cases_maps16(time|3,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1n16)    
  
        rainb2n16 = (raindN61_cases_maps16(time|4,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2n16)    
  
        rainb3n16 = (raindN61_cases_maps16(time|6,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3n16)    
  
        rainb4n16 = (raindN61_cases_maps16(time|7,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4n16) 
        
       rain_fn16 = (/rainb1n16,rainb2n16,rainb3n16,rainb4n16/)        
       copy_VarCoords_n(sst_6_s,rain_fn16,(/1,2,3/))
       rain_fn16!0 = "time"    
       printVarSummary(rain_fn16) 
;==============================================================
  non_canonical_diff16 = (dim_avg_n_Wrap(rain_fp16,0) - dim_avg_n_Wrap(rain_fn16,0))/(2.)
  copy_VarCoords(sst_6_s(0,:,:,:),non_canonical_diff16)                             ; copy coordinate arrays 
  printVarSummary(non_canonical_diff16)
;----------------------------------------------------------------------------------------------------------------------
 ; atl3_lag6_non =  dim_avg_n_Wrap(non_canonical_diff16(month|5:7,{lon|-20:0}, {lat|-3:3}),0)       ;atl3 index
  atl3_lag6_non =  non_canonical_diff16(month|5:7,{lon|-20:0}, {lat|-3:3})       ;atl3 index
 ; atl3_lag6_non =  non_canonical_diff16(month|5:7,lon|:, lat|:)       ;atl3 index
  printVarSummary(atl3_lag6_non)
  
  opt = 1
;---------------------------------------------------------------------------------------------------------------------- 
  Rmse1 = wgt_arearmse2(atl3_can, atl3_lag1_can, wgt, opt) 
  Rmse2 = wgt_arearmse2(atl3_can, atl3_lag2_can, wgt, opt) 
  Rmse3 = wgt_arearmse2(atl3_can, atl3_lag3_can, wgt, opt) 
  Rmse4 = wgt_arearmse2(atl3_can, atl3_lag4_can, wgt, opt) 
  Rmse5 = wgt_arearmse2(atl3_can, atl3_lag5_can, wgt, opt) 
  Rmse6 = wgt_arearmse2(atl3_can, atl3_lag6_can, wgt, opt) 
  
  rmse1 = dim_avg_n_Wrap(Rmse1,0)
  rmse2 = dim_avg_n_Wrap(Rmse2,0)
  rmse3 = dim_avg_n_Wrap(Rmse3,0)
  rmse4 = dim_avg_n_Wrap(Rmse4,0)
  rmse5 = dim_avg_n_Wrap(Rmse5,0)
  rmse6 = dim_avg_n_Wrap(Rmse6,0)

;----------------------------------------------------------------------------------------------------------------------
  RmseN1 = wgt_arearmse2(atl3_non, atl3_lag1_non, wgt, opt) 
  RmseN2 = wgt_arearmse2(atl3_non, atl3_lag2_non, wgt, opt) 
  RmseN3 = wgt_arearmse2(atl3_non, atl3_lag3_non, wgt, opt) 
  RmseN4 = wgt_arearmse2(atl3_non, atl3_lag4_non, wgt, opt)   
  RmseN5 = wgt_arearmse2(atl3_non, atl3_lag5_non, wgt, opt)   
  RmseN6 = wgt_arearmse2(atl3_non, atl3_lag6_non, wgt, opt)   
  
  rmse1N = dim_avg_n_Wrap(RmseN1,0)
  rmse2N = dim_avg_n_Wrap(RmseN2,0)
  rmse3N = dim_avg_n_Wrap(RmseN3,0)
  rmse4N = dim_avg_n_Wrap(RmseN4,0)
  rmse5N = dim_avg_n_Wrap(RmseN5,0)
  rmse6N = dim_avg_n_Wrap(RmseN6,0) 

  yAa = (/rmse1,rmse2,rmse3,rmse4,rmse5,rmse6/)
  yAa1 = (/rmse1N,rmse2N,rmse3N,rmse4N,rmse5N,rmse6N/)
  printVarSummary(yAa)
  printVarSummary(yAa1)
  zyAa = new((/6,2/),"float")
     do hh = 0,5
        zyAa(hh,0) = (/yAa(hh::6) /)
        zyAa(hh,1) = (/yAa1(hh::6) /)
     end do
     printVarSummary(zyAa)   
     
;***************************************************************************************************************          
;------------------------------------------------------------------     
;***************************************************************************************************************
  wks = gsn_open_wks("eps","Figure_4_rmse")         ; send graphics to PNG file
  plot = new(1,graphic)                ; create graphic array ; only needed if paneling

  res = True
  res@vpWidthF = 0.7
  res@vpHeightF = 0.5
  res@vpXF = 0.15

  res@trXMinF               = 0.50                    ; adds space on either end
  res@trXMaxF               = 6.0                    ; of the 1st and last bars
  
 
  res@trYMinF =  0.
  res@trYMaxF =  0.8
 
  res@gsnDraw = False
  res@gsnFrame = False
  
  res@gsnYRefLine = 0.0   ;-- create a reference line
  
  res@gsnXYBarChart = True
  res@gsnXYBarChartBarWidth = 0.15           ; change bar widths
  
  res@tmXTOn = False ; Turn off tickmarks and labels 
  res@tmXBOn = True   ;-- don't draw bottom x-axis
  res@tmXTLabelsOn = True   ;-- draw top x-axis labels
  
 ; res@tmYROn               = False    ; Turn on right tickmarks
 ; res@tmXBMajorOutwardLengthF = 0.0               ; draw tickmarks inward
  res@tmXBMinorOutwardLengthF = 0.0               ; draw minor ticks inward
  res@tmXBMajorLengthF = 0.01  
  res@tmXBMinorLengthF = 0.02  


  res@tmXBMode          = "Explicit"         ; explicit labels
  res@tmXBValues        = (/.78,1.78,2.78,3.78,4.78,5.78/)
  res@tmXBLabels = (/"lag-1","lag-2","lag-3","lag-4","lag-5","lag-6"/)
  res@tmXBLabelFontHeightF = 0.0205
  res@tmXTLabelFontHeightF = 0.0205
  res@tmYLLabelFontHeightF = 0.0225
  res@tiMainFontHeightF = 0.025
  res@tiMainFont = "helvetica"
  
  res@gsnRightString = ""
  res@tiYAxisString = "RMSE skill"
     
 ; res@tiMainString          = "Atlantic Nino"
  res@gsnXYBarChartColors = (/"blue"/)
  plotA1  = gsn_csm_xy(wks,fspan(0.7,5.7,6),zyAa(:,0),res)                  ; create plot
  res@trYMinF = -0.5
  res@trYMaxF = 0.5
  res@gsnXYBarChartColors = (/"red"/)
  plotB1  = gsn_csm_xy(wks,fspan(0.85,5.85,6),zyAa(:,1),res)                  ; create plot
  
  overlay(plotA1,plotB1)
  plot(0)  = plotA1
;**********************************************************
;---Create two individual labelbars
 ; lbres                    = True          ; labelbar only resources
 ; lbres@vpWidthF           = 0.3           ; labelbar width
 ; lbres@vpHeightF          = 0.1           ; labelbar height
 ; lbres@lbBoxMajorExtentF  = 0.36          ; puts space between color boxes
 ; lbres@lbMonoFillPattern  = True          ; Solid fill pattern
 ; lbres@lbFillColors       = (/"firebrick","red"/)
 ; lbres@lbLabelFontHeightF = 0.015         ; font height. default is small
 ; lbres@lbLabelJust        = "CenterLeft"  ; left justify labels
 ; lbres@lbPerimOn          = False
 ; lbres@lgPerimColor       = "white"
	
 ; lbres@lbFillColors    = (/"red","blue"/)
 ; labels                = (/"Non-canonical Events","Canonical events"/)
 ; gsn_labelbar_ndc(wks,2,labels,0.20,0.45,lbres)    ; draw left labelbar column

;===============================================================================================================
;***************************************************************************************************************************************************

  x_0 = 0.
  x_1 = 12.
  y_0 = 0.508795 
  y_1 = 0.508795 
  
  
  
 ;---Resources for the polyline
  pres1                  = True           ; polyline mods desired
  pres1@gsLineThicknessF = 0.8            ; line thickness   
  pres1@gsLineColor      = "black"         ; color of lines
  pres1@gsLineDashPattern = 2
  pres1@tfPolyDrawOrder   = "postDraw"
 ; pres1@gsLineColor      = "blue"         ; color of lines
 ; dumL0 = gsn_add_polyline(wks,plot(0), (/x_0,x_1/),(/y_0,y_1/) ,pres1)
 ; dumL2 = gsn_add_polyline(wks,plot2, (/x_0,x_1/),(/y_0,y_1/) ,pres1)
 ; dumL4 = gsn_add_polyline(wks,plot(2), (/x_0,x_1/),(/y_0,y_1/) ,pres1)
  
  x_2 = -1.
  x_3 = 12.
  
  y_2 = 0.6334947 
  y_3 = 0.6334947  
  pres1@gsLineColor      = "red"         ; color of lines
 ; dumL1 = gsn_add_polyline(wks,plot(0), (/x_2,x_3/),(/y_2,y_3/) ,pres1)
 ; dumL3 = gsn_add_polyline(wks,plot(1), (/x_2,x_3/),(/y_2,y_3/) ,pres1)
 ; dumL5 = gsn_add_polyline(wks,plot(2), (/x_2,x_3/),(/y_2,y_3/) ,pres1)
;===============================================================
;********************************************************** 
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
  amres5@amParallelPosF   =  -0.68 ; 0.5 is the right edge of the plot, so
  amres5@amOrthogonalPosF =  -0.55    ; This is the center of the plot.
  amres5@amJust           = "TopLeft"

 tt = (/"b","c","d","e","f","g","g"/)
create_tt = gsn_create_text(wks,tt,txres5);
; tt_p0 = new(nplots,graphic)                ; create graphic array ; only needed if paneling
tt_p0 = new(1,graphic)                ; create graphic array ; only needed if paneling
;---Attach the polyline
     do i = 0,0
	print(i)
        tt_p0(i) = gsn_add_annotation(plot(i),create_tt(i),amres5)
     end do
delete(i)
;=============================================================================  
;----------------------------------------------------
  resP                                  = True         ; modify the panel plot
  resP@gsnDraw                          = True      ; don't draw yet
  resP@gsnFrame                         = False        ; don't advance frame yet
  resP@gsnPanelLabelBar    	        = False       ; add common colorbar
  resP@gsnPanelYWhiteSpacePercent      	= 4.0    ; default is 1.0
  resP@gsnPanelXWhiteSpacePercent      	= 4.0    ; default is 1.0
  resP@gsnPanelTop                      = 1.0
  resP@gsnPanelBottom                   = 0.11
  resP@gsnPanelLeft                     = 0.01
  resP@gsnPanelRight                    = 0.99

  resP@lbTopMarginF                     = 0.5
  resP@gsnMaximize                      = True         ; large format
  resP@gsnPaperOrientation              = "portrait"   ; force portrait
  resP@gsnLeftString 			= ""
 
;*******************************************
   
  gsn_panel(wks,plot(0:),(/1,1/),resP)     ; now draw as one plot
  frame(wks)
end
