begin
;======================================
  yrStrt = 1982
  yrLast = 2022
  
  latS   = -60.
  latN   =  60.
  lonL   = -180.
  lonR   =  180.
  
  nmos = 12	;months in a year
  thres = 1.0	;set a threshold
  
  xlabel = ispan(0,nmos-1,1)
  ;print(xlabel)
  
  p95 = 0.05
 
  diri = "/home/ngyuk/Downloads/GEOMAR/Final_work/"
  
  sst_data = "HadISST_sst.nc"
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
  printVarSummary(hadi)       ; variable overview
  ;printMinMax(hadi,True)
;==========================
  hadiClm      = clmMonTLL(hadi)
  printVarSummary(hadiClm)                   ; (12,nlat,nlon)   
  
  hadi         = calcMonAnomTLL(hadi,hadiClm)    ; replace with anonamlies
  hadi@long_name = "ANOMALIES: "+hadi@long_name
  ;delete(hadiClm)
  printMinMax(hadi,True)
  printVarSummary(hadi)
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
  printVarSummary(hadi_s)
  ;printMinMax(hadi_s,True)
;---------------------------------------------------------------------------
  ahadi_s = wgt_areaave_Wrap(hadi_s(time|:,month|:,{lon|-20:0}, {lat|-3:3}),1.0, 1.0, 0)        ;atl3 index
  printVarSummary(ahadi_s)
  ;printMinMax(ahadi_s,True)
;----------------------------
  delete(hadi)
  delete(hadi_s)

  jja_ahadi = dim_avg_n_Wrap(ahadi_s(time|:,month|5:7),1)		;JJA, ncl counts from zero
  print(jja_ahadi)

   jja_ahadi = dim_standardize_n(jja_ahadi, 0, 0)  ;Normalise deviations
   print(jja_ahadi)
;----------------------------------------------------------------------------
; Positive cases
;============================================================================
  n_atlP = dim_num_n(jja_ahadi.ge.thres,0)        ;Number of cases
  print(n_atlP)

  valsP = ind(jja_ahadi.ge.thres)
  print(valsP)
  print(valsP+yrStrt)		;the years
;----------------------------------------------------------------------------
; Negative cases,
; note the "-" sign attached to the threshold here
;============================================================================
  n_atlN = dim_num_n(jja_ahadi.le.-thres,0)        ;Number of cases
  print(n_atlN)

  valsN = ind(jja_ahadi.le.-thres)
  ;print(valsN)
  print(valsN+yrStrt)             ;the years
;==============================================================
;==============================================================
  f      = addfile("/home/ngyuk/Downloads/GEOMAR/Final_work/sst_ensemble_data_results/sst1_mean.nc", "r")
  TIME   = f->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  sst_1   = short2flt(f->sst(iYYYY,:,:))
  
  sst_1 = lonFlip(sst_1)
  
  sst_1 = sst_1 - 273.15        
  
  sst_1!1 = "lat"
  sst_1!2 = "lon"
  printVarSummary(sst_1)       ; variable overview
		
  delete(f)		;	Not needed anymore
  delete(TIME)
  delete(YYYY)
  delete(iYYYY)

  ;printMinMax(sst_1,True)
 ; rain = where(rain.gt.-30.,rain,rain@_FillValue)	;removes padded values
 ; printVarSummary(rain)       ; variable overview
  ;printMinMax(rain,True)
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
  ;printVarSummary(sst_1_s)


  sst_1_s = dtrend_msg_n(sst_1_s&time,sst_1_s,True,False,0)        ;removes linear trends
  ;printVarSummary(sst_1_s)
  ;printMinMax(sst_1_s,True)
;==============================================================
;=================== Canonical
;==============================================================
;***************************************************************************************************************
  dimZ = dimsizes(sst_1_s)
  printVarSummary(dimZ)
      
        raindP_cases_maps  = new((/n_atlP,dimZ(1),dimZ(2),dimZ(3)/),float,sst_1_s@_FillValue)
        printVarSummary(raindP_cases_maps)
    
        raindP_cases_maps = sst_1_s(valsP,:,:,:)
        printVarSummary(raindP_cases_maps)
        
        rainb1p = (raindP_cases_maps(time|2,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1p)    
  
        rainb2p = (raindP_cases_maps(time|3,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2p)    
  
        rainb3p = (raindP_cases_maps(time|4,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3p)    
  
       rainb4p = (raindP_cases_maps(time|6,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4p) 
  
        rainb5p = (raindP_cases_maps(time|7,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb5p) 
         
       rain_fp = (/rainb1p,rainb2p,rainb3p,rainb4p,rainb5p/)        
       copy_VarCoords_n(sst_1_s,rain_fp,(/1,2,3/))
       rain_fp!0 = "time"
       printVarSummary(rain_fp)   
;===================================================================
        raindN_cases_maps  = new((/n_atlN,dimZ(1),dimZ(2),dimZ(3)/),float,sst_1_s@_FillValue)
        printVarSummary(raindN_cases_maps)
    
        raindN_cases_maps = sst_1_s(valsN,:,:,:)
        printVarSummary(raindN_cases_maps)
        
        rainb1n = (raindN_cases_maps(time|0,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1n)    
  
        rainb2n = (raindN_cases_maps(time|1,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2n)    
  
        rainb3n = (raindN_cases_maps(time|2,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3n)    
  
        rainb4n = (raindN_cases_maps(time|5,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4n) 
       
       rain_fn = (/rainb1n,rainb2n,rainb3n,rainb4n/)       
 
       copy_VarCoords_n(sst_1_s,rain_fn,(/1,2,3/))
       rain_fn!0 = "time"   
       printVarSummary(rain_fn) 
;==============================================================
  pc1SST = (dim_avg_n_Wrap(rain_fp,0) - dim_avg_n_Wrap(rain_fn,0))/(2.)                       ; Difference
  copy_VarCoords(sst_1_s(0,:,:,:),pc1SST)                             ; copy coordinate arrays 
  printVarSummary(pc1SST)
  printMinMax(pc1SST,True)
;===================================== 
  pc1SSTd = pc1SST  ; Transfer metdata for plotting sig. values        

  pc1SSTd = ttest((dim_avg_n_Wrap(rain_fp,0)),(dim_stddev_n_Wrap(rain_fp,0))^2,5\
               ,(dim_avg_n_Wrap(rain_fn,0)),(dim_stddev_n_Wrap(rain_fn,0))^2,4,True, False)
  pc1SSTd@long_name = "Probability: difference between means"
  printVarSummary(pc1SSTd)

  pc1SSTd = where(abs(pc1SSTd).le.p95,pc1SSTd,pc1SSTd@_FillValue) ; Sig.u winds, where UWND or VWND is significanlt
  ;print("==============================================================")  
  printVarSummary(pc1SSTd)
  printMinMax(pc1SSTd,True)
  
  d_hadiP6 = dim_avg_n_Wrap(pc1SST(month|:,{lat|-3:3}, lon|:),1)
  printVarSummary(d_hadiP6)
  d_hadidP1 = dim_avg_n_Wrap(pc1SSTd(month|:,{lat|-3:3}, lon|:),1)
  printVarSummary(d_hadidP1)
;==============================================================
;==========Non-canonical
;==============================================================
        rainb1p1 = (raindP_cases_maps(time|0,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1p)    
  
        rainb2p1 = (raindP_cases_maps(time|1,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2p1)    
  
        rainb3p1 = (raindP_cases_maps(time|5,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3p1)    
  
        rainb4p1 = (raindP_cases_maps(time|8,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4p1) 
         
       rain_fp1 = (/rainb1p1,rainb2p1,rainb3p1,rainb4p1/)       
       copy_VarCoords_n(sst_1_s,rain_fp1,(/1,2,3/))
       rain_fp1!0 = "time"    
       printVarSummary(rain_fp1) 
;=========================================================================================================        
        rainb1n1 = (raindN_cases_maps(time|3,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1n1)    
  
        rainb2n1 = (raindN_cases_maps(time|4,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2n1)    
  
        rainb3n1 = (raindN_cases_maps(time|6,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3n1)    
  
        rainb4n1 = (raindN_cases_maps(time|7,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4n1) 
      
       rain_fn1 = (/rainb1n1,rainb2n1,rainb3n1,rainb4n1/)       
       copy_VarCoords_n(sst_1_s,rain_fn1,(/1,2,3/))
       rain_fn1!0 = "time"    
       printVarSummary(rain_fn1) 
;==============================================================
  pc1SST111 = (dim_avg_n_Wrap(rain_fp1,0) - dim_avg_n_Wrap(rain_fn1,0))/(2.)
  copy_VarCoords(sst_1_s(0,:,:,:),pc1SST111)                             ; copy coordinate arrays 
  printVarSummary(pc1SST111)
  printMinMax(pc1SST111,True)
;===================================== 
  pc1SSTd111 = pc1SST111  ; Transfer metdata for plotting sig. values        

  pc1SSTd111 = ttest((dim_avg_n_Wrap(rain_fp1,0)),(dim_stddev_n_Wrap(rain_fp1,0))^2,4\
               ,(dim_avg_n_Wrap(rain_fn1,0)),(dim_stddev_n_Wrap(rain_fn1,0))^2,4,True, False)
  pc1SSTd111@long_name = "Probability: difference between means"
  printVarSummary(pc1SSTd111)

  pc1SSTd111 = where(abs(pc1SSTd111).le.p95,pc1SSTd111,pc1SSTd111@_FillValue) ; Sig.u winds, where UWND or VWND is significanlt
  ;print("==============================================================")  
  printVarSummary(pc1SSTd111)
  printMinMax(pc1SSTd111,True)
  
  d_hadiP6111 = dim_avg_n_Wrap(pc1SST111(month|:,{lat|-3:3}, lon|:),1)
  
  d_hadidP1111 = dim_avg_n_Wrap(pc1SSTd111(month|:,{lat|-3:3}, lon|:),1)
  
  printVarSummary(d_hadiP6)
  printVarSummary(d_hadidP1)
  
  printVarSummary(d_hadiP6111)
  printVarSummary(d_hadidP1111)
;==============================================================
  delete(sst_1)		;	Not needed anymore
;==============================================================
  f      = addfile("/home/ngyuk/Downloads/GEOMAR/Final_work/sst_ensemble_data_results/sst2_mean.nc", "r")
  TIME   = f->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  sst_2   = short2flt(f->sst(iYYYY,:,:))
  
  sst_2 = lonFlip(sst_2)
  
  sst_2 = sst_2 - 273.15        
  
  sst_2!1 = "lat"
  sst_2!2 = "lon"
  printVarSummary(sst_2)       ; variable overview
		
  delete(f)		;	Not needed anymore
  delete(TIME)
  delete(YYYY)
  delete(iYYYY)

  ;printMinMax(sst_2,True)
;==========================
  sst_2Clm      = clmMonTLL(sst_2)
  printVarSummary(sst_2Clm)                   ; (12,nlat,nlon)   
  
  sst_2         = calcMonAnomTLL(sst_2,sst_2Clm)    ; replace with anonamlies
  sst_2@long_name = "ANOMALIES: "+sst_2@long_name
  delete(sst_2Clm)
  ;printMinMax(sst_2,True)
  printVarSummary(sst_2)
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
;==============================================================
;***************************************************************************************************************
  dimZ = dimsizes(sst_2_s)
  printVarSummary(dimZ)
      
        raindP_cases_maps2  = new((/n_atlP,dimZ(1),dimZ(2),dimZ(3)/),float,sst_2_s@_FillValue)
        printVarSummary(raindP_cases_maps2)
    
        raindP_cases_maps2 = sst_2_s(valsP,:,:,:)
        printVarSummary(raindP_cases_maps2)
        
        rainb1p2 = (raindP_cases_maps2(time|2,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1p2)    
  
        rainb2p2 = (raindP_cases_maps2(time|3,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2p2)    
  
        rainb3p2 = (raindP_cases_maps2(time|4,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3p2)    
  
       rainb4p2 = (raindP_cases_maps(time|6,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4p2) 
  
        rainb5p2 = (raindP_cases_maps2(time|7,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb5p) 
         
       rain_fp2 = (/rainb1p2,rainb2p2,rainb3p2,rainb4p2,rainb5p2/)       
 
       copy_VarCoords_n(sst_2_s,rain_fp2,(/1,2,3/))
       rain_fp2!0 = "time"  
       printVarSummary(rain_fp2)  
;===================================================================
        raindN_cases_maps2  = new((/n_atlN,dimZ(1),dimZ(2),dimZ(3)/),float,sst_2_s@_FillValue)
        printVarSummary(raindN_cases_maps2)
    
        raindN_cases_maps2 = sst_2_s(valsN,:,:,:)
        printVarSummary(raindN_cases_maps2)
        
        rainb1n2 = (raindN_cases_maps2(time|0,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1n2)    
  
        rainb2n2 = (raindN_cases_maps2(time|1,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2n2)    
  
        rainb3n2 = (raindN_cases_maps2(time|2,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3n2)    
  
        rainb4n2 = (raindN_cases_maps2(time|5,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4n2) 
       
       rain_fn2 = (/rainb1n2,rainb2n2,rainb3n2,rainb4n2/)       
 
       copy_VarCoords_n(sst_2_s,rain_fn2,(/1,2,3/))
       rain_fn2!0 = "time"   
       printVarSummary(rain_fn2) 
;==============================================================
  pc1SST2 = (dim_avg_n_Wrap(rain_fp2,0) - dim_avg_n_Wrap(rain_fn2,0))/(2.)                       ; Difference
  copy_VarCoords(sst_2_s(0,:,:,:),pc1SST2)                             ; copy coordinate arrays 
  printVarSummary(pc1SST2)
  printMinMax(pc1SST2,True)
;===================================== 
  pc1SSTd2 = pc1SST2  ; Transfer metdata for plotting sig. values        

  pc1SSTd2 = ttest((dim_avg_n_Wrap(rain_fp2,0)),(dim_stddev_n_Wrap(rain_fp2,0))^2,5\
               ,(dim_avg_n_Wrap(rain_fn2,0)),(dim_stddev_n_Wrap(rain_fn2,0))^2,4,True, False)
  pc1SSTd2@long_name = "Probability: difference between means"
  printVarSummary(pc1SSTd2)

  pc1SSTd2 = where(abs(pc1SSTd2).le.p95,pc1SSTd2,pc1SSTd2@_FillValue) ; Sig.u winds, where UWND or VWND is significanlt
  ;print("==============================================================")  
  printVarSummary(pc1SSTd2)
  printMinMax(pc1SSTd2,True)
  
  d_hadiP62 = dim_avg_n_Wrap(pc1SST2(month|:,{lat|-3:3}, lon|:),1)
  printVarSummary(d_hadiP62)
  d_hadidP12 = dim_avg_n_Wrap(pc1SSTd2(month|:,{lat|-3:3}, lon|:),1)
  printVarSummary(d_hadidP12)
;==============================================================
;==========Non-canonical
;==============================================================
        rainb1p12 = (raindP_cases_maps2(time|0,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1p2)    
  
        rainb2p12 = (raindP_cases_maps2(time|1,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2p12)    
  
        rainb3p12 = (raindP_cases_maps2(time|5,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3p12)    
  
        rainb4p12 = (raindP_cases_maps2(time|8,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4p12) 
         
       rain_fp12 = (/rainb1p12,rainb2p12,rainb3p12,rainb4p12/)       
       copy_VarCoords_n(sst_2_s,rain_fp12,(/1,2,3/))
       rain_fp12!0 = "time"    
       printVarSummary(rain_fp12) 
;=========================================================================================================        
        rainb1n12 = (raindN_cases_maps2(time|3,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1n12)    
  
        rainb2n12 = (raindN_cases_maps2(time|4,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2n12)    
  
        rainb3n12 = (raindN_cases_maps2(time|6,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3n12)    
  
        rainb4n12 = (raindN_cases_maps2(time|7,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4n12) 
      
       rain_fn12 = (/rainb1n12,rainb2n12,rainb3n12,rainb4n12/)       
       copy_VarCoords_n(sst_2_s,rain_fn12,(/1,2,3/))
       rain_fn12!0 = "time"    
       printVarSummary(rain_fn12) 
;==============================================================
  pc1SST1112 = (dim_avg_n_Wrap(rain_fp12,0) - dim_avg_n_Wrap(rain_fn12,0))/(2.)
  copy_VarCoords(sst_2_s(0,:,:,:),pc1SST1112)                             ; copy coordinate arrays 
  printVarSummary(pc1SST1112)
  printMinMax(pc1SST1112,True)
;===================================== 
  pc1SSTd1112 = pc1SST1112  ; Transfer metdata for plotting sig. values        

  pc1SSTd1112 = ttest((dim_avg_n_Wrap(rain_fp12,0)),(dim_stddev_n_Wrap(rain_fp12,0))^2,4\
               ,(dim_avg_n_Wrap(rain_fn12,0)),(dim_stddev_n_Wrap(rain_fn12,0))^2,4,True, False)
  pc1SSTd1112@long_name = "Probability: difference between means"
  printVarSummary(pc1SSTd1112)

  pc1SSTd1112 = where(abs(pc1SSTd1112).le.p95,pc1SSTd1112,pc1SSTd1112@_FillValue) ; Sig.u winds, where UWND or VWND is significanlt
  ;print("==============================================================")  
  printVarSummary(pc1SSTd1112)
  printMinMax(pc1SSTd1112,True)
  
  d_hadiP61112 = dim_avg_n_Wrap(pc1SST1112(month|:,{lat|-3:3}, lon|:),1)
  
  d_hadidP11112 = dim_avg_n_Wrap(pc1SSTd1112(month|:,{lat|-3:3}, lon|:),1)
  
  printVarSummary(d_hadiP62)
  printVarSummary(d_hadidP12)
  
  printVarSummary(d_hadiP61112)
  printVarSummary(d_hadidP11112)
;==============================================================
 delete(sst_2)		;	Not needed anymore
;==============================================================
  f      = addfile("/home/ngyuk/Downloads/GEOMAR/Final_work/sst_ensemble_data_results/sst3_mean.nc", "r")
  TIME   = f->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  sst_3   = short2flt(f->sst(iYYYY,:,:))
  
  sst_3 = lonFlip(sst_3)
  
  sst_3 = sst_3 - 273.15        
  
  sst_3!1 = "lat"
  sst_3!2 = "lon"
  printVarSummary(sst_3)       ; variable overview
		
  delete(f)		;	Not needed anymore
  delete(TIME)
  delete(YYYY)
  delete(iYYYY)

  ;printMinMax(sst_3,True)
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
;==============================================================
;***************************************************************************************************************
  dimZ = dimsizes(sst_3_s)
  printVarSummary(dimZ)
      
        raindP_cases_maps3  = new((/n_atlP,dimZ(1),dimZ(2),dimZ(3)/),float,sst_3_s@_FillValue)
        printVarSummary(raindP_cases_maps3)
    
        raindP_cases_maps3 = sst_3_s(valsP,:,:,:)
        printVarSummary(raindP_cases_maps3)
        
        rainb1p3 = (raindP_cases_maps3(time|2,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1p3)    
  
        rainb2p3 = (raindP_cases_maps3(time|3,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2p3)    
  
        rainb3p3 = (raindP_cases_maps3(time|4,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3p3)    
  
       rainb4p3 = (raindP_cases_maps3(time|6,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4p3) 
  
        rainb5p3 = (raindP_cases_maps3(time|7,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb5p3) 
         
       rain_fp3 = (/rainb1p3,rainb2p3,rainb3p3,rainb4p3,rainb5p3/)        
       copy_VarCoords_n(sst_3_s,rain_fp3,(/1,2,3/))
       rain_fp3!0 = "time"    
       printVarSummary(rain_fp3) 
;===================================================================
        raindN_cases_maps3  = new((/n_atlN,dimZ(1),dimZ(2),dimZ(3)/),float,sst_3_s@_FillValue)
        printVarSummary(raindN_cases_maps3)
    
        raindN_cases_maps3 = sst_3_s(valsN,:,:,:)
        printVarSummary(raindN_cases_maps3)
        
        rainb1n3 = (raindN_cases_maps3(time|0,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1n3)    
  
        rainb2n3 = (raindN_cases_maps3(time|1,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2n3)    
  
        rainb3n3 = (raindN_cases_maps3(time|2,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3n3)    
  
        rainb4n3 = (raindN_cases_maps3(time|5,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4n3) 
       
       rain_fn3 = (/rainb1n3,rainb2n3,rainb3n3,rainb4n3/)        
       copy_VarCoords_n(sst_3_s,rain_fn3,(/1,2,3/))
       rain_fn3!0 = "time"    
       printVarSummary(rain_fn3) 
;==============================================================
  pc1SST3 = (dim_avg_n_Wrap(rain_fp3,0) - dim_avg_n_Wrap(rain_fn3,0))/(2.)                       ; Difference
  copy_VarCoords(sst_3_s(0,:,:,:),pc1SST3)                             ; copy coordinate arrays 
  printVarSummary(pc1SST3)
  printMinMax(pc1SST3,True)
;===================================== 
  pc1SSTd3 = pc1SST3  ; Transfer metdata for plotting sig. values        

  pc1SSTd3 = ttest((dim_avg_n_Wrap(rain_fp3,0)),(dim_stddev_n_Wrap(rain_fp3,0))^2,5\
               ,(dim_avg_n_Wrap(rain_fn3,0)),(dim_stddev_n_Wrap(rain_fn3,0))^2,4,True, False)
  pc1SSTd3@long_name = "Probability: difference between means"
  printVarSummary(pc1SSTd3)

  pc1SSTd3 = where(abs(pc1SSTd3).le.p95,pc1SSTd3,pc1SSTd3@_FillValue) ; Sig.u winds, where UWND or VWND is significanlt
  ;print("==============================================================")  
  printVarSummary(pc1SSTd3)
  printMinMax(pc1SSTd3,True)
  
  d_hadiP63 = dim_avg_n_Wrap(pc1SST3(month|:,{lat|-3:3}, lon|:),1)
  printVarSummary(d_hadiP63)
  d_hadidP13 = dim_avg_n_Wrap(pc1SSTd3(month|:,{lat|-3:3}, lon|:),1)
  printVarSummary(d_hadidP13)
;==============================================================
;==========Non-canonical
;==============================================================
        rainb1p13 = (raindP_cases_maps3(time|0,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1p3)    
  
        rainb2p13 = (raindP_cases_maps3(time|1,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2p13)    
  
        rainb3p13 = (raindP_cases_maps3(time|5,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3p13)    
  
        rainb4p13 = (raindP_cases_maps3(time|8,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4p13) 
         
       rain_fp13 = (/rainb1p13,rainb2p13,rainb3p13,rainb4p13/)       
       copy_VarCoords_n(sst_3_s,rain_fp13,(/1,2,3/))
       rain_fp13!0 = "time"    
       printVarSummary(rain_fp13) 
;=========================================================================================================        
        rainb1n13 = (raindN_cases_maps3(time|3,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1n13)    
  
        rainb2n13 = (raindN_cases_maps3(time|4,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2n13)    
  
        rainb3n13 = (raindN_cases_maps3(time|6,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3n13)    
  
        rainb4n13 = (raindN_cases_maps3(time|7,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4n13) 
      
       rain_fn13 = (/rainb1n13,rainb2n13,rainb3n13,rainb4n13/)       
       copy_VarCoords_n(sst_3_s,rain_fn13,(/1,2,3/))
       rain_fn13!0 = "time"    
       printVarSummary(rain_fn13) 
;==============================================================
  pc1SST1113 = (dim_avg_n_Wrap(rain_fp13,0) - dim_avg_n_Wrap(rain_fn13,0))/(2.)
  copy_VarCoords(sst_3_s(0,:,:,:),pc1SST1113)                             ; copy coordinate arrays 
  printVarSummary(pc1SST1113)
  printMinMax(pc1SST1113,True)
;===================================== 
  pc1SSTd1113 = pc1SST1113  ; Transfer metdata for plotting sig. values        

  pc1SSTd1113 = ttest((dim_avg_n_Wrap(rain_fp13,0)),(dim_stddev_n_Wrap(rain_fp13,0))^2,4\
               ,(dim_avg_n_Wrap(rain_fn13,0)),(dim_stddev_n_Wrap(rain_fn13,0))^2,4,True, False)
  pc1SSTd1113@long_name = "Probability: difference between means"
  printVarSummary(pc1SSTd1113)

  pc1SSTd1113 = where(abs(pc1SSTd1113).le.p95,pc1SSTd1113,pc1SSTd1113@_FillValue) ; Sig.u winds, where UWND or VWND is significanlt
  ;print("==============================================================")  
  printVarSummary(pc1SSTd1113)
  printMinMax(pc1SSTd1113,True)
  
  d_hadiP61113 = dim_avg_n_Wrap(pc1SST1113(month|:,{lat|-3:3}, lon|:),1)
  
  d_hadidP11113 = dim_avg_n_Wrap(pc1SSTd1113(month|:,{lat|-3:3}, lon|:),1)
  
  printVarSummary(d_hadiP63)
  printVarSummary(d_hadidP13)
  
  printVarSummary(d_hadiP61113)
  printVarSummary(d_hadidP11113)
;==============================================================
   delete(sst_3)		;	Not needed anymore
;==============================================================
  f      = addfile("/home/ngyuk/Downloads/GEOMAR/Final_work/sst_ensemble_data_results/sst4_mean.nc", "r")
  TIME   = f->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  sst_4   = short2flt(f->sst(iYYYY,:,:))
  
  sst_4 = lonFlip(sst_4)
  
  sst_4 = sst_4 - 273.15        
  
  sst_4!1 = "lat"
  sst_4!2 = "lon"
  printVarSummary(sst_4)       ; variable overview
		
  delete(f)		;	Not needed anymore
  delete(TIME)
  delete(YYYY)
  delete(iYYYY)

  ;printMinMax(sst_4,True)
;==========================
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
  ;printMinMax(sst_4_s,True)
;==============================================================
;=================== Canonical
;==============================================================
;***************************************************************************************************************
  dimZ = dimsizes(sst_4_s)
  printVarSummary(dimZ)
      
        raindP_cases_maps4  = new((/n_atlP,dimZ(1),dimZ(2),dimZ(3)/),float,sst_4_s@_FillValue)
        printVarSummary(raindP_cases_maps4)
    
        raindP_cases_maps4 = sst_4_s(valsP,:,:,:)
        printVarSummary(raindP_cases_maps4)
        
        rainb1p4 = (raindP_cases_maps4(time|2,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1p4)    
  
        rainb2p4 = (raindP_cases_maps4(time|3,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2p4)    
  
        rainb3p4 = (raindP_cases_maps4(time|4,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3p4)    
  
        rainb4p4 = (raindP_cases_maps4(time|6,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4p4) 
  
        rainb5p4 = (raindP_cases_maps4(time|7,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb5p4) 
         
       rain_fp4 = (/rainb1p4,rainb2p4,rainb3p4,rainb4p4,rainb5p4/)       
 
       copy_VarCoords_n(sst_4_s,rain_fp4,(/1,2,3/))
       rain_fp4!0 = "time"    
       printVarSummary(rain_fp4) 
;===================================================================
        raindN_cases_maps4  = new((/n_atlN,dimZ(1),dimZ(2),dimZ(3)/),float,sst_4_s@_FillValue)
        printVarSummary(raindN_cases_maps4)
    
        raindN_cases_maps4 = sst_4_s(valsN,:,:,:)
        printVarSummary(raindN_cases_maps4)
        
        rainb1n4 = (raindN_cases_maps4(time|0,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1n4)    
  
        rainb2n4 = (raindN_cases_maps4(time|1,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2n4)    
  
        rainb3n4 = (raindN_cases_maps4(time|2,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3n4)    
  
        rainb4n4 = (raindN_cases_maps4(time|5,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4n4) 
       
       rain_fn4 = (/rainb1n4,rainb2n4,rainb3n4,rainb4n4/)       
 
       copy_VarCoords_n(sst_4_s,rain_fn4,(/1,2,3/))
       rain_fn4!0 = "time"    
       printVarSummary(rain_fn4) 
;==============================================================
  pc1SST4 = (dim_avg_n_Wrap(rain_fp4,0) - dim_avg_n_Wrap(rain_fn4,0))/(2.)                       ; Difference
  copy_VarCoords(sst_4_s(0,:,:,:),pc1SST4)                             ; copy coordinate arrays 
  printVarSummary(pc1SST4)
  printMinMax(pc1SST4,True)
;===================================== 
  pc1SSTd4 = pc1SST4  ; Transfer metdata for plotting sig. values        

  pc1SSTd4 = ttest((dim_avg_n_Wrap(rain_fp4,0)),(dim_stddev_n_Wrap(rain_fp4,0))^2,5\
               ,(dim_avg_n_Wrap(rain_fn4,0)),(dim_stddev_n_Wrap(rain_fn4,0))^2,4,True, False)
  pc1SSTd4@long_name = "Probability: difference between means"
  printVarSummary(pc1SSTd4)

  pc1SSTd4 = where(abs(pc1SSTd4).le.p95,pc1SSTd4,pc1SSTd4@_FillValue) ; Sig.u winds, where UWND or VWND is significanlt
  ;print("==============================================================")  
  printVarSummary(pc1SSTd4)
  printMinMax(pc1SSTd4,True)
  
  d_hadiP64 = dim_avg_n_Wrap(pc1SST4(month|:,{lat|-3:3}, lon|:),1)
  printVarSummary(d_hadiP64)
  d_hadidP14 = dim_avg_n_Wrap(pc1SSTd4(month|:,{lat|-3:3}, lon|:),1)
  printVarSummary(d_hadidP14)
;==============================================================
;==========Non-canonical
;==============================================================
        rainb1p14 = (raindP_cases_maps4(time|0,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1p4)    
  
        rainb2p14 = (raindP_cases_maps4(time|1,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2p14)    
  
        rainb3p14 = (raindP_cases_maps4(time|5,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3p14)    
  
        rainb4p14 = (raindP_cases_maps4(time|8,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4p14) 
         
       rain_fp14 = (/rainb1p14,rainb2p14,rainb3p14,rainb4p14/)       
       copy_VarCoords_n(sst_4_s,rain_fp14,(/1,2,3/))
       rain_fp14!0 = "time"    
       printVarSummary(rain_fp14) 
;=========================================================================================================        
        rainb1n14 = (raindN_cases_maps4(time|3,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1n14)    
  
        rainb2n14 = (raindN_cases_maps4(time|4,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2n14)    
  
        rainb3n14 = (raindN_cases_maps4(time|6,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3n14)    
  
        rainb4n14 = (raindN_cases_maps4(time|7,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4n14) 
      
       rain_fn14 = (/rainb1n14,rainb2n14,rainb3n14,rainb4n14/)       
       copy_VarCoords_n(sst_4_s,rain_fn14,(/1,2,3/))
       rain_fn14!0 = "time"    
       printVarSummary(rain_fn14) 
;==============================================================
  pc1SST1114 = (dim_avg_n_Wrap(rain_fp14,0) - dim_avg_n_Wrap(rain_fn14,0))/(2.)
  copy_VarCoords(sst_4_s(0,:,:,:),pc1SST1114)                             ; copy coordinate arrays 
  printVarSummary(pc1SST1114)
  printMinMax(pc1SST1114,True)
;===================================== 
  pc1SSTd1114 = pc1SST1114  ; Transfer metdata for plotting sig. values        

  pc1SSTd1114 = ttest((dim_avg_n_Wrap(rain_fp14,0)),(dim_stddev_n_Wrap(rain_fp14,0))^2,4\
               ,(dim_avg_n_Wrap(rain_fn14,0)),(dim_stddev_n_Wrap(rain_fn14,0))^2,4,True, False)
  pc1SSTd1114@long_name = "Probability: difference between means"
  printVarSummary(pc1SSTd1114)

  pc1SSTd1114 = where(abs(pc1SSTd1114).le.p95,pc1SSTd1114,pc1SSTd1114@_FillValue) ; Sig.u winds, where UWND or VWND is significanlt
  ;print("==============================================================")  
  printVarSummary(pc1SSTd1114)
  printMinMax(pc1SSTd1114,True)
  
  d_hadiP61114 = dim_avg_n_Wrap(pc1SST1114(month|:,{lat|-3:3}, lon|:),1)
  
  d_hadidP11114 = dim_avg_n_Wrap(pc1SSTd1114(month|:,{lat|-3:3}, lon|:),1)
  
  printVarSummary(d_hadiP64)
  printVarSummary(d_hadidP14)
  
  printVarSummary(d_hadiP61114)
  printVarSummary(d_hadidP11114)
;==============================================================
  delete(sst_4)		;	Not needed anymore
;==============================================================
  f      = addfile("/home/ngyuk/Downloads/GEOMAR/Final_work/sst_ensemble_data_results/sst5_mean.nc", "r")
  TIME   = f->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  sst_5   = short2flt(f->sst(iYYYY,:,:))
  
  sst_5 = lonFlip(sst_5)
  
  sst_5 = sst_5 - 273.15        
  
  sst_5!1 = "lat"
  sst_5!2 = "lon"
  printVarSummary(sst_5)       ; variable overview
		
  delete(f)		;	Not needed anymore
  delete(TIME)
  delete(YYYY)
  delete(iYYYY)

  ;printMinMax(sst_5,True)
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
  ;printMinMax(sst_5_s,True)
;==============================================================
;=================== Canonical
;==============================================================
;***************************************************************************************************************
  dimZ = dimsizes(sst_5_s)
  printVarSummary(dimZ)
      
        raindP_cases_maps5  = new((/n_atlP,dimZ(1),dimZ(2),dimZ(3)/),float,sst_5_s@_FillValue)
        printVarSummary(raindP_cases_maps5)
    
        raindP_cases_maps5 = sst_5_s(valsP,:,:,:)
        printVarSummary(raindP_cases_maps5)
        
        rainb1p5 = (raindP_cases_maps5(time|2,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1p5)    
  
        rainb2p5 = (raindP_cases_maps5(time|3,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2p5)    
  
        rainb3p5 = (raindP_cases_maps5(time|4,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3p5)    
  
        rainb4p5 = (raindP_cases_maps5(time|6,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4p5) 
  
        rainb5p5 = (raindP_cases_maps5(time|7,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb5p5) 
         
       rain_fp5 = (/rainb1p5,rainb2p5,rainb3p5,rainb4p5,rainb5p5/)       
       copy_VarCoords_n(sst_5_s,rain_fp5,(/1,2,3/))
       rain_fp5!0 = "time"    
       printVarSummary(rain_fp5) 
;===================================================================
        raindN_cases_maps5  = new((/n_atlN,dimZ(1),dimZ(2),dimZ(3)/),float,sst_5_s@_FillValue)
        printVarSummary(raindN_cases_maps5)
    
        raindN_cases_maps5 = sst_5_s(valsN,:,:,:)
        printVarSummary(raindN_cases_maps5)
        
        rainb1n5 = (raindN_cases_maps5(time|0,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1n5)    
  
        rainb2n5 = (raindN_cases_maps5(time|1,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2n5)    
  
        rainb3n5 = (raindN_cases_maps5(time|2,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3n5)    
  
        rainb4n5 = (raindN_cases_maps5(time|5,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4n5) 
       
       rain_fn5 = (/rainb1n5,rainb2n5,rainb3n5,rainb4n5/)       
       copy_VarCoords_n(sst_5_s,rain_fn5,(/1,2,3/))
       rain_fn5!0 = "time"   
       printVarSummary(rain_fn5) 
;==============================================================
  pc1SST5 = (dim_avg_n_Wrap(rain_fp5,0) - dim_avg_n_Wrap(rain_fn5,0))/(2.)                       ; Difference
  copy_VarCoords(sst_5_s(0,:,:,:),pc1SST5)                             ; copy coordinate arrays 
  printVarSummary(pc1SST5)
  printMinMax(pc1SST5,True)
;===================================== 
  pc1SSTd5 = pc1SST5  ; Transfer metdata for plotting sig. values        

  pc1SSTd5 = ttest((dim_avg_n_Wrap(rain_fp5,0)),(dim_stddev_n_Wrap(rain_fp5,0))^2,5\
               ,(dim_avg_n_Wrap(rain_fn5,0)),(dim_stddev_n_Wrap(rain_fn5,0))^2,4,True, False)
  pc1SSTd5@long_name = "Probability: difference between means"
  printVarSummary(pc1SSTd5)

  pc1SSTd5 = where(abs(pc1SSTd5).le.p95,pc1SSTd5,pc1SSTd5@_FillValue) ; Sig.u winds, where UWND or VWND is significanlt
  ;print("==============================================================")  
  printVarSummary(pc1SSTd5)
  printMinMax(pc1SSTd5,True)
  
  d_hadiP65 = dim_avg_n_Wrap(pc1SST5(month|:,{lat|-3:3}, lon|:),1)
  printVarSummary(d_hadiP65)
  d_hadidP15 = dim_avg_n_Wrap(pc1SSTd5(month|:,{lat|-3:3}, lon|:),1)
  printVarSummary(d_hadidP15)
;==============================================================
;==========Non-canonical
;==============================================================
        rainb1p15 = (raindP_cases_maps5(time|0,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1p5)    
  
        rainb2p15 = (raindP_cases_maps5(time|1,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2p15)    
  
        rainb3p15 = (raindP_cases_maps5(time|5,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3p15)    
  
        rainb4p15 = (raindP_cases_maps5(time|8,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4p15) 
         
       rain_fp15 = (/rainb1p15,rainb2p15,rainb3p15,rainb4p15/)       
       copy_VarCoords_n(sst_5_s,rain_fp15,(/1,2,3/))
       rain_fp15!0 = "time"    
       printVarSummary(rain_fp15) 
;=========================================================================================================        
        rainb1n15 = (raindN_cases_maps5(time|3,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1n15)    
  
        rainb2n15 = (raindN_cases_maps5(time|4,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2n15)    
  
        rainb3n15 = (raindN_cases_maps5(time|6,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3n15)    
  
        rainb4n15 = (raindN_cases_maps5(time|7,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4n15) 
      
       rain_fn15 = (/rainb1n15,rainb2n15,rainb3n15,rainb4n15/)       
       copy_VarCoords_n(sst_5_s,rain_fn15,(/1,2,3/))
       rain_fn15!0 = "time"    
       printVarSummary(rain_fn15) 
;==============================================================
  pc1SST1115 = (dim_avg_n_Wrap(rain_fp15,0) - dim_avg_n_Wrap(rain_fn15,0))/(2.)
  copy_VarCoords(sst_5_s(0,:,:,:),pc1SST1115)                             ; copy coordinate arrays 
  printVarSummary(pc1SST1115)
  printMinMax(pc1SST1115,True)
;===================================== 
  pc1SSTd1115 = pc1SST1115  ; Transfer metdata for plotting sig. values        

  pc1SSTd1115 = ttest((dim_avg_n_Wrap(rain_fp15,0)),(dim_stddev_n_Wrap(rain_fp15,0))^2,4\
               ,(dim_avg_n_Wrap(rain_fn15,0)),(dim_stddev_n_Wrap(rain_fn15,0))^2,4,True, False)
  pc1SSTd1115@long_name = "Probability: difference between means"
  printVarSummary(pc1SSTd1115)

  pc1SSTd1115 = where(abs(pc1SSTd1115).le.p95,pc1SSTd1115,pc1SSTd1115@_FillValue) ; Sig.u winds, where UWND or VWND is significanlt
  ;print("==============================================================")  
  printVarSummary(pc1SSTd1115)
  printMinMax(pc1SSTd1115,True)
  
  d_hadiP61115 = dim_avg_n_Wrap(pc1SST1115(month|:,{lat|-3:3}, lon|:),1)
  
  d_hadidP11115 = dim_avg_n_Wrap(pc1SSTd1115(month|:,{lat|-3:3}, lon|:),1)
  
  printVarSummary(d_hadiP65)
  printVarSummary(d_hadidP15)
  
  printVarSummary(d_hadiP61115)
  printVarSummary(d_hadidP11115)
;==============================================================
  delete(sst_5)		;	Not needed anymore
;==============================================================
  f      = addfile("/home/ngyuk/Downloads/GEOMAR/Final_work/sst_ensemble_data_results/sst6_mean.nc", "r")
  TIME   = f->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  sst_6   = short2flt(f->sst(iYYYY,:,:))
  
  sst_6 = lonFlip(sst_6)
  
  sst_6 = sst_6 - 273.15        
  
  sst_6!1 = "lat"
  sst_6!2 = "lon"
  printVarSummary(sst_6)       ; variable overview
		
  delete(f)		;	Not needed anymore
  delete(TIME)
  delete(YYYY)
  delete(iYYYY)

  ;printMinMax(sst_6,True)
;==========================
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
;==============================================================
;***************************************************************************************************************
  dimZ = dimsizes(sst_6_s)
  printVarSummary(dimZ)
      
        raindP_cases_maps6  = new((/n_atlP,dimZ(1),dimZ(2),dimZ(3)/),float,sst_6_s@_FillValue)
        printVarSummary(raindP_cases_maps)
    
        raindP_cases_maps6 = sst_6_s(valsP,:,:,:)
        printVarSummary(raindP_cases_maps6)
        
        rainb1p6 = (raindP_cases_maps6(time|2,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1p6)    
  
        rainb2p6 = (raindP_cases_maps6(time|3,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2p6)    
  
        rainb3p6 = (raindP_cases_maps6(time|4,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3p6)    
  
       rainb4p6 = (raindP_cases_maps6(time|6,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4p6) 
  
        rainb5p6 = (raindP_cases_maps6(time|7,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb5p6) 
         
       rain_fp6 = (/rainb1p6,rainb2p6,rainb3p6,rainb4p6,rainb5p6/)       
       copy_VarCoords_n(sst_6_s,rain_fp6,(/1,2,3/))
       rain_fp6!0 = "time"   
       printVarSummary(rain_fp6) 
;===================================================================
        raindN_cases_maps6  = new((/n_atlN,dimZ(1),dimZ(2),dimZ(3)/),float,sst_6_s@_FillValue)
        printVarSummary(raindN_cases_maps6)
    
        raindN_cases_maps6 = sst_6_s(valsN,:,:,:)
        printVarSummary(raindN_cases_maps6)
        
        rainb1n6 = (raindN_cases_maps6(time|0,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1n6)    
  
        rainb2n6 = (raindN_cases_maps6(time|1,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2n6)    
  
        rainb3n6 = (raindN_cases_maps6(time|2,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3n6)    
  
        rainb4n6 = (raindN_cases_maps6(time|5,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4n6) 
       
       rain_fn6 = (/rainb1n6,rainb2n6,rainb3n6,rainb4n6/)       
       copy_VarCoords_n(sst_6_s,rain_fn6,(/1,2,3/))
       rain_fn6!0 = "time"    
       printVarSummary(rain_fn6) 
;==============================================================
  pc1SST6 = (dim_avg_n_Wrap(rain_fp6,0) - dim_avg_n_Wrap(rain_fn6,0))/(2.)                       ; Difference
  copy_VarCoords(sst_6_s(0,:,:,:),pc1SST6)                             ; copy coordinate arrays 
  printVarSummary(pc1SST6)
  printMinMax(pc1SST6,True)
;===================================== 
  pc1SSTd6 = pc1SST6  ; Transfer metdata for plotting sig. values        

  pc1SSTd6 = ttest((dim_avg_n_Wrap(rain_fp6,0)),(dim_stddev_n_Wrap(rain_fp6,0))^2,5\
               ,(dim_avg_n_Wrap(rain_fn6,0)),(dim_stddev_n_Wrap(rain_fn6,0))^2,4,True, False)
  pc1SSTd6@long_name = "Probability: difference between means"
  printVarSummary(pc1SSTd6)

  pc1SSTd6 = where(abs(pc1SSTd6).le.p95,pc1SSTd6,pc1SSTd6@_FillValue) ; Sig.u winds, where UWND or VWND is significanlt
  ;print("==============================================================")  
  printVarSummary(pc1SSTd6)
  printMinMax(pc1SSTd6,True)
  
  d_hadiP66 = dim_avg_n_Wrap(pc1SST6(month|:,{lat|-3:3}, lon|:),1)
  printVarSummary(d_hadiP66)
  d_hadidP16 = dim_avg_n_Wrap(pc1SSTd6(month|:,{lat|-3:3}, lon|:),1)
  printVarSummary(d_hadidP16)
;==============================================================
;==========Non-canonical
;==============================================================
        rainb1p16 = (raindP_cases_maps6(time|0,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1p6)    
  
        rainb2p16 = (raindP_cases_maps6(time|1,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2p16)    
  
        rainb3p16 = (raindP_cases_maps6(time|5,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3p16)    
  
        rainb4p16 = (raindP_cases_maps6(time|8,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4p16) 
         
       rain_fp16 = (/rainb1p16,rainb2p16,rainb3p16,rainb4p16/)       
       copy_VarCoords_n(sst_6_s,rain_fp16,(/1,2,3/))
       rain_fp16!0 = "time"    
       printVarSummary(rain_fp16) 
;=========================================================================================================        
        rainb1n16 = (raindN_cases_maps6(time|3,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1n16)    
  
        rainb2n16 = (raindN_cases_maps6(time|4,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2n16)    
  
        rainb3n16 = (raindN_cases_maps6(time|6,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3n16)    
  
        rainb4n16 = (raindN_cases_maps6(time|7,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4n16) 
      
       rain_fn16 = (/rainb1n16,rainb2n16,rainb3n16,rainb4n16/)       
       copy_VarCoords_n(sst_6_s,rain_fn16,(/1,2,3/))
       rain_fn16!0 = "time"    
       printVarSummary(rain_fn16) 
;==============================================================
  pc1SST1116 = (dim_avg_n_Wrap(rain_fp16,0) - dim_avg_n_Wrap(rain_fn16,0))/(2.)
  copy_VarCoords(sst_6_s(0,:,:,:),pc1SST1116)                             ; copy coordinate arrays 
  printVarSummary(pc1SST1116)
  printMinMax(pc1SST1116,True)
;===================================== 
  pc1SSTd1116 = pc1SST1116  ; Transfer metdata for plotting sig. values        

  pc1SSTd1116 = ttest((dim_avg_n_Wrap(rain_fp16,0)),(dim_stddev_n_Wrap(rain_fp16,0))^2,4\
               ,(dim_avg_n_Wrap(rain_fn16,0)),(dim_stddev_n_Wrap(rain_fn16,0))^2,4,True, False)
  pc1SSTd1116@long_name = "Probability: difference between means"
  printVarSummary(pc1SSTd1116)

  pc1SSTd1116 = where(abs(pc1SSTd1116).le.p95,pc1SSTd1116,pc1SSTd1116@_FillValue) ; Sig.u winds, where UWND or VWND is significanlt
  ;print("==============================================================")  
  printVarSummary(pc1SSTd1116)
  printMinMax(pc1SSTd1116,True)
  
  d_hadiP61116 = dim_avg_n_Wrap(pc1SST1116(month|:,{lat|-3:3}, lon|:),1)
  
  d_hadidP11116 = dim_avg_n_Wrap(pc1SSTd1116(month|:,{lat|-3:3}, lon|:),1)
  
  printVarSummary(d_hadiP66)
  printVarSummary(d_hadidP16)
  
  printVarSummary(d_hadiP61116)
  printVarSummary(d_hadidP11116)         
;***************************************************************************************************************
;***************************************************************************************************************
  wks = gsn_open_wks("eps","Figure_S5")         ; send graphics to PNG file
  plot = new(12,graphic)                ; create graphic array ; only needed if paneling
  
  ress                      = True         
  ress@gsnDraw              = False        ; don't draw yet
  ress@gsnFrame             = False        ; don't advance frame yet
  ress@gsnMaximize          = False  
  ress@gsnAddCyclic         = False
  
  ress@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels
  ress@cnMinLevelValF       = -1.0            ; set min contour level
  ress@cnMaxLevelValF       =  1.0             ; set max contour level
  ress@cnLevelSpacingF      =  0.1            ; set contour spacing   
  ress@cnFillOn             = True         ; turn on color fill
  ress@cnFillPalette        = "cmp_b2r"  
  ress@cnLevels             =  4.
  
  ress@cnFillOn             = True
  
  ress@cnSmoothingOn        = True
  ress@cnSmoothingDistanceF = 0.005
  ress@cnSmoothingTensionF  = 0.001
  
  ress@tfDoNDCOverlay        = True            ; necessary for correct overlay on map

  ress@cnLinesOn             = False            ; turn off contour lines
  ress@cnLineLabelsOn        = False            ; turn off contour labels
  
  ress@cnFillMode            = "AreaFill"
  ress@cnFillOpacityF        = 1.

  ress@cnInfoLabelOn        = False
  ress@lbLabelBarOn         = False       ; turn off individual lb's 
 
  ress@cnLineLabelInterval = 1           ; Default is 2
  
 ; ress@pmTickMarkDisplayMode = "Always"         ; turn on nicer tickmarks
  ress@tmXBLabelFontHeightF  = 0.020
  ress@tmYLLabelFontHeightF  = 0.020
 ; ress@tmYLLabelStride       = 2                ; label every other tickmark
  ress@tmXBLabelStride       = 2
 
  ress@tmXTOn                  = True 
  ress@tmYROn                  = True 
  ;ress@tmYRLalabelsbelsOn            = False
  ress@tmXTLabelsOn            = False
 
  ress@lbOrientation        = "vertical"
  ress@lbTitleFontHeightF  = 0.020 
  ress@lbLabelFontHeightF  = 0.020
  ress@lbLeftMarginF       = 0.10
  ress@lbRightMarginF      = 0.55
  ;ress@lbLabelOffsetF	  = 0.20
  ress@lbLabelFont   	  = "Helvetica-Bold"     ; label font
  ress@lbTitleString       = "SSTA (~S~o~N~C)"
  ;ress@lbLabelPosition	  = "Top"
  ress@lbTitleAngleF       = 270.
  ress@lbTitlePosition      = "Right"
  ress@lbTitleDirection     = "Across"
  
  ress@lbLabelStride        = 10 
  
  ress@tmXBMode   = "Explicit"
  ress@tmXBValues = ispan(-50,10,10)
  ress@tmXBLabels = (/"50~S~o~N~W","40~S~o~N~W","30~S~o~N~W","20~S~o~N~W","10~S~o~N~W","0~S~o~N~","10~S~o~N~E","20~S~o~N~E"/)
  ;ress@tmXBMinorValues  = ispan(-50,20,5)
  ress@tmXBLabelFontHeightF = 0.025

  ress@tmYLMode = "Explicit"
  ress@tmYLValues = (/"0", "1","2","3","4","5","6","7","8","9","10","11"/)
  ress@tmYLLabels = (/"J", "F","M","A","M","J","J","A","S","O","N","D"/)
  ;ress@tmXBLabels = (/"-5", "-4","-3","-2","-1","0","+1","+2","+3","+4","+5","+6"/)
  ress@tiYAxisString = "Month"                  ; yaxis
  ress@tmYLLabelFontHeightF = 0.025
  
  ress@tmXTOn = False ; Turn off top tickmarks and labels
 ; ress@tmYROn = False ; Turn off left tickmarks and labels

  ress@tiXAxisFontHeightF = 0.025
  ress@tiYAxisFontHeightF = 0.025
  
  ress@tiYAxisFontThicknessF    = 1.
  ress@tmXBMajorThicknessF      = 1.
  ress@tmXBMinorThicknessF      = 1.
  ress@tmYLMajorThicknessF      = 1.
  ress@tmYLMinorThicknessF      = 1.
  ress@tmYLLabelFontThicknessF  = 1.
  ress@tmBorderThicknessF       = 1.
  ress@gsnLeftStringFontHeightF = 0.035
  ress@vpHeightF 	  = 0.75        ; Changes the aspect ratio
  ress@vpWidthF 	  = 0.80
;-----------------------------------------------------------------------------
;=============================================================================
;-----------------------------------------------------------------------------
  res1                      =   True
  res1@gsnDraw              = False        ; don't draw yet
  res1@gsnFrame             = False        ; don't advance frame yet
 ; res1@mpFillDrawOrder      = "preDraw"
  ;res1@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels
 ; res1@cnMinLevelValF       =  -1.0            ; set min contour level
 ; res1@cnMaxLevelValF       =  1.0             ; set max contour level
 ; res1@cnLevelSpacingF      =  .05            ; set contour spacing  
  res1@gsnAddCyclic         = True
  res1@cnLinesOn            = False       ; True is default
  res1@cnLineLabelsOn       = False        ; True is default
  res1@cnInfoLabelOn        = False
  res1@gsnRightString   = ""
  res1@gsnLeftString   = ""
  res1@cnLineLabelDensityF =  5
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
  opt = True
  opt@gsnShadeFillType = "pattern"        ; color is the default
  opt@gsnShadeLow = 17                     ; use fill pattern #17
  opt@gsnShadeHigh = 17                     ; use fill pattern #17
  opt@gsnShadeFillScaleF = 1.0
  opt@gsnShadeFillDotSizeF = 0.001
  

   delete(ress@cnFillPalette)
   ress@cnFillPalette        	= "amwg_blueyellowred"       ;different colour bar
   
   
   ;ress@gsnCenterString  = "Canonical Events"   
   ress@gsnLeftString  		= "Lag1"
   plotA0 = gsn_csm_contour(wks,d_hadiP6(month|:,{lon|-52:11}),ress)  ;  
   plotS0 = gsn_csm_contour(wks,d_hadidP1(:,:),res1)       ;
   plotS0 = gsn_contour_shade(plotS0,p95,d_hadidP1@_FillValue,opt)
   overlay(plotA0,plotS0)
   plot(0) = plotA0
   
   ress@gsnLeftString  		= "Lag2"
   plotA1 = gsn_csm_contour(wks,d_hadiP6111(month|:,{lon|-52:11}),ress)  ;           ; 
   plotS1 = gsn_csm_contour(wks,d_hadidP1111(:,:),res1)       ;
   plotS1 = gsn_contour_shade(plotS1,p95,d_hadidP1111@_FillValue,opt)
   overlay(plotA1,plotS1)
   plot(1) = plotA1
   
   ress@gsnLeftString  		= "Lag3"
   plotA0 = gsn_csm_contour(wks,d_hadiP62(month|:,{lon|-52:11}),ress)  ;  
   plotS0 = gsn_csm_contour(wks,d_hadidP12(:,:),res1)       ;
   plotS0 = gsn_contour_shade(plotS0,p95,d_hadidP12@_FillValue,opt)
   overlay(plotA0,plotS0)
   plot(2) = plotA0
   
   ress@gsnLeftString  		= "Lag4"
   plotA1 = gsn_csm_contour(wks,d_hadiP61112(month|:,{lon|-52:11}),ress)  ;           ; 
   plotS1 = gsn_csm_contour(wks,d_hadidP11112(:,:),res1)       ;
   plotS1 = gsn_contour_shade(plotS1,p95,d_hadidP11112@_FillValue,opt)
   overlay(plotA1,plotS1)
   plot(3) = plotA1
   
   ress@gsnLeftString  		= "Lag5"
   plotA0 = gsn_csm_contour(wks,d_hadiP63(month|:,{lon|-52:11}),ress)  ;  
   plotS0 = gsn_csm_contour(wks,d_hadidP13(:,:),res1)       ;
   plotS0 = gsn_contour_shade(plotS0,p95,d_hadidP13@_FillValue,opt)
   overlay(plotA0,plotS0)
   plot(4) = plotA0
   
   ress@gsnLeftString  		= "Lag6"
   plotA1 = gsn_csm_contour(wks,d_hadiP61113(month|:,{lon|-52:11}),ress)  ;           ; 
   plotS1 = gsn_csm_contour(wks,d_hadidP11113(:,:),res1)       ;
   plotS1 = gsn_contour_shade(plotS1,p95,d_hadidP11113@_FillValue,opt)
   overlay(plotA1,plotS1)
   plot(5) = plotA1

   ress@gsnLeftString  		= "Lag1"
   plotA0 = gsn_csm_contour(wks,d_hadiP64(month|:,{lon|-52:11}),ress)  ;  
   plotS0 = gsn_csm_contour(wks,d_hadidP14(:,:),res1)       ;
   plotS0 = gsn_contour_shade(plotS0,p95,d_hadidP14@_FillValue,opt)
   overlay(plotA0,plotS0)
   plot(6) = plotA0
   
   ress@gsnLeftString  		= "Lag2"
   plotA1 = gsn_csm_contour(wks,d_hadiP61114(month|:,{lon|-52:11}),ress)  ;           ; 
   plotS1 = gsn_csm_contour(wks,d_hadidP11114(:,:),res1)       ;
   plotS1 = gsn_contour_shade(plotS1,p95,d_hadidP11114@_FillValue,opt)
   overlay(plotA1,plotS1)
   plot(7) = plotA1
   
   ress@gsnLeftString  		= "Lag3"
   plotA0 = gsn_csm_contour(wks,d_hadiP65(month|:,{lon|-52:11}),ress)  ;  
   plotS0 = gsn_csm_contour(wks,d_hadidP15(:,:),res1)       ;
   plotS0 = gsn_contour_shade(plotS0,p95,d_hadidP15@_FillValue,opt)
   overlay(plotA0,plotS0)
   plot(8) = plotA0
   
   ress@gsnLeftString  		= "Lag4"
   plotA1 = gsn_csm_contour(wks,d_hadiP61115(month|:,{lon|-52:11}),ress)  ;           ; 
   plotS1 = gsn_csm_contour(wks,d_hadidP11115(:,:),res1)       ;
   plotS1 = gsn_contour_shade(plotS1,p95,d_hadidP11115@_FillValue,opt)
   overlay(plotA1,plotS1)
   plot(9) = plotA1
   
   ress@gsnLeftString  		= "Lag5"
   plotA0 = gsn_csm_contour(wks,d_hadiP66(month|:,{lon|-52:11}),ress)  ;  
   plotS0 = gsn_csm_contour(wks,d_hadidP16(:,:),res1)       ;
   plotS0 = gsn_contour_shade(plotS0,p95,d_hadidP16@_FillValue,opt)
   overlay(plotA0,plotS0)
   plot(10) = plotA0
   
   ress@gsnLeftString  		= "Lag6"
   plotA1 = gsn_csm_contour(wks,d_hadiP61116(month|:,{lon|-52:11}),ress)  ;           ; 
   plotS1 = gsn_csm_contour(wks,d_hadidP11116(:,:),res1)       ;
   plotS1 = gsn_contour_shade(plotS1,p95,d_hadidP11116@_FillValue,opt)
   overlay(plotA1,plotS1)
   plot(11) = plotA1
   
 delete(ress@cnMinLevelValF)
 delete(ress@cnMaxLevelValF)
 delete(ress@cnLevelSpacingF)
;------------------------------------------------------------------------------------------------------------------------------------------------------------------------
;=========================================================================
;=====================================================================
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
  
  txres5@txFontHeightF = 0.038   ; Use same font height as left axis
  amres5@amParallelPosF   =  -0.63 ; 0.5 is the right edge of the plot, so
  amres5@amOrthogonalPosF =  -0.55    ; This is the center of the plot.
  amres5@amJust           = "TopLeft"

 tt = (/"a","b","c","d","e","f","g","h","i","j","k","l"/)
create_tt = gsn_create_text(wks,tt,txres5);
; tt_p0 = new(nplots,graphic)                ; create graphic array ; only needed if paneling
tt_p0 = new(12,graphic)                ; create graphic array ; only needed if paneling
;---Attach the polyline
     do i = 0,11
	print(i)
        tt_p0(i) = gsn_add_annotation(plot(i),create_tt(i),amres5)
     end do
delete(i)
;=============================================================================
;========================================================================
 ;---Resources for the polyline
  pres1                  = True           ; polyline mods desired
  pres1@gsLineThicknessF = 2.0            ; line thickness   
  pres1@gsLineColor      = "black"         ; color of lines
  pres1@gsLineDashPattern = 2
  pres1@tfPolyDrawOrder   = "postDraw"
  pres1@gsLineColor      = "black"         ; color of lines
  ;dumL0 = gsn_add_polyline(wks,plotA0, (/x_0,x_1/),(/y_0,y_1/) ,pres1)
  
  x_2 = -60
  x_3 = 20
  y_2 = 5. 
  y_3 = 5.  
  
  x_0 = -60.
  x_1 = 20.
  y_0 = 8. 
  y_1 = 8.  
  dumL1 = gsn_add_polyline(wks,plot(0), (/x_2,x_3/),(/y_2,y_3/) ,pres1)
  dumL0 = gsn_add_polyline(wks,plot(1), (/x_2,x_3/),(/y_2,y_3/) ,pres1)
  
  dumL2 = gsn_add_polyline(wks,plot(2), (/x_2,x_3/),(/y_2,y_3/) ,pres1)
  dumL3 = gsn_add_polyline(wks,plot(3), (/x_2,x_3/),(/y_2,y_3/) ,pres1)
  dumL4 = gsn_add_polyline(wks,plot(4), (/x_2,x_3/),(/y_2,y_3/) ,pres1)
  dumL5 = gsn_add_polyline(wks,plot(5), (/x_2,x_3/),(/y_2,y_3/) ,pres1)
  dumL6 = gsn_add_polyline(wks,plot(6), (/x_2,x_3/),(/y_2,y_3/) ,pres1)
  dumL7 = gsn_add_polyline(wks,plot(7), (/x_2,x_3/),(/y_2,y_3/) ,pres1)
  dumL8 = gsn_add_polyline(wks,plot(8), (/x_2,x_3/),(/y_2,y_3/) ,pres1)
  dumL9 = gsn_add_polyline(wks,plot(9), (/x_2,x_3/),(/y_2,y_3/) ,pres1)
  dumL10 = gsn_add_polyline(wks,plot(10), (/x_2,x_3/),(/y_2,y_3/) ,pres1)
  dumL11 = gsn_add_polyline(wks,plot(11), (/x_2,x_3/),(/y_2,y_3/) ,pres1)
;========================================================
;==============================================================================================
  resP                                  = True         ; modify the panel plot
  resP@gsnDraw                          = True      ; don't draw yet
  resP@gsnFrame                         = False        ; don't advance frame yet
  resP@gsnPanelLabelBar    	        = True        ; add common colorbar
  resP@gsnPanelYWhiteSpacePercent      	= 6.0    ; default is 1.0
  resP@gsnPanelXWhiteSpacePercent      	= 3.0    ; default is 1.0
  resP@gsnPanelTop                      = 1.0
  resP@gsnPanelBottom                   = 0.25
  resP@gsnPanelLeft                     = 0.01
  resP@gsnPanelRight                    = 0.99
  resP@lbTitleFontHeightF	        = 0.010 
  resP@lbLabelFontHeightF               = 0.010
 ; resP@lbLabelStride		        = 5
  resP@lbLabelFont   		        = "Helvetica-Bold"     ; label font
  resP@lbTitleString                  	= "SST anomaly (~S~o~N~C)"
  ;resP@lbTitleString                  	= "SST trend (~S~o~N~C per 35-yr)"
  resP@lbTopMarginF                     = 0.5
  resP@gsnMaximize                      = True         ; large format
  resP@gsnPaperOrientation              = "portrait"   ; force portrait
  resP@gsnLeftString 			= ""
  resP@txString          		= ""

;**************************************************************
;**************************************************************
   gsn_panel(wks,plot(0:),(/2,6/),resP)     ; now draw as one plot
;=========================================================================
;=============================================================================
frame(wks)
end
