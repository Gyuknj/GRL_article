begin

;======================================
;======================================
  yrStrt = 1982
  yrLast = 2022
  
  nmos = 12	;months in a year
  
  p95 = 0.05

  mxlag = 7
  lag = 1	;first-order autocorrelation

  latS   = -15.
  latN   =  15. 
  lonL   = -60.
  lonR   =  20.
;==============================================================
  f   = addfile("/home/gyuk/Downloads/GYUK/GEOMAR/Final_work/HadISST_sst.nc", "r")   
  TIME   = f->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  hadi   = f->sst(iYYYY,:,:)
  hadi!1 = "lat"
  hadi!2 = "lon"
 ; printVarSummary(hadi)       ; variable overview
 ; printMinMax(hadi, True)
	
  delete(f)		;	Not needed anymore
  delete(TIME)
  delete(YYYY)
  delete(iYYYY)

  hadi = where(hadi.gt.-30.,hadi,hadi@_FillValue)	;removes padded values
  ;printVarSummary(hadi)       ; variable overview
  ;printMinMax(hadi,True)

;==========================
  hadiClm      = clmMonTLL(hadi)
  ;printVarSummary(hadiClm)                   ; (12,nlat,nlon)   
  hadi         = calcMonAnomTLL(hadi,hadiClm)    ; replace with anonamlies
  hadi@long_name = "ANOMALIES: "+hadi@long_name
  delete(hadiClm)
 ; printVarSummary(hadi)
;=========================
;----------------------------
  hadisst_d = new(dimsizes(hadi),typeof(hadi),hadi@_FillValue)
  x   = hadi&time
  x1  = x(::12)
  do tt = 0,nmos-1
    work =hadi(tt::nmos,:,:)
    work = dtrend_msg_n(x1,work,True,False,0)
    hadisst_d(tt::nmos,:,:) = work
  end do

 ; printVarSummary(hadisst_d)
  hadi = hadisst_d

  delete([/tt,x,x1,hadisst_d,work/])

  nyears = dimsizes(hadi&time)/nmos
  ;print(nyears)

  hadi_s     = onedtond(ndtooned(hadi),(/nyears,nmos,dimsizes(hadi&lat),dimsizes(hadi&lon)/))
  copy_VarCoords(hadi(0,:,:),hadi_s(0,0,:,:))                             ; copy coordinate arrays
  hadi_s!0 = "time"
  hadi_s&time = fspan(yrStrt,yrLast,nyears)
  hadi_s!1 = "month"
  printVarSummary(hadi_s)
  
  delete(hadi)
;============================================================== ============================================  
  hadi_s = dtrend_msg_n(hadi_s&time,hadi_s,True,False,0)        ;removes linear trends
  printVarSummary(hadi_s)
;============================================================== ============================================ 
;---------------------------------------------------------------------------------------------------------------------------------------------------
  ahadi_s_A = wgt_areaave_Wrap(hadi_s(time|:,month|3,{lon|-20:0}, {lat|-3:3}),1.0, 1.0, 0)        ;atl3 index
 ; printVarSummary(ahadi_s_A)
  ;printMinMax(ahadi_s_A,True)
  
  ahadi_s_A = dim_standardize_n(ahadi_s_A,0,0)      ; normalize 
 ; printVarSummary(ahadi_s_A)
 ; printMinMax(ahadi_s_A,True)
;---------------------------------------------------------------------------------------------------------------------------------------------------
  ahadi_s_M = wgt_areaave_Wrap(hadi_s(time|:,month|4,{lon|-20:0}, {lat|-3:3}),1.0, 1.0, 0)        ;atl3 index
 ; printVarSummary(ahadi_s_M)
  ;printMinMax(ahadi_s_M,True)
  
  ahadi_s_M = dim_standardize_n(ahadi_s_M,0,0)      ; normalize 
 ; printVarSummary(ahadi_s_M)
 ; printMinMax(ahadi_s_M,True)    
;------------------------------------------------------------------------------------------------------------------------------------
;---------------------------------------------------------------------------------------------------------------------------------------------------
  ahadi_all = wgt_areaave_Wrap(hadi_s(time|:,month|:,{lon|-20:0}, {lat|-3:3}),1.0, 1.0, 0)        ;atl3 index
 ; printVarSummary(ahadi_all)
  ;printMinMax(ahadi_all,True)
;------------------------------------------------------------------------------------------------------------------------------------
;------------------------------------------------------------------------------------------------------------------------------------
  ahadi_J = wgt_areaave_Wrap(hadi_s(time|:,month|5,{lon|-20:0}, {lat|-3:3}),1.0, 1.0, 0)        ;atl3 index
 ; printVarSummary(ahadi_J)
 ; printMinMax(ahadi_J,True)  
;============================================================== ============================================  
  ahadi_J = dim_standardize_n(ahadi_J,0,0)      ; normalize 
 ; printVarSummary(ahadi_J)
 ; printMinMax(ahadi_J,True)
;=================================================================================================================  
  d_SST = regCoef_n(ahadi_J,hadi_s(:,5,:,:),0,0)
  copy_VarCoords(hadi_s(0,0,:,:),d_SST)
 ; printVarSummary(d_SST)
;==============================================================
;==============================================================
;--------------test for significance
  d_SSTd = d_SST
 
  df = d_SST

  lag_ts = esacr(ahadi_J,mxlag) ; acr(nlat,nlon,mxlat+1)
  printVarSummary(lag_ts)
 
  lag_ts_map   = esacr(hadi_s(lat|:,lon|:,month|5,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
  printVarSummary(lag_ts_map)

  tval = onedtond(d_SST@tval,dimsizes(d_SST))
      
    ;do ne = 0,neof-1,1
           df(:,:) = floattointeger(dimsizes(hadi_s(:,0,0,0))*(1. -(lag_ts(lag)*lag_ts_map(:,:,lag)))\
        /(1. + (lag_ts_map(:,:,lag)*lag_ts(lag)))) 
    ;end do

   printVarSummary(tval)
   printVarSummary(df)
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   probs = betainc(df/(df+tval^2),df/2.0,b)       ; probU(nlat,nlon)
  
   probs@long_name   = "regression coefficient"
   probs@long_name = "probability"
   printVarSummary(probs)

  d_SSTd = where(probs.le.p95,d_SSTd,d_SSTd@_FillValue) ; Significant
  printVarSummary(d_SSTd)
  printMinMax(probs,True)

  delete([/tval,df,b,probs/])
  delete(lag_ts)
  delete(lag_ts_map)
  delete(nyears)
;=================================================================================================================================
;==============================================================
  d_SST_A = regCoef_n(ahadi_s_A,hadi_s(:,3,:,:),0,0)
  copy_VarCoords(hadi_s(0,0,:,:),d_SST_A)
  printVarSummary(d_SST_A)
;==============================================================
;==============================================================
;--------------test for significance
  d_SSTd_A = d_SST_A

  df = d_SST_A

  lag_ts = esacr(ahadi_s_A,mxlag) ; acr(nlat,nlon,mxlat+1)
  printVarSummary(lag_ts)

  lag_ts_map   = esacr(hadi_s(lat|:,lon|:,month|3,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
  printVarSummary(lag_ts_map)

  tval = onedtond(d_SST_A@tval,dimsizes(d_SST_A))
      
    ;do ne = 0,neof-1,1
           df(:,:) = floattointeger(dimsizes(hadi_s(:,0,0,0))*(1. -(lag_ts(lag)*lag_ts_map(:,:,lag)))\
        /(1. + (lag_ts_map(:,:,lag)*lag_ts(lag)))) 
    ;end do

   printVarSummary(tval)
   printVarSummary(df)
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   probs = betainc(df/(df+tval^2),df/2.0,b)       ; probU(nlat,nlon)
  
   probs@long_name   = "regression coefficient"
   probs@long_name = "probability"
   printVarSummary(probs)

  d_SSTd_A = where(probs.le.p95,d_SSTd_A,d_SSTd_A@_FillValue) ; Significant
  printVarSummary(d_SSTd_A)
  printMinMax(probs,True)

  delete([/tval,df,b,probs/])
  delete(lag_ts)
  delete(lag_ts_map)
;=================================================================================================================================
;=================================================================================================================================
;==============================================================
  d_SST_M = regCoef_n(ahadi_s_M,hadi_s(:,4,:,:),0,0)
  copy_VarCoords(hadi_s(0,0,:,:),d_SST_M)
;  printVarSummary(d_SST_M)
;==============================================================
;==============================================================
;--------------test for significance
  d_SSTd_M = d_SST_M

  df = d_SST_M

  lag_ts = esacr(ahadi_s_M,mxlag) ; acr(nlat,nlon,mxlat+1)
  printVarSummary(lag_ts)

  lag_ts_map   = esacr(hadi_s(lat|:,lon|:,month|4,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
  printVarSummary(lag_ts_map)

  tval = onedtond(d_SST_A@tval,dimsizes(d_SST_M))
      
    ;do ne = 0,neof-1,1
           df(:,:) = floattointeger(dimsizes(hadi_s(:,0,0,0))*(1. -(lag_ts(lag)*lag_ts_map(:,:,lag)))\
        /(1. + (lag_ts_map(:,:,lag)*lag_ts(lag)))) 
    ;end do

   printVarSummary(tval)
   printVarSummary(df)
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   probs = betainc(df/(df+tval^2),df/2.0,b)       ; probU(nlat,nlon)
  
   probs@long_name   = "regression coefficient"
   probs@long_name = "probability"
   printVarSummary(probs)

  d_SSTd_M = where(probs.le.p95,d_SSTd_M,d_SSTd_M@_FillValue) ; Significant
 ; printVarSummary(d_SSTd_M)
 ; printMinMax(probs,True)

  delete([/tval,df,b,probs/])
  delete(lag_ts)
  delete(lag_ts_map)
  delete([/hadi_s,ahadi_all/])
;=================================================================================================================================
;==============================================================
;==============================================================
  
  fwnd1   = addfile("/home/gyuk/Downloads/GYUK/GEOMAR/Final_work/regrid_taux.nc", "r")   
  TIME   = fwnd1->time_counter
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  gpcc1   = fwnd1->sozotaux(iYYYY,:,:)
 ; gpcc1!0 = "time"
  gpcc1 = lonFlip(gpcc1)
  printVarSummary(gpcc1)       ; variable overview
 ; printMinMax(gpcc1,True)
  
  delete(fwnd1)             ;       Not needed anymore
  delete(TIME)
  delete(YYYY)
  delete(iYYYY)
;==========================
;----------------------------
  nyears = dimsizes(gpcc1&time_counter)/nmos
  ;print(nyears)
  
  gpcc1_d = new(dimsizes(gpcc1),typeof(gpcc1),gpcc1@_FillValue)
  x   = gpcc1&time_counter
  x1  = x(::12)
  do tt = 0,nmos-1
     work =gpcc1(tt::nmos,:,:)
     work = dtrend_msg_n(x1,work,False,False,0)
     gpcc1_d(tt::nmos,:,:) = work
  end do

  printVarSummary(gpcc1_d)
  gpcc1 = gpcc1_d
;===
  delete([/tt,x,x1,gpcc1_d,work/])
;==========================
;----------------------------
  nyears = dimsizes(gpcc1&time_counter)/nmos
  ;print(nyears) 

  gpcc1_s     = onedtond(ndtooned(gpcc1),(/nyears,nmos,dimsizes(gpcc1&lat),dimsizes(gpcc1&lon)/))
  copy_VarCoords(gpcc1(0,:,:),gpcc1_s(0,0,:,:))                             ; copy coordinate arrays
  gpcc1_s!0 = "time_counter"
  gpcc1_s&time_counter = fspan(yrStrt,yrLast,nyears)
  gpcc1_s!1 = "month"
 ; printVarSummary(gpcc1_s)
  
  delete(gpcc1)
  
  gpcc1_s = dtrend_msg_n(gpcc1_s&time_counter,gpcc1_s,True,False,0)        ;removes linear trends
  printVarSummary(gpcc1_s)
  ;printMinMax(gpcc1_s,True)
;==========================================================================================
;==========================================================================================
  d_USTR = regCoef_n(ahadi_J,gpcc1_s(:,:,:,:),0,0)
  copy_VarCoords(gpcc1_s(0,:,:,:),d_USTR)
  printVarSummary(d_USTR)
;==============================================================
;==============================================================
;--------------test for significance
  d_USTRd = d_USTR
 
  df = d_USTR

  lag_ts = esacr(ahadi_J,mxlag) ; acr(nlat,nlon,mxlat+1)
  printVarSummary(lag_ts)
 
  lag_ts_map   = esacr(gpcc1_s(month|:,lat|:,lon|:,time_counter|:),mxlag) ; acr(nlat,nlon,mxlat+1)
  printVarSummary(lag_ts_map)

  tval = onedtond(d_USTR@tval,dimsizes(d_USTR))
      
    ;do ne = 0,neof-1,1
           df(:,:,:) = floattointeger(dimsizes(gpcc1_s(:,0,0,0))*(1. -(lag_ts(lag)*lag_ts_map(:,:,:,lag)))\
        /(1. + (lag_ts_map(:,:,:,lag)*lag_ts(lag)))) 
    ;end do
   
   printVarSummary(tval)
   printVarSummary(df)
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   probU = betainc(df/(df+tval^2),df/2.0,b)       ; probU(nlat,nlon)


   
   probU@long_name   = "regression coefficient"
   probU@long_name = "probability"
   printVarSummary(probU)

  d_USTRd = where(probU.le.p95,d_USTRd,d_USTRd@_FillValue) ; Significant
  printVarSummary(d_USTRd)
  printMinMax(probU,True)

  delete([/tval,df,b,probU/])
  delete(lag_ts)
  delete(lag_ts_map)
 ; delete(gpcc1_s)
 ; delete(gpcc1)
 ; delete(nyears)
;=================================================================================================================================
;==========================================================================================
  d_USTR_A = regCoef_n(ahadi_s_A,gpcc1_s(:,:,:,:),0,0)
  copy_VarCoords(gpcc1_s(0,:,:,:),d_USTR_A)
  printVarSummary(d_USTR)
;==============================================================
;==============================================================
;--------------test for significance
  d_USTRd_A = d_USTR_A
 
  df = d_USTR_A

  lag_ts = esacr(ahadi_s_A,mxlag) ; acr(nlat,nlon,mxlat+1)
  printVarSummary(lag_ts)
 
  lag_ts_map   = esacr(gpcc1_s(month|:,lat|:,lon|:,time_counter|:),mxlag) ; acr(nlat,nlon,mxlat+1)
  printVarSummary(lag_ts_map)

  tval = onedtond(d_USTR_A@tval,dimsizes(d_USTR_A))
      
    ;do ne = 0,neof-1,1
           df(:,:,:) = floattointeger(dimsizes(gpcc1_s(:,0,0,0))*(1. -(lag_ts(lag)*lag_ts_map(:,:,:,lag)))\
        /(1. + (lag_ts_map(:,:,:,lag)*lag_ts(lag)))) 
    ;end do
   
   printVarSummary(tval)
   printVarSummary(df)
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   probU = betainc(df/(df+tval^2),df/2.0,b)       ; probU(nlat,nlon)


   
   probU@long_name   = "regression coefficient"
   probU@long_name = "probability"
   printVarSummary(probU)

  d_USTRd_A = where(probU.le.p95,d_USTRd_A,d_USTRd_A@_FillValue) ; Significant
  printVarSummary(d_USTRd_A)
  printMinMax(probU,True)

  delete([/tval,df,b,probU/])
  delete(lag_ts)
  delete(lag_ts_map)
;==========================================================================================
;==========================================================================================
  d_USTR_M = regCoef_n(ahadi_s_M,gpcc1_s(:,:,:,:),0,0)
  copy_VarCoords(gpcc1_s(0,:,:,:),d_USTR_M)
  printVarSummary(d_USTR_M)
;==============================================================
;==============================================================
;--------------test for significance
  d_USTRd_M = d_USTR_M
 
  df = d_USTR_M

  lag_ts = esacr(ahadi_s_M,mxlag) ; acr(nlat,nlon,mxlat+1)
  printVarSummary(lag_ts)
 
  lag_ts_map   = esacr(gpcc1_s(month|:,lat|:,lon|:,time_counter|:),mxlag) ; acr(nlat,nlon,mxlat+1)
  printVarSummary(lag_ts_map)

  tval = onedtond(d_USTR_M@tval,dimsizes(d_USTR_M))
      
    ;do ne = 0,neof-1,1
           df(:,:,:) = floattointeger(dimsizes(gpcc1_s(:,0,0,0))*(1. -(lag_ts(lag)*lag_ts_map(:,:,:,lag)))\
        /(1. + (lag_ts_map(:,:,:,lag)*lag_ts(lag)))) 
    ;end do
   
   printVarSummary(tval)
   printVarSummary(df)
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   probU = betainc(df/(df+tval^2),df/2.0,b)       ; probU(nlat,nlon)


   
   probU@long_name   = "regression coefficient"
   probU@long_name = "probability"
   printVarSummary(probU)

  d_USTRd_M = where(probU.le.p95,d_USTRd_M,d_USTRd_M@_FillValue) ; Significant
  printVarSummary(d_USTRd_M)
  printMinMax(probU,True)

  delete([/tval,df,b,probU,gpcc1_s/])
  delete(lag_ts)
  delete(lag_ts_map)
;==========================================================================================
;=================================================================================================================================
;=================================================================================================================================

  fwnd1    = addfile("/home/gyuk/Downloads/GYUK/GEOMAR/Final_work/regrid_tauy.nc", "r")   
  TIME     = fwnd1->time_counter
  YYYY     = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY    = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  gpcc1v   = fwnd1->sometauy(iYYYY,:,:)
  gpcc1v!0 = "time"
  gpcc1v = lonFlip(gpcc1v)
  printVarSummary(gpcc1v)       ; variable overview
  printMinMax(gpcc1v,True)
  
  delete(fwnd1)             ;       Not needed anymore
  delete(TIME)
  delete(YYYY)
  delete(iYYYY)
;==========================
;----------------------------
  nyears = dimsizes(gpcc1v&time)/nmos
  ;print(nyears)
  
  gpcc1v_d = new(dimsizes(gpcc1v),typeof(gpcc1v),gpcc1v@_FillValue)
  x   = gpcc1v&time
  x1  = x(::12)
  do tt = 0,nmos-1
     work =gpcc1v(tt::nmos,:,:)
     work = dtrend_msg_n(x1,work,False,False,0)
     gpcc1v_d(tt::nmos,:,:) = work
  end do

  printVarSummary(gpcc1v_d)
  gpcc1v = gpcc1v_d
;===
  delete([/tt,x,x1,gpcc1v_d,work/])
;==========================
;----------------------------
  nyears = dimsizes(gpcc1v&time)/nmos
  ;print(nyears) 
 
  gpcc1v_s     = onedtond(ndtooned(gpcc1v),(/nyears,nmos,dimsizes(gpcc1v&lat),dimsizes(gpcc1v&lon)/))
  copy_VarCoords(gpcc1v(0,:,:),gpcc1v_s(0,0,:,:))                             ; copy coordinate arrays
  gpcc1v_s!0 = "time"
  gpcc1v_s&time = fspan(yrStrt,yrLast,nyears)
  gpcc1v_s!1 = "month"
  printVarSummary(gpcc1v_s)

  gpcc1v_s = dtrend_msg_n(gpcc1v_s&time,gpcc1v_s,True,False,0)        ;removes linear trends
  printVarSummary(gpcc1v_s)
  ;printMinMax(gpcc1v_s,True)
  
;==========================================================================================
;==========================================================================================
  d_VSTR = regCoef_n(ahadi_J,gpcc1v_s(:,:,:,:),0,0)
  copy_VarCoords(gpcc1v_s(0,:,:,:),d_VSTR)
  printVarSummary(d_VSTR)
;==============================================================
;==============================================================
;--------------test for significance
  d_VSTRd = d_VSTR
 
  df = d_VSTR

  lag_ts = esacr(ahadi_J,mxlag) ; acr(nlat,nlon,mxlat+1)
  printVarSummary(lag_ts)
 
  lag_ts_map   = esacr(gpcc1v_s(month|:,lat|:,lon|:,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
  printVarSummary(lag_ts_map)

  tval = onedtond(d_VSTR@tval,dimsizes(d_VSTR))
      
    ;do ne = 0,neof-1,1
           df(:,:,:) = floattointeger(dimsizes(gpcc1v_s(:,0,0,0))*(1. -(lag_ts(lag)*lag_ts_map(:,:,:,lag)))\
        /(1. + (lag_ts_map(:,:,:,lag)*lag_ts(lag)))) 
    ;end do
   
   printVarSummary(tval)
   printVarSummary(df)
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   probV = betainc(df/(df+tval^2),df/2.0,b)       ; probU(nlat,nlon)
 
   probV@long_name   = "regression coefficient"
   probV@long_name = "probability"
   printVarSummary(probV)

  d_VSTRd = where(probV.le.p95,d_VSTRd,d_VSTRd@_FillValue) ; Significant
  printVarSummary(d_VSTRd)
  printMinMax(probV,True)

  delete([/tval,df,b,probV/])
  delete(lag_ts)
  delete(lag_ts_map)
  delete(gpcc1v)
  delete(nyears)
;==========================================================================================
  d_VSTR_A = regCoef_n(ahadi_s_A,gpcc1v_s(:,:,:,:),0,0)
  copy_VarCoords(gpcc1v_s(0,:,:,:),d_VSTR_A)
  printVarSummary(d_VSTR_A)
;==============================================================
;==============================================================
;--------------test for significance
  d_VSTRd_A = d_VSTR_A
 
  df = d_VSTR_A

  lag_ts = esacr(ahadi_s_A,mxlag) ; acr(nlat,nlon,mxlat+1)
  printVarSummary(lag_ts)
 
  lag_ts_map   = esacr(gpcc1v_s(month|:,lat|:,lon|:,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
  printVarSummary(lag_ts_map)

  tval = onedtond(d_VSTR_A@tval,dimsizes(d_VSTR_A))
      
    ;do ne = 0,neof-1,1
           df(:,:,:) = floattointeger(dimsizes(gpcc1v_s(:,0,0,0))*(1. -(lag_ts(lag)*lag_ts_map(:,:,:,lag)))\
        /(1. + (lag_ts_map(:,:,:,lag)*lag_ts(lag)))) 
    ;end do
   
   printVarSummary(tval)
   printVarSummary(df)
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   probV = betainc(df/(df+tval^2),df/2.0,b)       ; probU(nlat,nlon)
 
   probV@long_name   = "regression coefficient"
   probV@long_name = "probability"
   printVarSummary(probV)

  d_VSTRd_A = where(probV.le.p95,d_VSTRd_A,d_VSTRd_A@_FillValue) ; Significant
  printVarSummary(d_VSTRd_A)
  printMinMax(probV,True)

  delete([/tval,df,b,probV/])
  delete(lag_ts)
  delete(lag_ts_map)
;==========================================================================================
;==========================================================================================
  d_VSTR_M = regCoef_n(ahadi_s_M,gpcc1v_s(:,:,:,:),0,0)
  copy_VarCoords(gpcc1v_s(0,:,:,:),d_VSTR_M)
  printVarSummary(d_VSTR_M)
;==============================================================
;==============================================================
;--------------test for significance
  d_VSTRd_M = d_VSTR_M
 
  df = d_VSTR_M

  lag_ts = esacr(ahadi_s_M,mxlag) ; acr(nlat,nlon,mxlat+1)
  printVarSummary(lag_ts)
 
  lag_ts_map   = esacr(gpcc1v_s(month|:,lat|:,lon|:,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
  printVarSummary(lag_ts_map)

  tval = onedtond(d_VSTR_M@tval,dimsizes(d_VSTR_M))
      
    ;do ne = 0,neof-1,1
           df(:,:,:) = floattointeger(dimsizes(gpcc1v_s(:,0,0,0))*(1. -(lag_ts(lag)*lag_ts_map(:,:,:,lag)))\
        /(1. + (lag_ts_map(:,:,:,lag)*lag_ts(lag)))) 
    ;end do
   
   printVarSummary(tval)
   printVarSummary(df)
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   probV = betainc(df/(df+tval^2),df/2.0,b)       ; probU(nlat,nlon)
 
   probV@long_name   = "regression coefficient"
   probV@long_name = "probability"
   printVarSummary(probV)

  d_VSTRd_M = where(probV.le.p95,d_VSTRd_M,d_VSTRd_M@_FillValue) ; Significant
  printVarSummary(d_VSTRd_M)
  printMinMax(probV,True)

  delete([/tval,df,b,probV/])
  delete(lag_ts)
  delete(lag_ts_map)
;=============================================================================================================

  fssh   = addfile("/home/gyuk/Downloads/GYUK/GEOMAR/Final_work/sea_surface_height.nc", "r")   
  TIME   = fssh->time
  YYYY   = cd_calendar(TIME,-1)/100  
  iYYYY  = ind(YYYY.ge.yrStrt.and.YYYY.le.yrLast)
  ssh    = fssh->sshg(iYYYY,:,:)
 
  ssh = ssh*100.
  printVarSummary(ssh)       ; variable overview

  printVarSummary(ssh)               ; entire file
  printMinMax(ssh,True)
;==============================================================
; compute climatology and Anomalies
;==============================================================
  sshClm   = clmMonTLL(ssh)                             ; (12,lat,lon)
  printVarSummary(sshClm)

  ssh  = calcMonAnomTLL(ssh,sshClm)                ; (time, lat,lon)
  printVarSummary(ssh)         
  printMinMax(ssh,True)
;=================================================================
   ssh = lonFlip(ssh)
   
   delete(fssh)
   delete(TIME) 
   delete(YYYY)
   delete(iYYYY)            ; Re-set the time dimension
;=========================================================================================
; Calculate the linear trends and remove the linear trends
;=============================================================================================================
;=============================================================================================================
  nyears = dimsizes(ssh&time)/nmos
  print(nyears)
  ssh_d = new(dimsizes(ssh),typeof(ssh),ssh@_FillValue)
  x   = ssh&time
  x1  = x(::12)
  do tt = 0,nmos-1
    work =ssh(tt::nmos,:,:)
    work = dtrend_msg_n(x1,work,True,False,0)
    ssh_d(tt::nmos,:,:) = work
  end do

  printVarSummary(ssh_d)
  ssh = ssh_d

;====================================================================================================================
  delete([/tt,x,x1,ssh_d,work/])
;=========================================================================================
;=========================================================================================
  nyears = dimsizes(ssh&time)/nmos
  print(nyears)

  ssh_s     = onedtond(ndtooned(ssh),(/nyears,nmos,dimsizes(ssh&lat),dimsizes(ssh&lon)/))
  printVarSummary(ssh_s)

  copy_VarCoords(ssh(0,:,:),ssh_s(0,0,:,:))                             ; copy coordinate arrays from slp to newslp
  ssh_s!0 = "time"
  ssh_s!1 = "month"
  printVarSummary(ssh_s)
;==========================================================================================
;==========================================================================================
  d_SSH = regCoef_n(ahadi_J,ssh_s(:,5,:,:),0,0)
  copy_VarCoords(ssh_s(0,0,:,:),d_SSH)
  printVarSummary(d_SSH)
;==============================================================
;--------------test for significance
 d_SSHd = d_SSH
 
 df = d_SSH

 lag_ts = esacr(ahadi_J,mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts)
 
 lag_ts_map   = esacr(ssh_s(lat|:,lon|:,month|5,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts_map)

 tval = onedtond(d_SSH@tval,dimsizes(d_SSH))
      
    ;do ne = 0,neof-1,1
           df(:,:) = floattointeger(dimsizes(ssh_s(:,0,0,0))*(1. -(lag_ts(lag)*lag_ts_map(:,:,lag)))\
        /(1. + (lag_ts_map(:,:,lag)*lag_ts(lag)))) 
    ;end do
   
   printVarSummary(tval)
   printVarSummary(df)
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)       ; prob(nlat,nlon)


   
   prob@long_name   = "regression coefficient"
   prob@long_name = "probability"
   printVarSummary(prob)

 d_SSHd = where(prob.le.p95,d_SSHd,d_SSHd@_FillValue) ; Significant
 printVarSummary(d_SSHd)
 printMinMax(prob,True)

 delete([/tval,df,prob,b/])
 delete(lag_ts)
 delete(lag_ts_map)
 delete(nyears)
;=============================================================================================================
;==========================================================================================
  d_SSH_A = regCoef_n(ahadi_s_A,ssh_s(:,3,:,:),0,0)
  copy_VarCoords(ssh_s(0,0,:,:),d_SSH_A)
  printVarSummary(d_SSH_A)
;==============================================================
;--------------test for significance
 d_SSHd_A = d_SSH_A
 
 df = d_SSH_A

 lag_ts = esacr(ahadi_s_A,mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts)
 
 lag_ts_map   = esacr(ssh_s(lat|:,lon|:,month|3,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts_map)

 tval = onedtond(d_SSH_A@tval,dimsizes(d_SSH_A))
      
    ;do ne = 0,neof-1,1
           df(:,:) = floattointeger(dimsizes(ssh_s(:,0,0,0))*(1. -(lag_ts(lag)*lag_ts_map(:,:,lag)))\
        /(1. + (lag_ts_map(:,:,lag)*lag_ts(lag)))) 
    ;end do
   
   printVarSummary(tval)
   printVarSummary(df)
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)       ; prob(nlat,nlon)


   
   prob@long_name   = "regression coefficient"
   prob@long_name = "probability"
   printVarSummary(prob)

 d_SSHd_A = where(prob.le.p95,d_SSHd_A,d_SSHd_A@_FillValue) ; Significant
 printVarSummary(d_SSHd_A)
 printMinMax(prob,True)

 delete([/tval,df,prob,b/])
 delete(lag_ts)
 delete(lag_ts_map)
 ;delete(nyears)
;=============================================================================================================
;==========================================================================================
  d_SSH_M = regCoef_n(ahadi_s_M,ssh_s(:,4,:,:),0,0)
  copy_VarCoords(ssh_s(0,0,:,:),d_SSH_M)
  printVarSummary(d_SSH_M)
;==============================================================
;--------------test for significance
 d_SSHd_M = d_SSH_M
 
 df = d_SSH_M

 lag_ts = esacr(ahadi_s_M,mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts)
 
 lag_ts_map   = esacr(ssh_s(lat|:,lon|:,month|4,time|:),mxlag) ; acr(nlat,nlon,mxlat+1)
 printVarSummary(lag_ts_map)

 tval = onedtond(d_SSH_M@tval,dimsizes(d_SSH_M))
      
    ;do ne = 0,neof-1,1
           df(:,:) = floattointeger(dimsizes(ssh_s(:,0,0,0))*(1. -(lag_ts(lag)*lag_ts_map(:,:,lag)))\
        /(1. + (lag_ts_map(:,:,lag)*lag_ts(lag)))) 
    ;end do
   
   printVarSummary(tval)
   printVarSummary(df)
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)       ; prob(nlat,nlon)


   
   prob@long_name   = "regression coefficient"
   prob@long_name = "probability"
   printVarSummary(prob)

 d_SSHd_A = where(prob.le.p95,d_SSHd_M,d_SSHd_M@_FillValue) ; Significant
 printVarSummary(d_SSHd_M)
 printMinMax(prob,True)

 delete([/tval,df,prob,b/])
 delete(lag_ts)
 delete(lag_ts_map)
 ;delete(nyears)
;=================================================================================================================================
;=============================================================================================================
  
  f   = addfile("/home/gyuk/Downloads/GYUK/GEOMAR/Final_work/precip.mon.mean.nc", "r")   
  TIME   = f->time
  YYYY   = cd_calendar(TIME,-1)/100  
  iYYYY  = ind(YYYY.ge.yrStrt.and.YYYY.le.yrLast)
  precip    = f->precip(iYYYY,:,:)
 
  printVarSummary(precip)               ; entire file
  printMinMax(precip,True)
;==============================================================
; compute climatology and Anomalies
;==============================================================
  precipClm   = clmMonTLL(precip)                             ; (12,lat,lon)
  printVarSummary(precipClm)

  precip  = calcMonAnomTLL(precip,precipClm)                ; (time, lat,lon)
  printVarSummary(precip)         
  printMinMax(precip,True)
;=================================================================
;=================================================================
   precip = lonFlip(precip)
   
   delete(f)
   delete(TIME) 
   delete(YYYY)
   delete(iYYYY)            ; Re-set the time dimension
;=========================================================================================
; Calculate the linear trends and remove the linear trends
;=============================================================================================================
;=============================================================================================================
  nyears = dimsizes(precip&time)/nmos
  print(nyears)
  precip_d = new(dimsizes(precip),typeof(precip),precip@_FillValue)
  x   = precip&time
  x1  = x(::12)
  do tt = 0,nmos-1
    work = precip(tt::nmos,:,:)
    work = dtrend_msg_n(x1,work,True,False,0)
    precip_d(tt::nmos,:,:) = work
  end do

  printVarSummary(precip_d)
  precip = precip_d

  delete([/tt,x,x1,precip_d,work/])
;=========================================================================================
;=========================================================================================
  nyears = dimsizes(precip&time)/nmos
  print(nyears)

  precip_s     = onedtond(ndtooned(precip),(/nyears,nmos,dimsizes(precip&lat),dimsizes(precip&lon)/))
  printVarSummary(precip_s)

  copy_VarCoords(precip(0,:,:),precip_s(0,0,:,:))                             ; copy coordinate arrays from slp to newslp
  precip_s!0 = "time"
  precip_s!1 = "month"
  printVarSummary(precip_s)
;==========================================================================================0,1,2,3,4,5,6
;==========================================================================================
  d_precip = regCoef_n(ahadi_J,precip_s(:,5,:,:),0,0)
  copy_VarCoords(precip_s(0,0,:,:),d_precip)
  printVarSummary(d_precip)
;==============================================================
;==========================================================================================
  d_precip_A = regCoef_n(ahadi_s_A,precip_s(:,3,:,:),0,0)
  copy_VarCoords(precip_s(0,0,:,:),d_precip_A)
  printVarSummary(d_precip_A)
;==============================================================
  d_precip_M = regCoef_n(ahadi_s_M,precip_s(:,4,:,:),0,0)
  copy_VarCoords(precip_s(0,0,:,:),d_precip_M)
  printVarSummary(d_precip_M)
;==============================================================
;===============================================================================================
;===============================================================================================
  wks = gsn_open_wks("eps","Figure_S1")         ; send graphics to PNG file
  plot = new(6,graphic)                ; create graphic array ; only needed if paneling
  
  res                      = True         
  res@gsnDraw              = False        ; don't draw yet
  res@gsnFrame             = False        ; don't advance frame yet
  res@gsnAddCyclic             = False             ; regional data 
   
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
  ;res@pmLabelBarParallelPosF = 0.50
  ;res@pmLabelBarOrthogonalPosF = -0.11
  res@lbOrientation        = "vertical"
  res@lbTitleFontHeightF  = 0.012 
  res@lbLabelFontHeightF  = 0.012
  res@lbLeftMarginF       = 0.10
  res@lbRightMarginF      = 0.55
  ;res@lbLabelOffsetF	  = 0.20
  res@lbLabelFont   	  = "Helvetica-Bold"     ; label font
  res@lbTitleString       = "SSH response (~F10~e~S~~F21~-2~N~m)"
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
  res@pmTickMarkDisplayMode = "Never"
  ;res@tmXBLabelsOn         = False      ;no lon ticks
  ;res@tmYLLabelsOn         = False      ;no lon ticks
  res@tmYLMode   = "Explicit"
  res@tmYLValues = ispan(-15,15,15)
  res@tmYLLabels = (/"15~S~o~N~S", "0~S~o~N","15~S~o~N~N"/)
  res@tmYLMinorValues  = ispan(-15,15,5) 
  res@tmYLLabelFontHeightF = 0.017

  res@tmXBMode   = "Explicit"
  res@tmXBValues = ispan(-60,20,20)
  res@tmXBLabels = (/"60~S~o~N~W","40~S~o~N~W","20~S~o~N~W","0~S~o~N~","20~S~o~N~E"/)
  ;res@tmXBMinorValues  = ispan(-50,20,5)
  res@tmXBLabelFontHeightF = 0.017 
 
  res@mpGeophysicalLineColor = "grey" 
  res@mpFillDrawOrder      = "postDraw"

  res@vpHeightF 	  = 0.80        ; Changes the aspect ratio
  res@vpWidthF 		  = 0.80
  res@gsnCenterStringFontHeightF   = 0.030
  
  res@gsnCenterString  = ""   
  res@gsnLeftString  = ""  
  res@gsnRightString  = ""  
    
  res@vcRefAnnoOrthogonalPosF = -1.0             ; move ref vector up
  res@vcRefMagnitudeF         = 0.01             ; define vector ref mag
  res@vcRefAnnoString1          = "0.01 N m~S~-2~N~ "
  res@vcRefLengthF            = 0.04            ; define length of vec ref
 ; vcres@vcGlyphStyle            = "CurlyVector"    ; turn on curly vectors
  res@vcMinDistanceF          = 0.030
  res@vcLineArrowThicknessF   = 1.
  res@vcLineArrowColor        = "blue"           ; change vector color
  res@gsnScalarContour	      = True
  res@vcLineArrowHeadMaxSizeF = 0.05
  res@vcLineArrowHeadMinSizeF = 0.01
 ; res@vcMaxLevelCount         = 100
  
  res@tiYAxisFontThicknessF    = 1.
  res@tmXBMajorThicknessF      = 1.
  res@tmXBMinorThicknessF      = 1.
  res@tmYLMajorThicknessF      = 1.
  res@tmYLMinorThicknessF      = 1.
  res@tmYLLabelFontThicknessF  = 1.
  res@tmBorderThicknessF       = 1. 
;-----------------------------------------------------------------------------
  resn                      = True         
  resn@gsnDraw              = False        ; don't draw yet
  resn@gsnFrame             = False        ; don't advance frame yet
  resn@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels
  resn@cnMinLevelValF       = -1.            ; set min contour level
  resn@cnMaxLevelValF       =  1.             ; set max contour level
  resn@cnLevelSpacingF      =  .1            ; set contour spacing   
  resn@cnFillOn             = True         ; turn on color fill
  resn@cnFillPalette        = "cmp_b2r"  
 
  resn@cnLinesOn            = False        ; True is default
  resn@cnLineLabelsOn       = False        ; True is default
  resn@cnInfoLabelOn        = False
  resn@lbLabelBarOn         = True       ; turn off individual lb's 
  ;resn@pmLabelBarParallelPosF = 0.50
  ;resn@pmLabelBarOrthogonalPosF = -0.11
  resn@lbOrientation        = "vertical"
  resn@lbTitleFontHeightF  = 0.012 
  resn@lbLabelFontHeightF  = 0.012
  resn@lbLeftMarginF       = 0.10
  resn@lbRightMarginF      = 0.55
  ;resn@lbLabelOffsetF	  = 0.20
  resn@lbLabelFont   	  = "Helvetica-Bold"     ; label font
  resn@lbTitleString       = "SST response (~S~o~N~C)"
  ;resn@lbLabelPosition	  = "Top"
  resn@lbTitleAngleF       = 270.
  resn@lbTitlePosition      = "Right"
  resn@lbTitleDirection     = "Across"
   
  resn@lbLabelStride        = 10 
  
  resn@mpFillOn             = True        ; turn off map fill
  resn@mpMinLatF            = latS         ; zoom in on map
  resn@mpMaxLatF            = latN
  resn@mpMinLonF            = lonL
  resn@mpMaxLonF            = lonR
  resn@pmTickMarkDisplayMode = "Never"
  ;resn@tmXBLabelsOn         = False      ;no lon ticks
  ;resn@tmYLLabelsOn         = False      ;no lon ticks
  resn@tmYLMode   = "Explicit"
  resn@tmYLValues = ispan(-15,15,15)
  resn@tmYLLabels = (/"15~S~o~N~S", "0~S~o~N","15~S~o~N~N"/)
  resn@tmYLMinorValues  = ispan(-15,15,5) 
  resn@tmYLLabelFontHeightF = 0.017

  resn@tmXBMode   = "Explicit"
  resn@tmXBValues = ispan(-60,20,20)
  resn@tmXBLabels = (/"60~S~o~N~W","40~S~o~N~W","20~S~o~N~W","0~S~o~N~","20~S~o~N~E"/)
  ;resn@tmXBMinorValues  = ispan(-50,20,5)
  resn@tmXBLabelFontHeightF = 0.017 
 
  resn@mpGeophysicalLineColor = "grey" 
  resn@mpFillDrawOrder      = "postDraw"

  resn@vpHeightF 	  = 0.80        ; Changes the aspect ratio
  resn@vpWidthF 		  = 0.80
  resn@gsnCenterStringFontHeightF   = 0.030
  
  resn@gsnCenterString  = ""   
  resn@gsnLeftString  = ""  
  resn@gsnRightString  = ""  
   
  resn@tiYAxisFontThicknessF    = 1.
  resn@tmXBMajorThicknessF      = 1.
  resn@tmXBMinorThicknessF      = 1.
  resn@tmYLMajorThicknessF      = 1.
  resn@tmYLMinorThicknessF      = 1.
  resn@tmYLLabelFontThicknessF  = 1.
  resn@tmBorderThicknessF       = 1. 
;-----------------------------------------------------------------------------
  resm                      = True         
  resm@gsnDraw              = False        ; don't draw yet
  resm@gsnFrame             = False        ; don't advance frame yet
  resm@gsnAddCyclic             = False             ; regional data 
   
  resm@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels
  resm@cnMinLevelValF       = -3.            ; set min contour level
  resm@cnMaxLevelValF       =  3.             ; set max contour level
  resm@cnLevelSpacingF      =  .3            ; set contour spacing   
  resm@cnFillOn             = True         ; turn on color fill
  resm@cnFillPalette        = "cmp_b2r"  
 
  resm@cnLinesOn            = False        ; True is default
  resm@cnLineLabelsOn       = False        ; True is default
  resm@cnInfoLabelOn        = False
  resm@lbLabelBarOn         = True       ; turn off individual lb's 
  ;resm@pmLabelBarParallelPosF = 0.50
  ;resm@pmLabelBarOrthogonalPosF = -0.11
  resm@lbOrientation        = "vertical"
  resm@lbTitleFontHeightF  = 0.012 
  resm@lbLabelFontHeightF  = 0.012
  resm@lbLeftMarginF       = 0.10
  resm@lbRightMarginF      = 0.55
  ;resm@lbLabelOffsetF	  = 0.20
  resm@lbLabelFont   	  = "Helvetica-Bold"     ; label font
  resm@lbTitleString       = "SSH response (~F10~e~S~~F21~-2~N~m)"
  ;resm@lbLabelPosition	  = "Top"
  resm@lbTitleAngleF       = 270.
  resm@lbTitlePosition      = "Right"
  resm@lbTitleDirection     = "Across"
   
  resm@lbLabelStride        = 10 
  
  resm@mpFillOn             = True        ; turn off map fill
  resm@mpMinLatF            = latS         ; zoom in on map
  resm@mpMaxLatF            = latN
  resm@mpMinLonF            = lonL
  resm@mpMaxLonF            = lonR
  resm@pmTickMarkDisplayMode = "Never"
  ;resm@tmXBLabelsOn         = False      ;no lon ticks
  ;resm@tmYLLabelsOn         = False      ;no lon ticks
  resm@tmYLMode   = "Explicit"
  resm@tmYLValues = ispan(-15,15,15)
  resm@tmYLLabels = (/"15~S~o~N~S", "0~S~o~N","15~S~o~N~N"/)
  resm@tmYLMinorValues  = ispan(-15,15,5) 
  resm@tmYLLabelFontHeightF = 0.017

  resm@tmXBMode   = "Explicit"
  resm@tmXBValues = ispan(-60,20,20)
  resm@tmXBLabels = (/"60~S~o~N~W","40~S~o~N~W","20~S~o~N~W","0~S~o~N~","20~S~o~N~E"/)
  ;resm@tmXBMinorValues  = ispan(-50,20,5)
  resm@tmXBLabelFontHeightF = 0.017 
 
  resm@mpGeophysicalLineColor = "grey" 
  resm@mpFillDrawOrder      = "postDraw"

  resm@vpHeightF 	          = 0.80        ; Changes the aspect ratio
  resm@vpWidthF 		  = 0.80
  resm@gsnCenterStringFontHeightF   = 0.030
  
  resm@gsnCenterString  = ""   
  resm@gsnLeftString  = ""  
  resm@gsnRightString  = ""  
    
  resm@vcRefAnnoOrthogonalPosF = -1.0             ; move ref vector up
  resm@vcRefMagnitudeF         = 0.01             ; define vector ref mag
  resm@vcRefAnnoString1          = "0.01 N m~S~-2~N~ "
  resm@vcRefLengthF            = 0.045            ; define length of vec ref
 ; vcres@vcGlyphStyle            = "CurlyVector"    ; turn on curly vectors
  resm@vcMinDistanceF          = 0.045
  resm@vcLineArrowThicknessF   = 1.
  resm@vcLineArrowColor        = "blue"           ; change vector color
  resm@gsnScalarContour	      = True
  resm@vcLineArrowHeadMaxSizeF = 0.05
  resm@vcLineArrowHeadMinSizeF = 0.01
 ; resm@vcMaxLevelCount         = 100
  
  resm@tiYAxisFontThicknessF    = 1.
  resm@tmXBMajorThicknessF      = 1.
  resm@tmXBMinorThicknessF      = 1.
  resm@tmYLMajorThicknessF      = 1.
  resm@tmYLMinorThicknessF      = 1.
  resm@tmYLLabelFontThicknessF  = 1.
  resm@tmBorderThicknessF       = 1. 
;-----------------------------------------------------------------------------
  ress2                      = True
  ress2@gsnDraw              = False        ; don't draw yet
  ress2@gsnFrame             = False        ; don't advance frame yet
  ress2@gsnAddCyclic             = False             ; regional data 
  ress2@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  ress2@cnMinLevelValF       = -20.              ; set min contour level
  ress2@cnMaxLevelValF       =  60.              ; set max contour level
  ress2@cnLevelSpacingF      =  10.              ; set contour spacing
  ress2@cnLineLabelsOn       = False
  ;ress2@cnLinesOn            = False       ; True is default
  ress2@cnInfoLabelOn        = False
  ress2@cnLineColor	    = "white"
  ress2@gsnContourZeroLineThicknessF = 2. 	; doubles thickness of zero contour
  ress2@gsnContourNegLineDashPattern = 1 	; sets negative contours to dash pattern 1
  ress2@gsnAddCyclic = False
 
;-----------------------------------------------------------------------------
  res3                      =   True
  res3@gsnDraw              = False        ; don't draw yet
  res3@gsnFrame             = False        ; don't advance frame yet
  ;res3@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
  ;res3@cnLevels		   = (/10./)   ; set levels
  res3@cnMinLevelValF       = -5.            ; set min contour level
  res3@cnMaxLevelValF       =  5.             ; set max contour level
  res3@cnLevelSpacingF      =  0.5            ; set contour spacing   
  res3@cnMonoLineThickness  = True
  res3@cnLineThicknessF     = 1.5
  res3@cnMonoLineColor	    = True
  ;res3@cnLineColor	    = "purple"
  ;res3@cnLineColor	    = "dodgerblue2"
  res3@cnLineColor	    = "black"
 ; res3@cnLineLabelFontColor = res2@cnLineColor
  res3@gsnAddCyclic         = False   
  res3@cnLinesOn            = True       ; True is default
  res3@cnLineLabelsOn       = True       ; True is default
  res3@cnInfoLabelOn        = False
  res3@gsnRightString   = ""
  res3@gsnLeftString   = ""
 ;res3@cnFillScaleF        = 2.         ; add extra density
  res3@cnLineLabelDensityF =  1.0
  res3@cnLabelMasking	   = True
  res3@gsnContourZeroLineThicknessF = 0.
 ;res3@gsnContourLineThicknessesScale = 1.
  res3@gsnContourNegLineDashPattern = 16
  
 ; res3@mpFillDrawOrder      = "preDraw"

;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
  ress                          = True
  ress@gsnDraw                  = False        ; don't draw yet
  ress@gsnFrame                 = False        ; don't advance frame yet
  ress@gsnAddCyclic             = False             ; regional data 
  ress@mpFillOn                 = True        ; turn off map fill
  ress@cnFillOn                 = True            ; turn on color
  ress@cnFillPalette            = "cmp_b2r"  ; set color map
  
  
  ress@cnMonoFillColor          = False
  ress@lbBoxLinesOn             = True
  ress@lbBoxSeparatorLinesOn    = True
  ress@cnFillDrawOrder          = "PreDraw" ; draw contours first
  
  ress@cnLevelSelectionMode     = "ManualLevels"  ; set manual contour levels
  ress@cnMinLevelValF           = -2.0            ; set min contour level
  ress@cnMaxLevelValF           =  2.0             ; set max contour level
  ress@cnLevelSpacingF          =  .2            ; set contour spacing   
 ; ress@cnLevels                =  0.2
  
  ress@mpGeophysicalLineColor   = "transparent"
  ress@pmTickMarkDisplayMode    = "Never"
   
  ress@cnLinesOn                = False        ; True is default
  ress@cnLineLabelsOn           = False        ; True is default
  ress@lbLabelBarOn             = True     ; turn off individual lb's
  
  ress@lbOrientation            = "vertical"
  ress@lbTitleFontHeightF       = 0.012 
  ress@lbLabelFontHeightF       = 0.012
  ress@lbLeftMarginF            = 0.0
  ress@lbRightMarginF           = 0.60
  ress@lbLabelOffsetF	       = 0.20
  ress@lbLabelFont   	       = "Helvetica-Bold"     ; label font
  ress@lbTitleString           = "Precipitation Response (mmday~S~-1~N~)"
  ;ress@lbLabelPosition	       = "Top"
  ress@lbTitleAngleF            = 270.
  ress@lbTitlePosition          = "Right"
  ress@lbTitleDirection         = "Across"

  ress@lbLabelStride            = 5
    
  ress@mpFillOn             = True        ; turn off map fill
  ress@mpMinLatF            = latS         ; zoom in on map
  ress@mpMaxLatF            = latN
  ress@mpMinLonF            = lonL
  ress@mpMaxLonF            = lonR
  ress@pmTickMarkDisplayMode = "Never"
  
  ress@lbTopMarginF             = 0.20
 
  ress@tmYLMode   = "Explicit"
  ress@tmYLValues = ispan(-15,15,15)
  ress@tmYLLabels = (/"15~S~o~N~S", "0~S~o~N","15~S~o~N~N"/)
  ress@tmYLMinorValues  = ispan(-15,15,5) 
  ress@tmYLLabelFontHeightF = 0.017

  ress@tmXBMode   = "Explicit"
  ress@tmXBValues = ispan(-60,20,20)
  ress@tmXBLabels = (/"60~S~o~N~W","40~S~o~N~W","20~S~o~N~W","0~S~o~N~","20~S~o~N~E"/)
  ;ress@tmXBMinorValues  = ispan(-50,20,5)
  ress@tmXBLabelFontHeightF = 0.017
 
  ;ress@mpGeophysicalLineColor  = "foreground"  
  ress@mpGeophysicalLineColor = "grey" 
  ress@mpFillDrawOrder      = "postDraw"

  ress@gsnCenterStringFontHeightF   = 0.030
  

  ress@gsnCenterString  = ""   
 ; ress@gsnLeftString  = "U-stress indices vs SSHA (cm/~F10~e~S~~F21~2~N~Nm~S~-2~N~)"  
  ress@gsnRightString  = ""  
 
  
  ress@tmXTOn = True ; Turn off top tickmarks and labels
  ress@tmYROn = True ; Turn off left tickmarks and labels
  
  ress@vpHeightF 	  = 0.80        ; Changes the aspect ratio
  ress@vpWidthF 		  = 0.80
 
  ress@tiYAxisFontThicknessF    = 1.
  ress@tmXBMajorThicknessF      = 1.
  ress@tmXBMinorThicknessF      = 1.
  ress@tmYLMajorThicknessF      = 1.
  ress@tmYLMinorThicknessF      = 1.
  ress@tmYLLabelFontThicknessF  = 1.
  ress@tmBorderThicknessF       = 1.
  
;************************************************
;-----------------------------------------------------------------------------
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
  opt@gsnShadeLow = 17                     ; use fill pattern #4
  opt@gsnShadeHigh = 17                  ; use pattern #2   
  ;opt@gsnShadeMid = 2   
  
  delete(resn@cnFillPalette)
  resn@cnFillPalette        = "amwg_blueyellowred"                        ;different colour bar  
;-----------------------------------------------------------------------------
  resn@gsnLeftString  = "Coupling of convection and SST anomalies"  
  plotA10 = gsn_csm_contour_map_ce(wks, d_SST, resn)
  plotB10 = gsn_csm_contour(wks,d_SSTd, res1)
  plotB10 = gsn_contour_shade(plotB10,p95,d_SSTd@_FillValue,opt)
  plotC10 = gsn_csm_contour(wks,d_precip,res3)
  overlay (plotA10, plotB10)
  overlay (plotA10,plotC10) 
  plot(4) = plotA10
;============================================================================= 
  delete(res@cnFillPalette)
  res@cnFillPalette        = "amwg_blueyellowred"                        ;different colour bar  
  res@gsnLeftString  = "Thermocline feedback"   
  plotA0 = gsn_csm_vector_scalar_map_ce(wks,d_USTRd(month|5,lat|:,lon|:),d_VSTRd(month|5,lat|:,lon|:),d_SSH,res)  ; create plot
  plotB0 = gsn_csm_contour(wks,d_SSHd,res1)
  plotB0 = gsn_contour_shade(plotB0,p95,d_SSHd@_FillValue,opt)
 
  overlay(plotA0,plotB0)
  plot(5)  =   plotA0 
;=======================================================================================================================
  delete(res@cnFillPalette)
  res@cnFillPalette        = "amwg_blueyellowred"                        ;different colour bar  
  res@gsnLeftString  = "Thermocline feedback"   
  plotA0J = gsn_csm_vector_scalar_map_ce(wks,d_USTRd_A(month|3,lat|:,lon|:),d_VSTRd_A(month|3,lat|:,lon|:),d_SSH_A,res)  ; create plot
  plotB0J = gsn_csm_contour(wks,d_SSHd_A,res1)
  plotB0J = gsn_contour_shade(plotB0J,p95,d_SSHd_A@_FillValue,opt)

  overlay(plotA0J,plotB0J)
  plot(1)  =   plotA0J   
;=====================================================================d_SSTd_A = d_SST_A
; create plots
  delete(ress@cnFillPalette)
  ress@cnFillPalette        = "amwg_blueyellowred"                        ;different colour bar  
  resn@gsnLeftString  = "Coupling of convection and SST anomalies"  
  plotA1 = gsn_csm_contour_map_ce(wks,d_SST_A, resn)
  plotB1 = gsn_csm_contour(wks,d_SSTd_A, res1)
  plotB1 = gsn_contour_shade(plotB1,p95,d_SSTd_A@_FillValue,opt)
  plotB11 = gsn_csm_contour(wks,d_precip_A,res3)
  overlay (plotA1, plotB1)
  overlay (plotA1, plotB11)
  plot(0) = plotA1
  
  resn@gsnLeftString  = "Coupling of convection and SST anomalies"  
  plotA1J = gsn_csm_contour_map_ce(wks,d_SST_M, resn)
  plotB1J = gsn_csm_contour(wks,d_SSTd_M, res1)
  plotB1J = gsn_contour_shade(plotB1J,p95,d_SSTd_M@_FillValue,opt)
  plotB11J = gsn_csm_contour(wks,d_precip_M,res3)
  overlay (plotA1J, plotB1J)
  overlay (plotA1J, plotB11J)
  plot(2) = plotA1J
;=====================================================================
  delete(ress@cnFillPalette)
  ress@cnFillPalette        = "amwg_blueyellowred"                        ;different colour bar  
 ; resn@gsnLeftString  = "Coupling of convection and SST anomalies"  
;  plotA1J = gsn_csm_contour_map_ce(wks,d_SST_JJA, resn)
;  plotB1J = gsn_csm_contour(wks,d_SSTd_JJA, res1)
;  plotB1J = gsn_contour_shade(plotB1J,p95,d_SSTd_JJA@_FillValue,opt)
;  plotB11J = gsn_csm_contour(wks,d_precip_JJA,res3)
 
;  overlay (plotA1J, plotB1J)
 ; overlay (plotA1J, plotB11J)
 ; plot(4) = plotA1J  
;==========================================================================
  delete(resm@cnFillPalette)
  resm@cnFillPalette        = "amwg_blueyellowred"                        ;different colour bar  
  
  resm@gsnLeftString  = "Thermocline feedback"   
  plotG10 = gsn_csm_vector_scalar_map_ce(wks,d_USTRd_M(month|4,lat|:,lon|:),d_VSTRd_M(month|4,lat|:,lon|:),d_SSH_M,resm)  ; create plot
  plotE10 = gsn_csm_contour(wks,d_SSHd_A,res1)
  plotE10 = gsn_contour_shade(plotE10,p95,d_SSHd_A@_FillValue,opt)
 
  overlay(plotG10,plotE10)
  plot(3)  =   plotG10 
  
 ; resm@gsnLeftString  = "Thermocline feedback"   
 ; plotG10 = gsn_csm_vector_scalar_map_ce(wks,d_USTRd_MJJA(lat|:,lon|:),d_VSTRd_MJJA(lat|:,lon|:),d_SSH_MJJA,resm)  ; create plot
 ; plotE10 = gsn_csm_contour(wks,d_SSHd_MJJA,res1)
 ; plotE10 = gsn_contour_shade(plotE10,p95,d_SSHd_MJJA@_FillValue,opt)
 
 ; overlay(plotG10,plotE10)
 ; plot(1)  =   plotG10 
;=====================================================================  e 3a shows the ∇ASSTA
  ; Draw polylines
  lon_values1 = (/-20., 0., 0., -20.,-20./)
  lat_values1 = (/ -3., -3., 3., 3., -3./)
  
; add the box
;************************************************
  resp                  = True                      ; polyline mods desired
  resp@gsLineColor      = "grey"                     ; color of lines
  resp@gsLineThicknessF = 2.0                       ; thickness of lines
  
; create array of dummy graphic variables. This is required, b/c each line
; must be associated with a unique dummy variable.  

  dum1 = new(4,graphic) 
  dum2 = new(4,graphic)  
  dum3 = new(4,graphic)  
  dum4 = new(4,graphic)  
  dum5 = new(4,graphic)  
  dum6 = new(4,graphic)  

; draw each line separately. Each line must contain two points.
  do i = 0 , 3
    dum1(i)=gsn_add_polyline(wks,plotA0,lon_values1(i:i+1),lat_values1(i:i+1),resp) 
    dum2(i)=gsn_add_polyline(wks,plotA1,lon_values1(i:i+1),lat_values1(i:i+1),resp) 
    dum3(i)=gsn_add_polyline(wks,plotA10,lon_values1(i:i+1),lat_values1(i:i+1),resp)        
    dum4(i)=gsn_add_polyline(wks,plotG10,lon_values1(i:i+1),lat_values1(i:i+1),resp)           
    dum5(i)=gsn_add_polyline(wks,plotA1J,lon_values1(i:i+1),lat_values1(i:i+1),resp)        
    dum6(i)=gsn_add_polyline(wks,plotA0J,lon_values1(i:i+1),lat_values1(i:i+1),resp)           
  end do  
;----------------------------------------------------------------------------- 
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

 tt = (/"a","b","c","d","e","f","g","h"/)
create_tt = gsn_create_text(wks,tt,txres5);
; tt_p0 = new(nplots,graphic)                ; create graphic array ; only needed if paneling
tt_p0 = new(6,graphic)                ; create graphic array ; only needed if paneling
;---Attach the polyline
     do i = 0,5
	print(i)
        tt_p0(i) = gsn_add_annotation(plot(i),create_tt(i),amres5)
     end do
delete(i)
;=============================================================================   
;===============================================================
  resP                                  = True         ; modify the panel plot
  resP@gsnDraw                          = True      ; don't draw yet
  resP@gsnFrame                         = False        ; don't advance frame yet
  resP@gsnPanelLabelBar    	        = False        ; add common colorbar
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
   resP@gsnPanelDebug      		  = True
   resP@amJust   			  = "BottomLeft"
   resP@gsnPanelFigureStringsPerimOn     = False
   resP@gsnPanelFigureStringsFontHeightF = 0.010
   resP@gsnPanelFigureStringsBackgroundFillColor = "transparent"
   resP@gsnPanelFigureStrings= (/"April","April","May","May","June","June","g","h","i","j","k","l"/) ; add strings to panel
   resP@gsnMaximize = True                        ; blow up plot
   
  gsn_panel(wks,plot(0:),(/3,2/),resP)               ; create panel plots
;*******************************************
;===============================================================  

  frame(wks)
end  

