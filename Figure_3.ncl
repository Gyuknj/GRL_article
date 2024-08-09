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
  
  yyyymm = yyyymm_time(yrStrt,yrLast,"integer")
 ; print(yyyymm) 
 
  xlabel = ispan(0,nmos-1,1)
  ;print(xlabel)
  
  p90 = 0.05

  mxlag = 7
  lag = 1	;first-order autocorrelation
  
  xy_colors_dheating        = (/"red","deepskyblue","darksalmon","black","deeppink4","yellow2","blue","seagreen3"/)
  xy_colors_thermocline        = (/"deeppink4","yellow2","deepskyblue","red","grey75","springgreen4"/)
;======================================================================================================================================
  f      = addfile("/home/gyuk/Downloads/GYUK/GEOMAR/Final_work/HadISST_sst.nc", "r")
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

  hadi = where(hadi.gt.-30.,hadi,hadi@_FillValue)	;removes padded values
  printVarSummary(hadi)       ; variable overview
  ;printMinMax(hadi,True)
;===================================================================================================================================
  hadiClm      = clmMonTLL(hadi)
  printVarSummary(hadiClm)                   ; (12,nlat,nlon)   

  hadi         = calcMonAnomTLL(hadi,hadiClm)    ; replace with anonamlies
  hadi@long_name = "ANOMALIES: "+hadi@long_name
  delete(hadiClm)
  ;printMinMax(hadi,True)
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
  delete(hadi_s)
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
;------------------------------------------------------------------------------------------------------------------------------------
  valsN = ind(jja_ahadi.le.-thres)
  print(valsN)
  print(valsN+yrStrt)             ;the years
;***************************************************************************************************************
;======================================================================================================================================
  f      = addfile("/home/gyuk/Downloads/q_net/q_net_regrid_1979_2022.nc", "r")
  TIME   = f->time_counter
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  Qnet   = f->sohefldo(iYYYY,:,:)
  Qnet!0 = "time"
  Qnet!1 = "lat"
  Qnet!2 = "lon"
  
  Qnet = lonFlip(Qnet)
  printVarSummary(Qnet)       ; variable overview
  printMinMax(Qnet, True)

  delete(f)		;	Not needed anymore
  delete(TIME)
  delete(YYYY)
  delete(iYYYY)
;===================================================================================================================================
  nyears = dimsizes(Qnet&time)/nmos
  ;print(nyears)
;===================================================================================================================================  
  rho        = 1000.  
  cw         = 4000.
  
  constant   = rho * cw
;======================================================================================================================================
;-------------------------------------------------------------------------------  
  fd1   = addfile("/home/gyuk/Downloads/mld_01/mld_regrid_01.nc", "r")   
  TIME   = fd1->time_counter
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  mld   = fd1->somxl010(iYYYY,:,:)
  mld = lonFlip(mld)
  
  mld!0 = "time"
  mld!1 = "lat"
  mld!2 = "lon"
  
  delete(fd1)             ;       Not needed anymore
  delete(TIME)
  delete(YYYY)
  delete(iYYYY)
  
  printVarSummary(mld)       ; variable overview
  ;printMinMax(mld,True)
  
;===================================================================================================================================
;====================================================================================================================================
   heatflux = [(Qnet)/(constant*mld)]
   copy_VarCoords(mld(:,:,:),heatflux(:,:,:))                             ; copy coordinate arrays
   printVarSummary(heatflux)
   printMinMax(heatflux, 0)
;======================================================================================================================================
;======================================================================================================================================
;-------------------------------------------------------------------------------  
   aheatflux_s = wgt_areaave_Wrap(heatflux(time|:,{lon|-20:0}, {lat|-3:3}),1.0, 1.0, 0)        ;atl3 index
   printVarSummary(aheatflux_s)
   printMinMax(aheatflux_s,True)
;-------------------------------------------------------------------------------  
   aheatflux_s = aheatflux_s*1000000.
   printVarSummary(aheatflux_s)
   printMinMax(aheatflux_s,True)

  oras5_sq     = onedtond(ndtooned(aheatflux_s),(/nyears,nmos/))
  oras5_sq!0 = "time"
  oras5_sq!1 = "month"
  printVarSummary(oras5_sq)
  
  delete(aheatflux_s)
  delete(Qnet)
  
  oras5_sq = dtrend_n(oras5_sq,False,0)
  printVarSummary(oras5_sq)
;======================================================================================================================================
;-------------------------------------------------------------------------------  
  xxqP = oras5_sq(valsP,:)
  printVarSummary(xxqP)
  
  xxqN = oras5_sq(valsN,:)
  printVarSummary(xxqN)
  
  delete(oras5_sq)
;======================================================================================================================================
  xxqP1 = (xxqP(time|2,month|:))
  xxqP2 = (xxqP(time|3,month|:))
  xxqP3 = (xxqP(time|4,month|:))
  xxqP4 = (xxqP(time|6,month|:))
  xxqP5 = (xxqP(time|7,month|:))
  
  canonical_qP = (/xxqP1,xxqP2,xxqP3,xxqP4,xxqP5/)       
  copy_VarCoords_n(xxqP,canonical_qP,(/0,1/))
  canonical_qP!0 = "nyears"    
  printVarSummary(canonical_qP)
;-------------------------------------------------------------------------------  
  xxqN1 = (xxqN(time|0,month|:))
  xxqN2 = (xxqN(time|1,month|:))
  xxqN3 = (xxqN(time|2,month|:))
  xxqN4 = (xxqN(time|5,month|:))
  
  canonical_qN = (/xxqN1,xxqN2,xxqN3,xxqN4/)       
  copy_VarCoords_n(xxqN,canonical_qN,(/0,1/))
  canonical_qN!0 = "nyears"    
  printVarSummary(canonical_qN) 
  
  delete([/xxqP1, xxqP2, xxqP3, xxqP4, xxqP5, xxqN1, xxqN2, xxqN3, xxqN4/])
;===================================== 
  canonical_qP_s = dim_avg_n_Wrap(canonical_qP,0)
  printVarSummary(canonical_qP_s) 
;***************************************************************************************************************
;     Difference in Canonical events
;*************************************************************************************************************** 
  canonical_qP_difference = (dim_avg_n_Wrap(canonical_qP,0) - dim_avg_n_Wrap(canonical_qN,0))/(2.)
  copy_VarCoords(canonical_qP_s,canonical_qP_difference)                             ; copy coordinate arrays 
  printVarSummary(canonical_qP_difference)
;*************************************************************************************************************** 
  pc1q_4_sd = canonical_qP_difference  ; Transfer metdata for plotting sig. values        

  pc1q_4_sd = ttest((dim_avg_n_Wrap(canonical_qP,0)),(dim_stddev_n_Wrap(canonical_qP,0))^2,5\
               ,(dim_avg_n_Wrap(canonical_qN,0)),(dim_stddev_n_Wrap(canonical_qN,0))^2,4,True, False)
  pc1q_4_sd@long_name = "Probability: difference between means"
  printVarSummary(pc1q_4_sd)

  pc1q_4_sd = where(abs(pc1q_4_sd).le.p90,canonical_qP_difference,pc1q_4_sd@_FillValue) ; Sig.u winds, where UWND or VWND is significanlt
  printVarSummary(pc1q_4_sd)
  
  delete([/canonical_qP, canonical_qN, canonical_qP_s/])
;======================================================================================================================================
;-------------------------------------------------------------------------------  
;======================================================================================================================================
  xxqP6 = (xxqP(time|0,month|:))
  xxqP7 = (xxqP(time|1,month|:))
  xxqP8 = (xxqP(time|5,month|:))
  xxqP9 = (xxqP(time|8,month|:))
  
  canonical_qPN = (/xxqP6,xxqP7,xxqP8,xxqP9/)       
  copy_VarCoords_n(xxqP,canonical_qPN,(/0,1/))
  canonical_qPN!0 = "nyears"    
  printVarSummary(canonical_qPN)
;-------------------------------------------------------------------------------  
  xxqN5 = (xxqN(time|3,month|:))
  xxqN6 = (xxqN(time|4,month|:))
  xxqN7 = (xxqN(time|6,month|:))
  xxqN8 = (xxqN(time|7,month|:))

  canonical_qNN = (/xxqN5,xxqN6,xxqN7,xxqN8/)       
  copy_VarCoords_n(xxqN,canonical_qNN,(/0,1/))
  canonical_qNN!0 = "nyears"    
  printVarSummary(canonical_qNN) 
  
  delete([/xxqP6, xxqP7, xxqP8, xxqP9, xxqN5, xxqN6, xxqN7, xxqN8/])
;===================================== 
  canonical_qPN_s = dim_avg_n_Wrap(canonical_qPN,0)
  printVarSummary(canonical_qPN_s) 
;***************************************************************************************************************
;     Difference in Canonical events
;*************************************************************************************************************** 
  canonical_qPN_difference = (dim_avg_n_Wrap(canonical_qPN,0) - dim_avg_n_Wrap(canonical_qNN,0))/(2.)
  copy_VarCoords(canonical_qPN_s,canonical_qPN_difference)                             ; copy coordinate arrays 
  printVarSummary(canonical_qPN_difference)
;*************************************************************************************************************** 
  pc1qN_4_sd = canonical_qPN_difference  ; Transfer metdata for plotting sig. values        

  pc1qN_4_sd = ttest((dim_avg_n_Wrap(canonical_qPN,0)),(dim_stddev_n_Wrap(canonical_qPN,0))^2,5\
               ,(dim_avg_n_Wrap(canonical_qNN,0)),(dim_stddev_n_Wrap(canonical_qNN,0))^2,4,True, False)
  pc1qN_4_sd@long_name = "Probability: difference between means"
  printVarSummary(pc1qN_4_sd)

  pc1qN_4_sd = where(abs(pc1qN_4_sd).le.p90,canonical_qPN_difference,pc1qN_4_sd@_FillValue) ; Sig.u winds, where UWND or VWND is significanlt
  printVarSummary(pc1qN_4_sd)
  
  delete([/canonical_qPN, canonical_qNN, canonical_qPN_s/])
  delete(xxqP)
  delete(xxqN)
;*************************************************************************************************************** 
;======================================================================================================================================
;-------------------------------------------------------------------------------  
  fwnd1   = addfile("/home/gyuk/Downloads/New_heat_budget/w_mlbase_oras5_1979_2023.nc", "r")   
  TIME   = fwnd1->TIME_COUNTER
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  Vel_base   = fwnd1->W_MLBASE(iYYYY,:,:)
  
  Vel_base = lonFlip(Vel_base)
  
  Vel_base!0 = "time"
  Vel_base!1 = "lat"
  Vel_base!2 = "lon"
  
  delete(fwnd1)             ;       Not needed anymore
  delete(TIME)
  delete(YYYY)
  delete(iYYYY)
  
  printVarSummary(Vel_base)       ; variable overview
  printMinMax(Vel_base,True)
;======================================================================================================================================
  f      = addfile("/home/gyuk/Downloads/New_heat_budget/TEMP_K_MLD_global_oras5.nc", "r")
  TIME   = f->TIME_COUNTER
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  temp_mld   = f->TEMP_K_MLD(iYYYY,:,:)
  
  temp_mld!0 = "time"
  temp_mld!1 = "lat"
  temp_mld!2 = "lon"
  
  temp_mld = lonFlip(temp_mld)
  printVarSummary(temp_mld)       ; variable overview
  printMinMax(temp_mld, True)

  delete(f)		;	Not needed anymore
  delete(TIME)
  delete(YYYY)
  delete(iYYYY)
;===================================================================================================================================
  nyears = dimsizes(temp_mld&time)/nmos
  ;print(nyears)
;===================================================================================================================================
;====================================================================================================================================
  f      = addfile("/home/gyuk/Downloads/New_heat_budget/temp_mlbase_global_oras5.nc", "r")
  TIME   = f->TIME_COUNTER
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  temp_mlbase   = f->TEMP_MLBASE(iYYYY,:,:)
  
  temp_mlbase = lonFlip(temp_mlbase)
  printVarSummary(temp_mlbase)       ; variable overview
  delete(f)             ;       Not needed anymore
  delete(TIME)
  delete(YYYY)
  delete(iYYYY)
  
  temp_mlbase!0 = "time"
  temp_mlbase!1 = "lat"
  temp_mlbase!2 = "lon"
;======================================================================================================================================
;====================================================================================================================================== 
   entrainment_f   = (1./mld)*[temp_mld - temp_mlbase]*(-Vel_base)
   copy_VarCoords(Vel_base(:,:,:), entrainment_f)
   printVarSummary(entrainment_f)
   printMinMax(entrainment_f, 0)
   
   delete(mld)
   delete(temp_mlbase) 
   delete(Vel_base)
;======================================================================================================================================
;======================================================================================================================================
;-------------------------------------------------------------------------------  
   aentrainment_s = wgt_areaave_Wrap(entrainment_f(time|:,{lon|-20:0}, {lat|-3:3}),1.0, 1.0, 0)        ;atl3 index
   printVarSummary(aentrainment_s)
   printMinMax(aentrainment_s,True)
;-------------------------------------------------------------------------------  
   aentrainment_s = aentrainment_s*1000000.
   printVarSummary(aentrainment_s)
   printMinMax(aentrainment_s,True)

  oras5_se     = onedtond(ndtooned(aentrainment_s),(/nyears,nmos/))
  oras5_se!0 = "time"
  oras5_se!1 = "month"
  printVarSummary(oras5_se)
  
  delete(aentrainment_s)
  
  oras5_se = dtrend_n(oras5_se,False,0)
  printVarSummary(oras5_se)
;======================================================================================================================================
;-------------------------------------------------------------------------------  
  xxeP = oras5_se(valsP,:)
  printVarSummary(xxeP)
  
  xxeN = oras5_se(valsN,:)
  printVarSummary(xxeN)
  
  delete(oras5_se)
;======================================================================================================================================
  xxeP1 = (xxeP(time|2,month|:))
  xxeP2 = (xxeP(time|3,month|:))
  xxeP3 = (xxeP(time|4,month|:))
  xxeP4 = (xxeP(time|6,month|:))
  xxeP5 = (xxeP(time|7,month|:))
  
  canonical_eP = (/xxeP1,xxeP2,xxeP3,xxeP4,xxeP5/)       
  copy_VarCoords_n(xxeP,canonical_eP,(/0,1/))
  canonical_eP!0 = "nyears"    
  printVarSummary(canonical_eP)
;-------------------------------------------------------------------------------  
  xxeN1 = (xxeN(time|0,month|:))
  xxeN2 = (xxeN(time|1,month|:))
  xxeN3 = (xxeN(time|2,month|:))
  xxeN4 = (xxeN(time|5,month|:))
  
  canonical_eN = (/xxeN1,xxeN2,xxeN3,xxeN4/)       
  copy_VarCoords_n(xxeN,canonical_eN,(/0,1/))
  canonical_eN!0 = "nyears"    
  printVarSummary(canonical_eN) 
  
  delete([/xxeP1, xxeP2, xxeP3, xxeP4, xxeP5, xxeN1, xxeN2, xxeN3, xxeN4/])
;===================================== 
  canonical_eP_s = dim_avg_n_Wrap(canonical_eP,0)
  printVarSummary(canonical_eP_s) 
;***************************************************************************************************************
;     Difference in Canonical events
;*************************************************************************************************************** 
  canonical_eP_difference = (dim_avg_n_Wrap(canonical_eP,0) - dim_avg_n_Wrap(canonical_eN,0))/(2.)
  copy_VarCoords(canonical_eP_s,canonical_eP_difference)                             ; copy coordinate arrays 
  printVarSummary(canonical_eP_difference)
;*************************************************************************************************************** 
  pc1e_4_sd = canonical_eP_difference  ; Transfer metdata for plotting sig. values        

  pc1e_4_sd = ttest((dim_avg_n_Wrap(canonical_eP,0)),(dim_stddev_n_Wrap(canonical_eP,0))^2,5\
               ,(dim_avg_n_Wrap(canonical_eN,0)),(dim_stddev_n_Wrap(canonical_eN,0))^2,4,True, False)
  pc1e_4_sd@long_name = "Probability: difference between means"
  printVarSummary(pc1e_4_sd)

  pc1e_4_sd = where(abs(pc1e_4_sd).le.p90,canonical_eP_difference,pc1e_4_sd@_FillValue) ; Sig.u winds, where UWND or VWND is significanlt
  printVarSummary(pc1e_4_sd)
  
  delete([/canonical_eP, canonical_eN, canonical_eP_s/])
;======================================================================================================================================
;======================================================================================================================================
;======================================================================================================================================
  xxeP6 = (xxeP(time|0,month|:))
  xxeP7 = (xxeP(time|4,month|:))
  xxeP8 = (xxeP(time|5,month|:))
  xxeP9 = (xxeP(time|8,month|:))
 
  canonical_ePN = (/xxeP6,xxeP7,xxeP8,xxeP9/)       
  copy_VarCoords_n(xxeP,canonical_ePN,(/0,1/))
  canonical_ePN!0 = "nyears"    
  printVarSummary(canonical_ePN)
;-------------------------------------------------------------------------------  
  xxeN5 = (xxeN(time|3,month|:))
  xxeN6 = (xxeN(time|4,month|:))
  xxeN7 = (xxeN(time|6,month|:))
  xxeN8 = (xxeN(time|7,month|:))
  
  canonical_eNN = (/xxeN5,xxeN6,xxeN7,xxeN8/)       
  copy_VarCoords_n(xxeN,canonical_eNN,(/0,1/))
  canonical_eNN!0 = "nyears"    
  printVarSummary(canonical_eNN) 
  
  delete([/xxeP6, xxeP7, xxeP8, xxeP9, xxeN5, xxeN6, xxeN7, xxeN8/])
;===================================== 
  canonical_ePN_s = dim_avg_n_Wrap(canonical_ePN,0)
  printVarSummary(canonical_ePN_s) 
;***************************************************************************************************************
;     Difference in Canonical events
;*************************************************************************************************************** 
  canonical_ePN_difference = (dim_avg_n_Wrap(canonical_ePN,0) - dim_avg_n_Wrap(canonical_eNN,0))/(2.)
  copy_VarCoords(canonical_ePN_s,canonical_ePN_difference)                             ; copy coordinate arrays 
  printVarSummary(canonical_ePN_difference)
;*************************************************************************************************************** 
  pc1eN_4_sd = canonical_ePN_difference  ; Transfer metdata for plotting sig. values        

  pc1eN_4_sd = ttest((dim_avg_n_Wrap(canonical_ePN,0)),(dim_stddev_n_Wrap(canonical_ePN,0))^2,5\
               ,(dim_avg_n_Wrap(canonical_eNN,0)),(dim_stddev_n_Wrap(canonical_eNN,0))^2,4,True, False)
  pc1eN_4_sd@long_name = "Probability: difference between means"
  printVarSummary(pc1eN_4_sd)

  pc1eN_4_sd = where(abs(pc1eN_4_sd).le.p90,canonical_ePN_difference,pc1eN_4_sd@_FillValue) ; Sig.u winds, where UWND or VWND is significanlt
  printVarSummary(pc1eN_4_sd)
  
  delete([/canonical_ePN, canonical_eNN, canonical_ePN_s/])
;======================================================================================================================================
;======================================================================================================================================
  f      = addfile("/home/gyuk/Downloads/New_heat_budget/U_MLD_global_oras5.nc", "r")
  TIME   = f->TIME_COUNTER
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  u_mld   = f->U_MLD(iYYYY,:,:)
  u_mld = lonFlip(u_mld)
  printVarSummary(u_mld)       ; variable overview
  delete(f)             ;       Not needed anymore
  delete(TIME)
  delete(YYYY)
  delete(iYYYY)
  
  u_mld!0 = "time"
  u_mld!1 = "lat"
  u_mld!2 = "lon"
;======================================================================================================================================
  fwnd1   = addfile("/home/gyuk/Downloads/New_heat_budget/V_MLD_global_oras5.nc", "r")   
  TIME   = fwnd1->TIME_COUNTER
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  v_mld   = fwnd1->V_MLD(iYYYY,:,:)
  v_mld = lonFlip(v_mld)
  
  v_mld!0 = "time"
  v_mld!1 = "lat"
  v_mld!2 = "lon"
  
  delete(fwnd1)             ;       Not needed anymore
  delete(TIME)
  delete(YYYY)
  delete(iYYYY)
  
  printVarSummary(v_mld)       ; variable overview
  ;printMinMax(v_mld,True)
;==============================================================
   cyclic    = True   ; Grid is global
   opt_adv   = 0     ; Return only the advected variable

   long_name = "temp advection: Global"
   units     = "K/s"

   advection   = advect_variable_cfd(u_mld,v_mld,temp_mld, temp_mld&lat,temp_mld&lon, cyclic, long_name, units, opt_adv)
  ; copy_VarCoords(Vel_base(:,:,:), advection)
   printVarSummary(advection)
   printMinMax(advection, 0)
   
   delete(u_mld)
   delete(v_mld)
;======================================================================================================================================
;======================================================================================================================================
   aadvection_s = wgt_areaave_Wrap(advection(time|:,{lon|-20:0}, {lat|-3:3}),1.0, 1.0, 0)        ;atl3 index
   printVarSummary(aadvection_s)
   printMinMax(aadvection_s,True)
;-------------------------------------------------------------------------------  
   aadvection_s = aadvection_s*1000000.
   printVarSummary(aadvection_s)
   printMinMax(aadvection_s,True)

  oras5_sa     = onedtond(ndtooned(aadvection_s),(/nyears,nmos/))
  oras5_sa!0 = "time"
  oras5_sa!1 = "month"
  printVarSummary(oras5_sa)
  
  delete(aadvection_s)
  
  oras5_sa = dtrend_n(oras5_sa,False,0)
  printVarSummary(oras5_sa)
;======================================================================================================================================
;-------------------------------------------------------------------------------  
  xxaP = oras5_sa(valsP,:)
  printVarSummary(xxaP)
  
  xxaN = oras5_sa(valsN,:)
  printVarSummary(xxaN)
  
  delete(oras5_sa)
;======================================================================================================================================
  xxaP1 = (xxaP(time|2,month|:))
  xxaP2 = (xxaP(time|3,month|:))
  xxaP3 = (xxaP(time|4,month|:))
  xxaP4 = (xxaP(time|6,month|:))
  xxaP5 = (xxaP(time|7,month|:))
  
  canonical_aP = (/xxaP1,xxaP2,xxaP3,xxaP4,xxaP5/)       
  copy_VarCoords_n(xxaP,canonical_aP,(/0,1/))
  canonical_aP!0 = "nyears"    
  printVarSummary(canonical_aP)
;-------------------------------------------------------------------------------  
  xxaN1 = (xxaN(time|0,month|:))
  xxaN2 = (xxaN(time|1,month|:))
  xxaN3 = (xxaN(time|2,month|:))
  xxaN4 = (xxaN(time|5,month|:))
  
  canonical_aN = (/xxaN1,xxaN2,xxaN3,xxaN4/)       
  copy_VarCoords_n(xxaN,canonical_aN,(/0,1/))
  canonical_aN!0 = "nyears"    
  printVarSummary(canonical_aN) 
  
  delete([/xxaP1, xxaP2, xxaP3, xxaP4, xxaP5, xxaN1, xxaN2, xxaN3, xxaN4/])
;===================================== 
  canonical_aP_s = dim_avg_n_Wrap(canonical_aP,0)
  printVarSummary(canonical_aP_s) 
;***************************************************************************************************************
;     Difference in Canonical events
;*************************************************************************************************************** 
  canonical_aP_difference = (dim_avg_n_Wrap(canonical_aP,0) - dim_avg_n_Wrap(canonical_aN,0))/(2.)
  copy_VarCoords(canonical_aP_s,canonical_aP_difference)                             ; copy coordinate arrays 
  printVarSummary(canonical_aP_difference)
;*************************************************************************************************************** 
  pc1a_4_sd = canonical_aP_difference  ; Transfer metdata for plotting sig. values        

  pc1a_4_sd = ttest((dim_avg_n_Wrap(canonical_aP,0)),(dim_stddev_n_Wrap(canonical_aP,0))^2,5\
               ,(dim_avg_n_Wrap(canonical_aN,0)),(dim_stddev_n_Wrap(canonical_aN,0))^2,4,True, False)
  pc1a_4_sd@long_name = "Probability: difference between means"
  printVarSummary(pc1a_4_sd)

  pc1a_4_sd = where(abs(pc1a_4_sd).le.p90,canonical_aP_difference,pc1a_4_sd@_FillValue) ; Sig.u winds, where UWND or VWND is significanlt
  printVarSummary(pc1a_4_sd)
  
  delete([/canonical_aP, canonical_aN, canonical_aP_s/])
;======================================================================================================================================
;======================================================================================================================================
;======================================================================================================================================
  xxaP6 = (xxaP(time|0,month|:))
  xxaP7 = (xxaP(time|1,month|:))
  xxaP8 = (xxaP(time|5,month|:))
  xxaP9 = (xxaP(time|8,month|:))
  
  canonical_aPN = (/xxaP6,xxaP7,xxaP8,xxaP9/)       
  copy_VarCoords_n(xxaP,canonical_aPN,(/0,1/))
  canonical_aPN!0 = "nyears"    
  printVarSummary(canonical_aPN)
;-------------------------------------------------------------------------------  
  xxaN5 = (xxaN(time|3,month|:))
  xxaN6 = (xxaN(time|4,month|:))
  xxaN7 = (xxaN(time|6,month|:))
  xxaN8 = (xxaN(time|7,month|:))
  
  canonical_aNN = (/xxaN5,xxaN6,xxaN7,xxaN8/)       
  copy_VarCoords_n(xxaN,canonical_aNN,(/0,1/))
  canonical_aNN!0 = "nyears"    
  printVarSummary(canonical_aNN) 
  
  delete([/xxaP6, xxaP7, xxaP8, xxaP9, xxaN5, xxaN6, xxaN7, xxaN8/])
;===================================== 
  canonical_aPN_s = dim_avg_n_Wrap(canonical_aPN,0)
  printVarSummary(canonical_aPN_s) 
;***************************************************************************************************************
;     Difference in Canonical events
;*************************************************************************************************************** 
  canonical_aPN_difference = (dim_avg_n_Wrap(canonical_aPN,0) - dim_avg_n_Wrap(canonical_aNN,0))/(2.)
  copy_VarCoords(canonical_aPN_s,canonical_aPN_difference)                             ; copy coordinate arrays 
  printVarSummary(canonical_aPN_difference)
;*************************************************************************************************************** 
  pc1aN_4_sd = canonical_aPN_difference  ; Transfer metdata for plotting sig. values        

  pc1aN_4_sd = ttest((dim_avg_n_Wrap(canonical_aPN,0)),(dim_stddev_n_Wrap(canonical_aPN,0))^2,5\
               ,(dim_avg_n_Wrap(canonical_aNN,0)),(dim_stddev_n_Wrap(canonical_aNN,0))^2,4,True, False)
  pc1aN_4_sd@long_name = "Probability: difference between means"
  printVarSummary(pc1aN_4_sd)

  pc1aN_4_sd = where(abs(pc1aN_4_sd).le.p90,canonical_aPN_difference,pc1aN_4_sd@_FillValue) ; Sig.u winds, where UWND or VWND is significanlt
  printVarSummary(pc1aN_4_sd)
  
  delete([/canonical_aPN, canonical_aNN, canonical_aPN_s/])
;======================================================================================================================================
;======================================================================================================================================
   dTdt = center_finite_diff_n((temp_mld),temp_mld&time,False,0,0)
   copy_VarCoords(temp_mld(:,:,:),dTdt(:,:,:))                             ; copy coordinate arrays
   printVarSummary(dTdt)
   printMinMax(dTdt, 0)  
;======================================================================================================================================
   adTdt_s = wgt_areaave_Wrap(dTdt(time|:,{lon|-20:0}, {lat|-3:3}),1.0, 1.0, 0)        ;atl3 index
   printVarSummary(adTdt_s)
   printMinMax(adTdt_s,True)
;======================================================================================================================================
;-------------------------------------------------------------------------------  
   adTdt_s = adTdt_s*1000000.
   printVarSummary(adTdt_s)
   printMinMax(adTdt_s,True)

  oras5_s     = onedtond(ndtooned(adTdt_s),(/nyears,nmos/))
  oras5_s!0 = "time"
  oras5_s!1 = "month"
  printVarSummary(oras5_s)
  
  delete(adTdt_s)
  
  oras5_s = dtrend_n(oras5_s,False,0)
  printVarSummary(oras5_s)
;======================================================================================================================================
;-------------------------------------------------------------------------------  
  xxdTP = oras5_s(valsP,:)
  printVarSummary(xxaP)
  
  xxdTN = oras5_s(valsN,:)
  printVarSummary(xxdTN)
  
  delete(oras5_s)
;======================================================================================================================================
  xxdTP1 = (xxdTP(time|2,month|:))
  xxdTP2 = (xxdTP(time|3,month|:))
  xxdTP3 = (xxdTP(time|4,month|:))
  xxdTP4 = (xxdTP(time|6,month|:))
  xxdTP5 = (xxdTP(time|7,month|:))
  
  canonical_dTP = (/xxdTP1,xxdTP2,xxdTP3,xxdTP4,xxdTP5/)       
  copy_VarCoords_n(xxdTP,canonical_dTP,(/0,1/))
  canonical_dTP!0 = "nyears"    
  printVarSummary(canonical_dTP)
;-------------------------------------------------------------------------------  
  xxdTN1 = (xxdTN(time|0,month|:))
  xxdTN2 = (xxdTN(time|1,month|:))
  xxdTN3 = (xxdTN(time|2,month|:))
  xxdTN4 = (xxdTN(time|5,month|:))
  
  canonical_dTN = (/xxdTN1,xxdTN2,xxdTN3,xxdTN4/)       
  copy_VarCoords_n(xxdTN,canonical_dTN,(/0,1/))
  canonical_dTN!0 = "nyears"    
  printVarSummary(canonical_dTN)
  
  delete([/xxdTP1, xxdTP2, xxdTP3, xxdTP4, xxdTP5, xxdTN1, xxdTN2, xxdTN3, xxdTN4/]) 
;===================================== 
  canonical_dTP_s = dim_avg_n_Wrap(canonical_dTP,0)
  printVarSummary(canonical_dTP_s) 
;***************************************************************************************************************
;     Difference in Canonical events
;*************************************************************************************************************** 
  canonical_dTP_difference = (dim_avg_n_Wrap(canonical_dTP,0) - dim_avg_n_Wrap(canonical_dTN,0))/(2.)
  copy_VarCoords(canonical_dTP_s,canonical_dTP_difference)                             ; copy coordinate arrays 
  printVarSummary(canonical_dTP_difference)
;*************************************************************************************************************** 
  pc1dT_4_sd = canonical_dTP_difference  ; Transfer metdata for plotting sig. values        

  pc1dT_4_sd = ttest((dim_avg_n_Wrap(canonical_dTP,0)),(dim_stddev_n_Wrap(canonical_dTP,0))^2,5\
               ,(dim_avg_n_Wrap(canonical_dTN,0)),(dim_stddev_n_Wrap(canonical_dTN,0))^2,4,True, False)
  pc1dT_4_sd@long_name = "Probability: difference between means"
  printVarSummary(pc1dT_4_sd)

  pc1dT_4_sd = where(abs(pc1dT_4_sd).le.p90,canonical_dTP_difference,pc1dT_4_sd@_FillValue) ; Sig.u winds, where UWND or VWND is significanlt
  printVarSummary(pc1dT_4_sd)
  
  delete([/canonical_dTP, canonical_dTN/]) 
;=======================================================================================================================================
;=======================================================================================================================================
;======================================================================================================================================
  xxdTP6 = (xxdTP(time|0,month|:))
  xxdTP7 = (xxdTP(time|1,month|:))
  xxdTP8 = (xxdTP(time|5,month|:))
  xxdTP9 = (xxdTP(time|8,month|:))
  
  canonical_dTPN = (/xxdTP6,xxdTP7,xxdTP8,xxdTP9/)       
  copy_VarCoords_n(xxdTP,canonical_dTPN,(/0,1/))
  canonical_dTPN!0 = "nyears"    
  printVarSummary(canonical_dTPN)
;-------------------------------------------------------------------------------  
  xxdTN5 = (xxdTN(time|3,month|:))
  xxdTN6 = (xxdTN(time|4,month|:))
  xxdTN7 = (xxdTN(time|6,month|:))
  xxdTN8 = (xxdTN(time|7,month|:))
  
  canonical_dTNN = (/xxdTN5,xxdTN6,xxdTN7,xxdTN8/)       
  copy_VarCoords_n(xxdTN,canonical_dTNN,(/0,1/))
  canonical_dTNN!0 = "nyears"    
  printVarSummary(canonical_dTNN) 
  
  delete([/xxdTP6, xxdTP7, xxdTP8, xxdTP9, xxdTN5, xxdTN6, xxdTN7, xxdTN8/]) 
;===================================== 
  canonical_dTPN_s = dim_avg_n_Wrap(canonical_dTPN,0)
  printVarSummary(canonical_dTPN_s) 
;***************************************************************************************************************
;     Difference in Canonical events
;*************************************************************************************************************** 
  canonical_dTPN_difference = (dim_avg_n_Wrap(canonical_dTPN,0) - dim_avg_n_Wrap(canonical_dTNN,0))/(2.)
  copy_VarCoords(canonical_dTPN_s,canonical_dTPN_difference)                             ; copy coordinate arrays 
  printVarSummary(canonical_dTPN_difference)
;*************************************************************************************************************** 
  pc1dTN_4_sd = canonical_dTPN_difference  ; Transfer metdata for plotting sig. values        

  pc1dTN_4_sd = ttest((dim_avg_n_Wrap(canonical_dTPN,0)),(dim_stddev_n_Wrap(canonical_dTPN,0))^2,5\
               ,(dim_avg_n_Wrap(canonical_dTNN,0)),(dim_stddev_n_Wrap(canonical_dTNN,0))^2,4,True, False)
  pc1dTN_4_sd@long_name = "Probability: difference between means"
  printVarSummary(pc1dTN_4_sd)

  pc1dTN_4_sd = where(abs(pc1dTN_4_sd).le.p90,canonical_dTPN_difference,pc1dTN_4_sd@_FillValue) ; Sig.u winds, where UWND or VWND is significanlt
  printVarSummary(pc1dTN_4_sd)
  
  delete([/canonical_dTPN, canonical_dTNN, canonical_dTPN_s/]) 
;=======================================================================================================================================
;=======================================================================================================================================

  x =   (canonical_dTP_difference) - (canonical_qP_difference + canonical_eP_difference + canonical_aP_difference)	;residual
 ; printVarSummary(x1)
 
  dT_dt_terms1 = (/canonical_qP_difference,canonical_eP_difference,canonical_aP_difference,canonical_dTP_difference,x/)
  copy_VarCoords_n(canonical_dTP_s,dT_dt_terms1,(/0,1/))
  dT_dt_terms1!0 = "lev"
  dT_dt_terms1!1 = "month"
  printVarSummary(dT_dt_terms1)
  
  x_sd =   (pc1dT_4_sd) - (pc1q_4_sd + pc1e_4_sd + pc1a_4_sd)		;residual
 ; printVarSummary(x_sd)
 
  dT_dt_terms_sd1 = (/pc1q_4_sd,pc1e_4_sd,pc1a_4_sd,pc1dT_4_sd,x_sd/)
  copy_VarCoords_n(canonical_dTP_s,dT_dt_terms_sd1,(/0,1/))
  dT_dt_terms_sd1!0 = "lev"
  dT_dt_terms_sd1!1 = "month"
  printVarSummary(dT_dt_terms_sd1)
;=======================================================================================================================================
;=======================================================================================================================================

  xN =   (canonical_dTPN_difference) - (canonical_qPN_difference + canonical_ePN_difference + canonical_aPN_difference)	;residual
  ;printVarSummary(x1N)
 
  dT_dt_terms1N = (/canonical_qPN_difference,canonical_ePN_difference,canonical_aPN_difference,canonical_dTPN_difference,xN/)
  copy_VarCoords_n(canonical_dTP_s,dT_dt_terms1N,(/0,1/))
  dT_dt_terms1N!0 = "lev"
  dT_dt_terms1N!1 = "month"
  printVarSummary(dT_dt_terms1N)
  
  x_sdN =   (pc1dTN_4_sd) - (pc1qN_4_sd + pc1eN_4_sd + pc1aN_4_sd)		;residual
 ; printVarSummary(x_sdN)
 
  dT_dt_terms_sd1N = (/pc1qN_4_sd,pc1eN_4_sd,pc1aN_4_sd,pc1dTN_4_sd,x_sdN/)
  copy_VarCoords_n(canonical_dTP_s,dT_dt_terms_sd1N,(/0,1/))
  dT_dt_terms_sd1N!0 = "lev"
  dT_dt_terms_sd1N!1 = "month"
  printVarSummary(dT_dt_terms_sd1N)

;=================================================================================================================================
   wks = gsn_open_wks("eps","Fig_3_oras5")         ; send graphics to PNG file
;***************************************************************************************************************************************************
   plot = new(2,graphic)                ; create graphic array
;***************************************************************************************************************************************************
   r2                             = True                        ; plot mods desired
   r2@gsnDraw                     = False                       ; do not draw
   r2@gsnFrame                    = False                       ; do not advance frame
   r2@tiMainString                = ""                       ; title
   r2@tiXAxisFontHeightF          = 0.015
   r2@tiYAxisFontHeightF          = 0.015
   r2@tiYAxisString               = "~F18~s~F21~T/~F18~s~F21~t~N~ terms (~F10~e~S~~F21~-6~N~ K s~S~-1~N~)"                  ; yaxis
   r2@tiXAxisString               = "Calendar month"                  ; yaxis
  
   r2@tmXBLabelFontHeightF        = 0.017
   r2@tmYLLabelFontHeightF        = 0.017
   r2@trYLog                      = False                ; log scaling
   r2@trXLog                      = False               ; log scaling

   r2@trYMinF                     =  -0.6                 ; manually set lower limit
   r2@trYMaxF                     =   0.6           

   r2@trXMinF                     = min(xlabel)        ; set minimum X-axis value
   r2@trXMaxF                     = max(xlabel)        ; set maximum X-axis value
   r2@xyLineColors                = xy_colors_thermocline
  ; r2@xyLineColors        =	 (/"red","slateblue","deepskyblue","darksalmon","black","deeppink4","yellow2","blue","seagreen3"/)
   r2@tmXTBorderOn                =  True
   r2@tmXTOn                      =  True
      
   r2@tmYLFormat                  = "@*+^sg"	; use "e" symbol, don't force 0's (see documentation)
  
   r2@vpWidthF                    = 0.55
   r2@vpHeightF                   = 0.30
 
   r2@tmYRMajorOutwardLengthF     = 0.0               ; draw tickmarks inward
   r2@tmYRMinorOutwardLengthF     = 0.0               ; draw minor ticsk inward
   r2@tmXTOn                      = True ; Turn off tickmarks and labels 

   r2@xyLineThicknesses           = (/2.,2.,2.,2.,2.,2.,2.,2.,2./)       ; Define line thicknesses 
   r2@xyDashPatterns              = (/0.,0.,0.,0.,0.,0.,0.,0.,0./)           ; Dash patterns 
   r2@tmYLLabelFontColor          = "foreground"   ; and right axes are same.
   r2@tiYAxisFontColor            = "foreground"
  
 
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
   r2@pmLegendParallelPosF   = 0.80                 ; move units right
   r2@pmLegendOrthogonalPosF = -0.51                ; move units down
   r2@pmLegendWidthF         = 0.12                ; Change width and
   r2@pmLegendHeightF        = 0.10                ; height of legend.
   r2@lgLabelFontHeightF     = 0.010                 ; change font height
   r2@lgTitleOn              = False               ; turn on legend title
   r2@lgTitleString          = "SST indices"           ; create legend title
   r2@lgTitleFontHeightF     = .003                ; font of legend title
   r2@lgPerimOn              = False               ; turn off box around
   r2@xyExplicitLegendLabels = (/"~F10~Q~B~~F21~net~N~/~F8~r~F10~C~B~~F21~w~N~~F10~h","1/~F10~h~F21~[~F10~T -T~B~~F21~-h~N~].~F10~w~B~~F21~e~N~",\
				"~F10~u~F18~s~F10~T/~F18~s~F10~x~N~ + v~F18~s~F10~T/~F18~s~F10~y~N~","~F18~s~F10~T/~F18~s~F10~t~N~","~F18~n~N~"/)        
   
   r2@tmXTOn = False ; Turn off top tickmarks and labels
  ; r2@tmYROn = False ; Turn off left tickmarks and labels
   				
    r2@tiXAxisString = "Calendar month"                  ; yaxis
    
;---For higher quality large images
  r2@tiYAxisFontThicknessF  = 1.
  r2@tmXBMajorThicknessF    = 1.
  r2@tmXBMinorThicknessF    = 1.
  r2@tmYLMajorThicknessF    = 1.
  r2@tmYLMinorThicknessF    = 1.
  r2@tmYLLabelFontThicknessF= 1.
  r2@tmBorderThicknessF     = 1.
  
;----------------------------------------------------------------------------------------------------------------   
   plotAS0  =  gsn_csm_xy(wks,xlabel,dT_dt_terms1(lev|0:,month|:),r2)
;*************************************************************************************************************** 
   plotDS0  =  gsn_csm_xy(wks,xlabel,dT_dt_terms1N(lev|0:,month|:),r2)
;*************************************************************************************************************** 
;===============================================================================================================
  ; r2@trYMinF             =  -0.2                 ; manually set lower limit
  ; r2@trYMaxF             =   0.2           

   r2@pmLegendDisplayMode    = "Always"            ; turn on legend
   r2@xyDashPatterns      = (/0.,0.,0.,0.,0.,0.,0.,0.,0./)           ; Dash patterns 
 
   r2@pmLegendDisplayMode    = "Never"            ; turn on legend
 

 ;  delete(r2@xyLineColors)
   r2@xyMarkLineMode	 =   "Markers"
   r2@xyMarkers     =  (/16.,16.,16.,16.,16.,16.,16.,16.,16./)
   r2@xyMarkerColors =  xy_colors_thermocline
   r2@xyMarkerSizeF   = 0.010
   r2@pmLegendDisplayMode    = "Never"            ; turn on legend
   r2@gsnStringFontHeightF = 0.015
   r2@gsnLeftString  		= "Canonical Events"
 ;  r2@gsnCenterString  = "Convective Events"  
;***************************************************************************************************************************************************  
   plotBS0  =  gsn_csm_xy(wks,xlabel,dT_dt_terms_sd1(lev|0:,month|:),r2)
;*************************************************************************************************************** 
   r2@gsnLeftString  		= "Non-Canonical Events"
   plotBSS0  =  gsn_csm_xy(wks,xlabel,dT_dt_terms_sd1N(lev|0:,month|:),r2)
;***************************************************************************************************************************************************  
  overlay(plotAS0,plotBS0)
  plot(0) = plotAS0
;*************************************************************************************************************** 
;=============================================================================================================== 
  overlay(plotDS0,plotBSS0)
  plot(1) = plotDS0
;===============================================================================================================
;***************************************************************************************************************************************************
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
  resP                                  = True         ; modify the panel plot
  resP@gsnDraw                          = True      ; don't draw yet
  resP@gsnFrame                         = False        ; don't advance frame yet
  resP@gsnPanelLabelBar    	        = False       ; add common colorbar
  resP@gsnPanelYWhiteSpacePercent      	= 0.0    ; default is 1.0
  resP@gsnPanelXWhiteSpacePercent      	= 0.0    ; default is 1.0
  resP@gsnPanelTop                      = 1.0
  resP@gsnPanelBottom                   = 0.51
  resP@gsnPanelLeft                     = 0.01
  resP@gsnPanelRight                    = 0.99

  resP@lbTopMarginF                     = 0.5
  resP@gsnMaximize                      = True         ; large format
  resP@gsnPaperOrientation              = "portrait"   ; force portrait
  resP@gsnLeftString 			= ""
  resP@gsnPanelFigureStrings            = (/"ORAS5","ORAS5","c","d","e","f"/) ; add strings to panel
  resP@amJust                           = "TopLeft"
  resP@gsnPanelFigureStringsFontHeightF =  0.012
  resP@gsnPanelFigureStringsFont        =  "Helvetica-Bold"

  resP@txFont                           = "Helvetica-Bold"

  gsn_panel(wks,plot(0:),(/2,1/),resP)     ; now draw as one plot 
    
  frame(wks)
end
