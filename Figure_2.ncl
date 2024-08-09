begin
;
;======================================
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
  
  xy_colors_dheating        = (/"black"/)
  p90 = 0.10
  
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
  ;printVarSummary(hadi)       ; variable overview
		
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
  delete(hadiClm)
  ;printMinMax(hadi,True)
;=========================
;----------------------------
  nyears = dimsizes(hadi&time)/nmos
  print(nyears)

  hadi_s     = onedtond(ndtooned(hadi),(/nyears,nmos,dimsizes(hadi&lat),dimsizes(hadi&lon)/))
  copy_VarCoords(hadi(0,:,:),hadi_s(0,0,:,:))                             ; copy coordinate arrays
  hadi_s!0 = "time"
  hadi_s&time = fspan(yrStrt,yrLast,nyears)
  hadi_s!1 = "month"
  printVarSummary(hadi_s)

  hadi_s = dtrend_msg_n(hadi_s&time,hadi_s,True,False,0)        ;removes linear trends
 ; printVarSummary(hadi_s)
  ;printMinMax(hadi_s,True)
;---------------------------------------------------------------------------
  ahadi_s = wgt_areaave_Wrap(hadi_s(time|:,month|:,{lon|-20:0}, {lat|-3:3}),1.0, 1.0, 0)        ;atl3 index
 ; printVarSummary(ahadi_s)
  ;printMinMax(ahadi_s,True)
;----------------------------
  delete(hadi)
;  delete(hadi_s)

  jja_ahadi = dim_avg_n_Wrap(ahadi_s(time|:,month|5:7),1)		;JJA, ncl counts from zero
 ; print(jja_ahadi)

  jja_ahadi = dim_standardize_n(jja_ahadi, 0, 0)  ;Normalise deviations
 ; print(jja_ahadi)
;----------------------------------------------------------------------------
; Positive cases
;============================================================================
  n_atlP = dim_num_n(jja_ahadi.ge.thres,0)        ;Number of cases
 ; print(n_atlP)

  valsP = ind(jja_ahadi.ge.thres)
  print(valsP)  
  ;print(valsP+yrStrt)		;the years
;----------------------------------------------------------------------------
; Negative cases,
; note the "-" sign attached to the threshold here
;============================================================================
  n_atlN = dim_num_n(jja_ahadi.le.-thres,0)        ;Number of cases
;  print(n_atlN)

  valsN = ind(jja_ahadi.le.-thres)
  ;print(valsN)
;  print(valsN+yrStrt)             ;the years 
;==============================================================
;==============================================================
  dimZ = dimsizes(hadi_s)
 ; printVarSummary(dimZ)
      
        saodP_cases_maps  = new((/n_atlP,dimZ(1),dimZ(2),dimZ(3)/),float,hadi_s@_FillValue)
        printVarSummary(saodP_cases_maps)
    
        saodP_cases_maps = hadi_s(valsP,:,:,:)
        printVarSummary(saodP_cases_maps)
        
        sstb1p = (saodP_cases_maps(time|2,month|:,lat|:,lon|:))	;now based on the basin average
       ; printVarSummary(sstb1p)    

        sstb2p = (saodP_cases_maps(time|3,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(sstb2p)    
  
        sstb3p = (saodP_cases_maps(time|4,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(sstb3p)    
  
        sstb4p = (saodP_cases_maps(time|6,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(sstb4p) 
  
        sstb5p = (saodP_cases_maps(time|7,month|:,lat|:,lon|:))	;now based on the basin average
       ; printVarSummary(sstb5p) 
        
       sst_fp = (/sstb1p,sstb2p,sstb3p,sstb4p,sstb5p/)       
       copy_VarCoords_n(hadi_s,sst_fp,(/1,2,3/))
       sst_fp!0 = "time"
       printVarSummary(sst_fp)    

;==============================================================
;==============================================================     
        saodN_cases_maps  = new((/n_atlN,dimZ(1),dimZ(2),dimZ(3)/),float,hadi_s@_FillValue)
       ; printVarSummary(saodN_cases_maps)
    
        saodN_cases_maps = hadi_s(valsN,:,:,:)
      ;  printVarSummary(saodN_cases_maps)
        
        sstb1n = (saodN_cases_maps(time|0,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(sstb1n)    
  
        sstb2n = (saodN_cases_maps(time|1,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(sstb2n)    
  
        sstb3n = (saodN_cases_maps(time|2,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(sstb3n)    
  
        sstb4n = (saodN_cases_maps(time|5,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(sstb4n) 
       
       sst_fn = (/sstb1n,sstb2n,sstb3n,sstb4n/)       
       copy_VarCoords_n(hadi_s,sst_fn,(/1,2,3/))
       sst_fn!0 = "time"    
       printVarSummary(sst_fn)    
;==============================================================
  pc1SST_C = (dim_avg_n_Wrap(sst_fp,0) - dim_avg_n_Wrap(sst_fn,0))/(2.)                       ; Difference
  copy_VarCoords(hadi_s(0,:,:,:),pc1SST_C)                             ; copy coordinate arrays 
  printVarSummary(pc1SST_C)
;===================================== 
  pc1SSTd_C = pc1SST_C  ; Transfer metdata for plotting sig. values        

  pc1SSTd_C = ttest((dim_avg_n_Wrap(sst_fp,0)),(dim_stddev_n_Wrap(sst_fp,0))^2,5\
               ,(dim_avg_n_Wrap(sst_fn,0)),(dim_stddev_n_Wrap(sst_fn,0))^2,4,True, False)
  pc1SSTd_C@long_name = "Probability: difference between means"
 ; printVarSummary(pc1SSTd_C)

  pc1SSTd_C = where(abs(pc1SSTd_C).le.p90,pc1SSTd_C,pc1SSTd_C@_FillValue) ; Sig.u winds, where UWND or VWND is significanlt
  printVarSummary(pc1SSTd_C)
  printMinMax(pc1SSTd_C,True)
    
  d_SSTP_C = dim_avg_n_Wrap(pc1SST_C(month|:,{lat|-3:3}, lon|:),1) 
  d_SSTdP_C = dim_avg_n_Wrap(pc1SSTd_C(month|:,{lat|-3:3}, lon|:),1)
  delete([/sstb1n, sstb2n, sstb3n, sstb4n/])
  delete([/sstb1p, sstb2p, sstb3p, sstb4p, sstb5p/]) 
;=============================================================================================================================================================      
        sstb1p1 = (saodP_cases_maps(time|0,month|:,lat|:,lon|:))	;now based on the basin average
       ; printVarSummary(sstb1p1)    
  
        sstb2p1 = (saodP_cases_maps(time|1,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(sstb2p1)    
  
        sstb3p1 = (saodP_cases_maps(time|5,month|:,lat|:,lon|:))	;now based on the basin average
       ; printVarSummary(sstb3p1)    
  
        sstb4p1 = (saodP_cases_maps(time|8,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(sstb4p1) 
         
       sst_fp1 = (/sstb1p1,sstb2p1,sstb3p1,sstb4p1/)        
       copy_VarCoords_n(hadi_s,sst_fp1,(/1,2,3/))
       sst_fp1!0 = "time"    
       printVarSummary(sst_fp1) 
;===================================================================
;===========================================================================================================================       
        sstb1n1 = (saodN_cases_maps(time|3,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(sstb1n1)    
  
        sstb2n1 = (saodN_cases_maps(time|4,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(sstb2n1)    
  
        sstb3n1 = (saodN_cases_maps(time|6,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(sstb3n1)    
  
        sstb4n1 = (saodN_cases_maps(time|7,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(sstb4n1) 
        
       sst_fn1 = (/sstb1n1,sstb2n1,sstb3n1,sstb4n1/)       
       copy_VarCoords_n(hadi_s,sst_fn1,(/1,2,3/))
       sst_fn1!0 = "time"    
    ;   printVarSummary(sst_fn1) 
;==============================================================
  pc1SST_N = (dim_avg_n_Wrap(sst_fp1,0) - dim_avg_n_Wrap(sst_fn1,0))/(2.)
  copy_VarCoords(hadi_s(0,:,:,:),pc1SST_N)                             ; copy coordinate arrays 
  printVarSummary(pc1SST_N)
;===================================== 
  pc1SSTd_N = pc1SST_N  ; Transfer metdata for plotting sig. values        

  pc1SSTd_N = ttest((dim_avg_n_Wrap(sst_fp1,0)),(dim_stddev_n_Wrap(sst_fp1,0))^2,4\
               ,(dim_avg_n_Wrap(sst_fn1,0)),(dim_stddev_n_Wrap(sst_fn1,0))^2,4,True, False)
  pc1SSTd_N@long_name = "Probability: difference between means"
 ; printVarSummary(pc1SSTd_N)

  pc1SSTd_N = where(abs(pc1SSTd_N).le.p90,pc1SSTd_N,pc1SSTd_N@_FillValue) ; Sig.u winds, where UWND or VWND is significanlt
  ;print("==============================================================")  
 ; printVarSummary(pc1SSTd_N)
  
  d_SSTP_N = dim_avg_n_Wrap(pc1SST_N(month|:,{lat|-3:3}, lon|:),1)
  
  d_SSTdP_N = dim_avg_n_Wrap(pc1SSTd_N(month|:,{lat|-3:3}, lon|:),1)
  
  delete([/sstb1p1, sstb2p1, sstb3p1, sstb4p1,saodN_cases_maps, pc1SSTd_N, sst_fn, sst_fp/])
  delete([/sstb1n1, sstb2n1, sstb3n1, sstb4n1,saodP_cases_maps, pc1SST_N, sst_fp1, sst_fn1/])
  delete(hadi_s)
;==============================================================         
;==============================================================
  f      = addfile("/home/gyuk/Downloads/GYUK/GEOMAR/Final_work/potential_temperature.nc", "r")
  TIME   = f->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  pot_t   = f->pottmp(iYYYY,:,:,:)
  printVarSummary(pot_t)       ; variable overview
  delete(f)             ;       Not needed anymore
  delete(TIME)
  delete(YYYY)
  delete(iYYYY)
;***************************************************************************************************************
;***************************************************************************************************************
; Interpolate to specific [constant] TEMP levels
; The default returned vertical coordinate is Z_T but change to 'tlev'
;***************************************************************************************************************
  level   = pot_t&level
  
  tlev    = (/20/)    ; same units [here, C] as TEMP
  tlev@units = pot_t@units
  tlev!0     = "tlev"
 
  depth   = conform(pot_t, level, 1)
  copy_VarCoords(pot_t, depth)
  printVarSummary(depth)
;=================================================================== 
 ; depth_xaxis = depth&lon
 ; copy_VarCoords_n(depth,depth_xaxis,(/3/))
 ; depth_xaxis1 = depth_xaxis({lon|-41:11})
 ; printVarSummary(depth_xaxis1)

  isot    = int2p_n_Wrap(pot_t,depth, tlev, 0,1)
  printVarSummary(isot)  
  printMinMax(isot, True)   

  iso1 = dim_avg_n_Wrap(isot(time|:,tlev|:,lat|:, lon|:),0)
  printVarSummary(iso1)
  
  iso2 = dim_avg_n_Wrap(iso1(tlev|0,{lat|-3:3}, lon|:),0)
  printVarSummary(iso2)
  printMinMax(iso2, True)    

;========================================================================================
  pot_tClm      = clmMonTLLL(pot_t)
  printVarSummary(pot_tClm)                   ; (12,nlat,nlon)   
  pot_t         = calcMonAnomTLLL(pot_t,pot_tClm)    ; replace with anonamlies
  pot_t@long_name = "ANOMALIES: "+pot_t@long_name
  delete(pot_tClm)
;***************************************************************************************************************
;***************************************************************************************************************  
  nyears = dimsizes(pot_t&time)/nmos
  ;print(nyears)
  pot_t_s     = onedtond(ndtooned(pot_t),(/nyears,nmos,dimsizes(pot_t&level),dimsizes(pot_t&lat),dimsizes(pot_t&lon)/))
  copy_VarCoords(pot_t(0,:,:,:),pot_t_s(0,0,:,:,:))                             ; copy coordinate arrays
  pot_t_s!0 = "time"
  pot_t_s&time = fspan(yrStrt,yrLast,nyears)
  pot_t_s!1 = "month"
  printVarSummary(pot_t_s)
  printMinMax(pot_t_s,True)

  pot_t_s = dtrend_msg_n(pot_t_s&time,pot_t_s,True,False,0)        ;removes linear trends
 ; printVarSummary(pot_t_s)
 ; printMinMax(pot_t_s,True)
;=================== Canonical==============================================================================
;***************************************************************************************************************
  ;  delete(dimZ)       
    dimZ = dimsizes(pot_t_s)
    ; printVarSummary(dimZ)
    
        pot_tdP_cases_maps  = new((/n_atlP,dimZ(1),dimZ(2),dimZ(3),dimZ(4)/),float,pot_t_s@_FillValue)
        printVarSummary(pot_tdP_cases_maps)
    
        pot_tdP_cases_maps = pot_t_s(valsP,:,:,:,:)
        printVarSummary(pot_tdP_cases_maps)
        
        pot_tb1p = (pot_tdP_cases_maps(time|1,month|:,level|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(pot_tb1p)    
  
        pot_tb2p = (pot_tdP_cases_maps(time|3,month|:,level|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(pot_tb2p)    
  
        pot_tb3p = (pot_tdP_cases_maps(time|4,month|:,level|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(pot_tb3p)    
  
        pot_tb4p = (pot_tdP_cases_maps(time|6,month|:,level|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(pot_tb4p) 
  
        pot_tb5p = (pot_tdP_cases_maps(time|7,month|:,level|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(pot_tb5p) 
         
       pot_t_fp = (/pot_tb1p,pot_tb2p,pot_tb3p,pot_tb4p,pot_tb5p/)       
       copy_VarCoords_n(pot_t_s,pot_t_fp,(/1,2,3,4/))
       pot_t_fp!0 = "time"    
       printVarSummary(pot_t_fp) 
;==============================================================
;==============================================================     
        pot_tdN_cases_maps  = new((/n_atlN,dimZ(1),dimZ(2),dimZ(3),dimZ(4)/),float,pot_t_s@_FillValue)
        printVarSummary(pot_tdN_cases_maps)
    
        pot_tdN_cases_maps = pot_t_s(valsN,:,:,:,:)
        printVarSummary(pot_tdN_cases_maps)
        
        pot_tb1n = (pot_tdN_cases_maps(time|0,month|:,level|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(pot_tb1n)    
  
        pot_tb2n = (pot_tdN_cases_maps(time|1,month|:,level|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(pot_tb2n)    
  
        pot_tb3n = (pot_tdN_cases_maps(time|2,month|:,level|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(pot_tb3n)    
  
        pot_tb4n = (pot_tdN_cases_maps(time|5,month|:,level|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(pot_tb4n) 
       
       pot_t_fn = (/pot_tb1n,pot_tb2n,pot_tb3n,pot_tb4n/)        
       copy_VarCoords_n(pot_t_s,pot_t_fn,(/1,2,3,4/))
       pot_t_fn!0 = "time"    
       printVarSummary(pot_t_fn) 
;==============================================================
  pc1pot_t_C = (dim_avg_n_Wrap(pot_t_fp,0) - dim_avg_n_Wrap(pot_t_fn,0))/(2.)                       ; Difference
  copy_VarCoords(pot_t_s(0,:,:,:,:),pc1pot_t_C)                             ; copy coordinate arrays 
  printVarSummary(pc1pot_t_C)
 
  d_pot_tP_C = dim_avg_n_Wrap(pc1pot_t_C(month|5:7,level|:,lat|:, lon|:),0)
  printVarSummary(d_pot_tP_C)
;=====================================================================================================
  pc1pot_td_C = pc1pot_t_C  ; Transfer metdata for plotting sig. values        

  pc1pot_td_C = ttest((dim_avg_n_Wrap(pot_t_fp,0)),(dim_stddev_n_Wrap(pot_t_fp,0))^2,5\
               ,(dim_avg_n_Wrap(pot_t_fn,0)),(dim_stddev_n_Wrap(pot_t_fn,0))^2,4,True, False)
  pc1pot_td_C@long_name = "Probability: difference between means"
  printVarSummary(pc1pot_td_C)

  pc1pot_td_C = where(abs(pc1pot_td_C).le.p90,pc1pot_td_C,pc1pot_td_C@_FillValue)
  printVarSummary(pc1pot_td_C)
 
  sub_SST_C = dim_avg_n_Wrap(pc1pot_td_C(month|5:7,level|:,lat|:, lon|:),0)            ; SSTA difference in events 
;===================================================================================================== 
  sub_C = dim_avg_n_Wrap(d_pot_tP_C(level|:,{lat|-3:3}, lon|:),1)
  printVarSummary(sub_C)
  sub_SSTdP_C = dim_avg_n_Wrap(sub_SST_C(level|:,{lat|-3:3}, lon|:),1)  
  printVarSummary(sub_SSTdP_C)
;==============================================================
;==========Non-canonical========================================
;==============================================================
        pot_tb1p1 = (pot_tdP_cases_maps(time|0,month|:,level|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(pot_tb1p)    
  
        pot_tb2p1 = (pot_tdP_cases_maps(time|2,month|:,level|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(pot_tb2p1)    
  
        pot_tb3p1 = (pot_tdP_cases_maps(time|5,month|:,level|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(pot_tb3p1)    
  
        pot_tb4p1 = (pot_tdP_cases_maps(time|8,month|:,level|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(pot_tb4p1) 
       
       pot_t_fp1 = (/pot_tb1p1,pot_tb2p1,pot_tb3p1,pot_tb4p1/)       
       copy_VarCoords_n(pot_t_s,pot_t_fp1,(/1,2,3,4/))
       pot_t_fp1!0 = "time"   
       printVarSummary(pot_t_fp1) 
;===================================================================
;==============================================================        
        pot_tb1n1 = (pot_tdN_cases_maps(time|3,month|:,level|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(pot_tb1n1)    
  
        pot_tb2n1 = (pot_tdN_cases_maps(time|4,month|:,level|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(pot_tb2n1)    
  
        pot_tb3n1 = (pot_tdN_cases_maps(time|6,month|:,level|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(pot_tb3n1)    
  
        pot_tb4n1 = (pot_tdN_cases_maps(time|7,month|:,level|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(pot_tb4n1) 
      
       pot_t_fn1 = (/pot_tb1n1,pot_tb2n1,pot_tb3n1,pot_tb4n1/)       
       copy_VarCoords_n(pot_t_s,pot_t_fn1,(/1,2,3,4/))
       pot_t_fn1!0 = "time"   
       printVarSummary(pot_t_fn1) 
;==================================================================================================================

  pc1pot_t_N = (dim_avg_n_Wrap(pot_t_fp1,0) - dim_avg_n_Wrap(pot_t_fn1,0))/(2.)
  copy_VarCoords(pot_t_s(0,:,:,:,:),pc1pot_t_N)                             ; copy coordinate arrays 
  printVarSummary(pc1pot_t_N)

  d_pot_tP_N = dim_avg_n_Wrap(pc1pot_t_N(month|5:7,level|:,lat|:, lon|:),0)
  printVarSummary(d_pot_tP_N)
;***************************************************************************************************************
;=====================================================================================================
  pc1pot_td_N = pc1pot_t_N  ; Transfer metdata for plotting sig. values        

  pc1pot_td_N = ttest((dim_avg_n_Wrap(pot_t_fp1,0)),(dim_stddev_n_Wrap(pot_t_fp1,0))^2,5\
               ,(dim_avg_n_Wrap(pot_t_fn1,0)),(dim_stddev_n_Wrap(pot_t_fn1,0))^2,4,True, False)
  pc1pot_td_N@long_name = "Probability: difference between means"
  printVarSummary(pc1pot_td_N)

  pc1pot_td_N = where(abs(pc1pot_td_N).le.p90,pc1pot_td_N,pc1pot_td_N@_FillValue)
  printVarSummary(pc1pot_td_N)
 
  sub_SST_N = dim_avg_n_Wrap(pc1pot_td_N(month|5:7,level|:,lat|:, lon|:),0)            ; SSTA difference in events 
;===================================================================================================== 
  sub_N = dim_avg_n_Wrap(d_pot_tP_N(level|:,{lat|-3:3}, lon|:),1)
  printVarSummary(sub_N)  
  
  sub_SSTdP_N = dim_avg_n_Wrap(sub_SST_N(level|:,{lat|-3:3}, lon|:),1)  
  printVarSummary(sub_SSTdP_N)
;========================================================================================================
;=========================================================================================================
  f      = addfile("/home/gyuk/Downloads/GYUK/GEOMAR/Final_work/sea_surface_height.nc", "r")
  TIME   = f->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  ssh   = f->sshg(iYYYY,:,:)
  printVarSummary(ssh)       ; variable overview
  delete(f)             ;       Not needed anymore
  delete(TIME)
  delete(YYYY)
  delete(iYYYY)
  
  ssh = lonFlip(ssh)
  ssh = ssh*10.     ;converts to cm
  
  ssh!1 = "lat"
  ssh!2 = "lon"
;==========================================================================================================
  sshClm      = clmMonTLL(ssh)
  printVarSummary(sshClm)                   ; (12,nlat,nlon)   
  ssh         = calcMonAnomTLL(ssh,sshClm)    ; replace with anonamlies
  ssh@long_name = "ANOMALIES: "+ssh@long_name
  delete(sshClm)
;----------------------------
  nyears = dimsizes(ssh&time)/nmos
  ;print(nyears)
  ssh_s     = onedtond(ndtooned(ssh),(/nyears,nmos,dimsizes(ssh&lat),dimsizes(ssh&lon)/))
  copy_VarCoords(ssh(0,:,:),ssh_s(0,0,:,:))                             ; copy coordinate arrays
  ssh_s!0 = "time"
  ssh_s&time = fspan(yrStrt,yrLast,nyears)
  ssh_s!1 = "month"
  printVarSummary(ssh_s)

  ssh_s = dtrend_msg_n(ssh_s&time,ssh_s,True,False,0)        ;removes linear trends
  ;printVarSummary(ssh_s)
  ;printMinMax(ssh_s,True)
;=================== Canonical================================================================================
;==============================================================
;***************************************************************************************************************
  dimZ = dimsizes(ssh_s)
 ; printVarSummary(dimZ)
      
        sshdP_cases_maps  = new((/n_atlP,dimZ(1),dimZ(2),dimZ(3)/),float,ssh_s@_FillValue)
     ;   printVarSummary(sshdP_cases_maps)
    
        sshdP_cases_maps = ssh_s(valsP,:,:,:)
     ;   printVarSummary(sshdP_cases_maps)
        
        sshb1p = (sshdP_cases_maps(time|2,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(sshb1p)    
  
        sshb2p = (sshdP_cases_maps(time|3,month|:,lat|:,lon|:))	;now based on the basin average
     ;   printVarSummary(sshb2p)    
  
        sshb3p = (sshdP_cases_maps(time|4,month|:,lat|:,lon|:))	;now based on the basin average
     ;   printVarSummary(sshb3p)    
  
       sshb4p = (sshdP_cases_maps(time|6,month|:,lat|:,lon|:))	;now based on the basin average
     ;   printVarSummary(sshb4p) 
  
        sshb5p = (sshdP_cases_maps(time|7,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(sshb5p)         
       ssh_fp = (/sshb1p,sshb2p,sshb3p,sshb4p,sshb5p/)       
       copy_VarCoords_n(ssh_s,ssh_fp,(/1,2,3/))
       ssh_fp!0 = "time"
       printVarSummary(ssh_fp) 
;==============================================================      
        sshdN_cases_maps  = new((/n_atlN,dimZ(1),dimZ(2),dimZ(3)/),float,ssh_s@_FillValue)
      ;  printVarSummary(sshdN_cases_maps)
    
        sshdN_cases_maps = ssh_s(valsN,:,:,:)
       ; printVarSummary(sshdN_cases_maps)
        
        sshb1n = (sshdN_cases_maps(time|0,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(sshb1n)    
  
        sshb2n = (sshdN_cases_maps(time|1,month|:,lat|:,lon|:))	;now based on the basin average
       ; printVarSummary(sshb2n)    
  
        sshb3n = (sshdN_cases_maps(time|2,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(sshb3n)    
  
        sshb4n = (sshdN_cases_maps(time|5,month|:,lat|:,lon|:))	;now based on the basin average
       ; printVarSummary(sshb4n) 
       
       ssh_fn = (/sshb1n,sshb2n,sshb3n,sshb4n/)        
       copy_VarCoords_n(ssh_s,ssh_fn,(/1,2,3/))
       ssh_fn!0 = "time"
       printVarSummary(ssh_fn) 
;==============================================================

  pc1ssh_C = (dim_avg_n_Wrap(ssh_fp,0) - dim_avg_n_Wrap(ssh_fn,0))/(2.)                       ; Difference
  copy_VarCoords(ssh_s(0,:,:,:),pc1ssh_C)                             ; copy coordinate arrays 
  printVarSummary(pc1ssh_C)
;===================================== 
  pc1sshd_C = pc1ssh_C  ; Transfer metdata for plotting sig. values        

  pc1sshd_C = ttest((dim_avg_n_Wrap(ssh_fp,0)),(dim_stddev_n_Wrap(ssh_fp,0))^2,5\
               ,(dim_avg_n_Wrap(ssh_fn,0)),(dim_stddev_n_Wrap(ssh_fn,0))^2,4,True, False)
  pc1sshd_C@long_name = "Probability: difference between means"
  printVarSummary(pc1sshd_C)

  pc1sshd_C = where(abs(pc1sshd_C).le.p90,pc1sshd_C,pc1sshd_C@_FillValue) ; Sig.u winds, where UWND or VWND is significanlt
  ;print("==============================================================")  
  printVarSummary(pc1sshd_C) 
  
  d_sshP_C = dim_avg_n_Wrap(pc1ssh_C(month|:,{lat|-3:3}, lon|:),1)
  d_sshdP_C = dim_avg_n_Wrap(pc1sshd_C(month|:,{lat|-3:3}, lon|:),1)
  
  delete([/sshb1p,sshb2p,sshb3p,sshb4p,sshb5p,pc1sshd_C/])
  delete([/sshb1n, sshb2n, sshb3n, sshb4n/])
;==============================================================
;==========Non-canonical=========================================================================
;==============================================================
        sshb1p1 = (sshdP_cases_maps(time|0,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(sshb1p)    
  
        sshb2p1 = (sshdP_cases_maps(time|1,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(sshb2p1)    
  
        sshb3p1 = (sshdP_cases_maps(time|5,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(sshb3p1)    
  
        sshb4p1 = (sshdP_cases_maps(time|8,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(sshb4p1)   
       ssh_fp1 = (/sshb1p1,sshb2p1,sshb3p1,sshb4p1/)       
       copy_VarCoords_n(ssh_s,ssh_fp1,(/1,2,3/))
       ssh_fp1!0 = "time"   
       printVarSummary(ssh_fp1) 
;===================================================================
;==============================================================       
        sshb1n1 = (sshdN_cases_maps(time|3,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(sshb1n1)    
  
        sshb2n1 = (sshdN_cases_maps(time|4,month|:,lat|:,lon|:))	;now based on the basin average
       ; printVarSummary(sshb2n1)    
  
        sshb3n1 = (sshdN_cases_maps(time|6,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(sshb3n1)    
  
        sshb4n1 = (sshdN_cases_maps(time|7,month|:,lat|:,lon|:))	;now based on the basin average
      ;  printVarSummary(sshb4n1) 
       
       ssh_fn1 = (/sshb1n1,sshb2n1,sshb3n1,sshb4n1/)       
       copy_VarCoords_n(ssh_s,ssh_fn1,(/1,2,3/))
       ssh_fn1!0 = "time"   
       printVarSummary(ssh_fn1) 
;==============================================================
  pc1ssh_N = (dim_avg_n_Wrap(ssh_fp1,0) - dim_avg_n_Wrap(ssh_fn1,0))/(2.)
  copy_VarCoords(ssh_s(0,:,:,:),pc1ssh_N)                             ; copy coordinate arrays 
  printVarSummary(pc1ssh_N)
;===================================== 
  pc1sshd_N = pc1ssh_N  ; Transfer metdata for plotting sig. values        

  pc1sshd_N = ttest((dim_avg_n_Wrap(ssh_fp1,0)),(dim_stddev_n_Wrap(ssh_fp1,0))^2,4\
               ,(dim_avg_n_Wrap(ssh_fn1,0)),(dim_stddev_n_Wrap(ssh_fn1,0))^2,4,True, False)
  pc1sshd_N@long_name = "Probability: difference between means"
  printVarSummary(pc1sshd_N)

  pc1sshd_N = where(abs(pc1sshd_N).le.p90,pc1sshd_N,pc1sshd_N@_FillValue) ; Sig.u winds, where UWND or VWND is significanlt
  ;print("==============================================================")  
  printVarSummary(pc1sshd_N)

  d_sshP_N = dim_avg_n_Wrap(pc1ssh_N(month|:,{lat|-3:3}, lon|:),1)
  
  d_sshdP_N = dim_avg_n_Wrap(pc1sshd_N(month|:,{lat|-3:3}, lon|:),1)
  
  delete([/sshb1p1,sshb2p1,sshb3p1,sshb4p1,sshdN_cases_maps,ssh_fp1,ssh_fn/])
  delete([/sshb1n1,sshb2n1,sshb3n1,sshb4n1,sshdP_cases_maps,ssh_fn1,ssh_fp/])
;==============================================================
;==============================================================
  f      = addfile("/home/gyuk/Downloads/GYUK/GEOMAR/Final_work/precip.mon.mean.nc", "r")
  TIME   = f->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  rain   = f->precip(iYYYY,:,:)
  rain = lonFlip(rain)
  printVarSummary(rain)       ; variable overview
  delete(f)             ;       Not needed anymore
  delete(TIME)
  delete(YYYY)
  delete(iYYYY)

  rain!1 = "lat"
  rain!2 = "lon"
 
; printMinMax(rain,True)
;==========================
  rainClm      = clmMonTLL(rain)
  printVarSummary(rainClm)                   ; (12,nlat,nlon)   
  rain         = calcMonAnomTLL(rain,rainClm)    ; replace with anonamlies
  rain@long_name = "ANOMALIES: "+rain@long_name
  delete(rainClm)
;----------------------------
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
;==============================================================
;=================== Canonical Precipitation=================================
;==============================================================================================================
;***************************************************************************************************************
  delete(dimZ)
  dimZ = dimsizes(rain_s)
  printVarSummary(dimZ)
      
        raindP_cases_maps  = new((/n_atlP,dimZ(1),dimZ(2),dimZ(3)/),float,rain_s@_FillValue)
        printVarSummary(raindP_cases_maps)
    
        raindP_cases_maps = rain_s(valsP,:,:,:)
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
       copy_VarCoords_n(rain_s,rain_fp,(/1,2,3/))
       rain_fp!0 = "time"    
       printVarSummary(rain_fp) 
;===================================================================================================================     
        raindN_cases_maps  = new((/n_atlN,dimZ(1),dimZ(2),dimZ(3)/),float,rain_s@_FillValue)
        printVarSummary(raindN_cases_maps)
    
        raindN_cases_maps = rain_s(valsN,:,:,:)
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
       copy_VarCoords_n(rain_s,rain_fn,(/1,2,3/))
       rain_fn!0 = "time"    
       printVarSummary(rain_fn) 
;==============================================================
  pc1rain_C = (dim_avg_n_Wrap(rain_fp,0) - dim_avg_n_Wrap(rain_fn,0))/(2.)                       ; Difference
  copy_VarCoords(rain_s(0,:,:,:),pc1rain_C)                             ; copy coordinate arrays 
  printVarSummary(pc1rain_C)
;======================================================================================================== 
  pc1raind_C = pc1rain_C  ; Transfer metdata for plotting sig. values        

  pc1raind_C = ttest((dim_avg_n_Wrap(rain_fp,0)),(dim_stddev_n_Wrap(rain_fp,0))^2,5\
               ,(dim_avg_n_Wrap(rain_fn,0)),(dim_stddev_n_Wrap(rain_fn,0))^2,4,True, False)
  pc1raind_C@long_name = "Probability: difference between means"
  printVarSummary(pc1raind_C)

  pc1raind_C = where(abs(pc1raind_C).le.p90,pc1raind_C,pc1raind_C@_FillValue)
  ;print("==============================================================")  
  printVarSummary(pc1raind_C)
;======================================================================================================================= 
  rain_C = dim_avg_n_Wrap(pc1rain_C(month|:,{lat|-3:3}, lon|:),1)

  raindP_C = dim_avg_n_Wrap(pc1raind_C(month|:,{lat|-3:3}, lon|:),1)
;=======================================================================================================================
;=================== Non-Canonical Precipitation======================================================  
;==========================================================================================================================
        rainb1p1 = (raindP_cases_maps(time|0,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1p)    
  
        rainb2p1 = (raindP_cases_maps(time|1,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2p1)    
  
        rainb3p1 = (raindP_cases_maps(time|5,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3p1)    
  
        rainb4p1 = (raindP_cases_maps(time|8,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4p1) 
        
       rain_fp1 = (/rainb1p1,rainb2p1,rainb3p1,rainb4p1/)       
       copy_VarCoords_n(rain_s,rain_fp1,(/1,2,3/))
       rain_fp1!0 = "time"   
       printVarSummary(rain_fp1) 
;====================================================================================================================     

        rainb1n1 = (raindN_cases_maps(time|3,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb1n1)    
  
        rainb2n1 = (raindN_cases_maps(time|4,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb2n1)    
  
        rainb3n1 = (raindN_cases_maps(time|6,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb3n1)    
  
        rainb4n1 = (raindN_cases_maps(time|7,month|:,lat|:,lon|:))	;now based on the basin average
        printVarSummary(rainb4n1) 
        
       rain_fn1 = (/rainb1n1,rainb2n1,rainb3n1,rainb4n1/)       
       copy_VarCoords_n(rain_s,rain_fn1,(/1,2,3/))
       rain_fn1!0 = "time"    
       printVarSummary(rain_fn1) 
;======================================================================================================================
  pc1rain_N = (dim_avg_n_Wrap(rain_fp1,0) - dim_avg_n_Wrap(rain_fn1,0))/(2.)
  copy_VarCoords(rain_s(0,:,:,:),pc1rain_N)                             ; copy coordinate arrays 
  printVarSummary(pc1rain_N)
;====================================================================================================================== 
  pc1raind_N = pc1rain_N  ; Transfer metdata for plotting sig. values        

  pc1raind_N = ttest((dim_avg_n_Wrap(rain_fp1,0)),(dim_stddev_n_Wrap(rain_fp1,0))^2,4\
               ,(dim_avg_n_Wrap(rain_fn1,0)),(dim_stddev_n_Wrap(rain_fn1,0))^2,4,True, False)
  pc1raind_N@long_name = "Probability: difference between means"
  printVarSummary(pc1raind_N)

  pc1raind_N = where(abs(pc1raind_N).le.p90,pc1raind_N,pc1raind_N@_FillValue)
  ;print("==============================================================")  
  printVarSummary(pc1raind_N)
 
  rain_N  = dim_avg_n_Wrap(pc1rain_N(month|:,{lat|-3:3}, lon|:),1)

  raindP_N = dim_avg_n_Wrap(pc1raind_N(month|:,{lat|-3:3}, lon|:),1)
;***************************************************************************************************************
;***************************************************************************************************************  
;***************************************************************************************************************
  wks = gsn_open_wks("eps","Figure_2")         ; send graphics to PNG file
  plot = new(8,graphic)                ; create graphic array ; only needed if paneling

  res                      = True         
  res@gsnDraw              = False        ; don't draw yet
  res@gsnFrame             = False        ; don't advance frame yet
  res@trYReverse            = True
  
  res@gsnMaximize          = True  
  
  res@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels
  res@cnMinLevelValF       = -2.            ; set min contour level
  res@cnMaxLevelValF       =  2.            ; set max contour level
  res@cnLevelSpacingF      =  0.2            ; set contour spacing  
 
 ; res@pmTickMarkDisplayMode = "Always"
 
  res@tmXBMode   = "Explicit"
  res@tmXBValues = ispan(-50,10,10)
  res@tmXBLabels = (/"50~S~o~N~W","40~S~o~N~W","30~S~o~N~W","20~S~o~N~W","10~S~o~N~W","0~S~o~N~","10~S~o~N~E","20~S~o~N~E"/)
 ;res@tmXBMinorValues  = ispan(-50,20,5)
  res@tmXBLabelFontHeightF = 0.030

;  res@tmYLMode = "Explicit"
 ; res@tmYLValues = (/"0", "50","100","150","200","250","300"/)
 ; res@tmYLLabels = (/"0", "50","100","150","200","250","300"/)
  ;res@tmXBLabels = (/"-5", "-4","-3","-2","-1","0","+1","+2","+3","+4","+5","+6"/)
  res@tiYAxisString = "Ocean Depth (m)"                  ; yaxis
  res@tmYLLabelFontHeightF = 0.040

 
  res@tmXBLabelFontHeightF = 0.040
      
  res@gsnAddCyclic         = False
  res@cnFillOn             = True         ; turn on color fill
  res@cnLinesOn            = False        ; True is default
  res@cnLineLabelsOn       = False        ; True is default
  res@cnMissingValFillColor = "grey" ; comment this line out 

  res@cnFillPalette = "ncl_default"
  
  res@lbLabelBarOn         =  True       ; turn off individual lb's
  res@lbOrientation        = "vertical"
  res@lbTitleFontHeightF   = 0.020 
  res@lbLabelFontHeightF   = 0.020
  res@lbLeftMarginF         = 0.10
  res@lbRightMarginF        = 0.55
  res@lbLabelFont   		     = "Helvetica-Bold"     ; label font
  res@lbTitleString        = "Ocean Temperature (~S~o~N~C)"
  res@lbTitleAngleF        = 270.
  res@lbTitlePosition      = "Right"
  res@lbTitleDirection     = "Across"
  
  res@tmYLLabelFontHeightF = 0.040
  res@tmXBLabelFontHeightF = 0.040
  
  res@lbLabelStride        = 10
  
  res@tmXTOn = True ; Turn off top tickmarks and labels
  res@tmYROn = True ; Turn off left tickmarks and labels
  
  res@tiXAxisFontHeightF = 0.040
  res@tiYAxisFontHeightF = 0.040
  
  res@vpHeightF 	  = 0.80        ; Changes the aspect ratio
  res@vpWidthF 		  = 0.80

  ;---For higher quality large images
  res@tiYAxisFontThicknessF  = 1.
  res@tmXBMajorThicknessF    = 1.
  res@tmXBMinorThicknessF    = 1.
  res@tmYLMajorThicknessF    = 1.
  res@tmYLMinorThicknessF    = 1.
  res@tmYLLabelFontThicknessF= 1.
  res@tmBorderThicknessF     = 1.
;=============================================================================
  ress                      = True         
  ress@gsnDraw              = False        ; don't draw yet
  ress@gsnFrame             = False        ; don't advance frame yet
  ress@gsnMaximize          = True  
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

  ress@cnInfoLabelOn         = False
  ress@lbLabelBarOn          = True       ; turn off individual lb's 
 
  ress@cnLineLabelInterval   = 1           ; Default is 2
  
 ; ress@pmTickMarkDisplayMode = "Always"         ; turn on nicer tickmarks
  ress@tmXBLabelFontHeightF  = 0.040
  ress@tmYLLabelFontHeightF  = 0.040
 ; ress@tmYLLabelStride       = 2                ; label every other tickmark
  ress@tmXBLabelStride       = 2
 
  ress@tmXTOn                  = True 
  ress@tmYROn                  = True 
  ;ress@tmYRLalabelsbelsOn            = False
  ress@tmXTLabelsOn            = False
 
  ress@lbOrientation           = "vertical"
  ress@lbTitleFontHeightF      = 0.040 
  ress@lbLabelFontHeightF      = 0.040
  ress@lbLeftMarginF           = 0.10
  ress@lbRightMarginF          = 0.55
  ;ress@lbLabelOffsetF	       = 0.20
  ress@lbLabelFont   	       = "Helvetica-Bold"     ; label font
  ress@lbTitleString           = "SSTA (~S~o~N~C)"
  ;ress@lbLabelPosition	       = "Top"
  ress@lbTitleAngleF           = 270.
  ress@lbTitlePosition         = "Right"
  ress@lbTitleDirection        = "Across"
  
  ress@lbLabelStride        = 10 
  
  ress@tmXBMode   = "Explicit"
  ress@tmXBValues = ispan(-50,10,10)
  ress@tmXBLabels = (/"50~S~o~N~W","40~S~o~N~W","30~S~o~N~W","20~S~o~N~W","10~S~o~N~W","0~S~o~N~","10~S~o~N~E","20~S~o~N~E"/)
  ;ress@tmXBMinorValues  = ispan(-50,20,5)
  ress@tmXBLabelFontHeightF = 0.030

  ress@tmYLMode = "Explicit"
  ress@tmYLValues = (/"0", "1","2","3","4","5","6","7","8","9","10","11"/)
  ress@tmYLLabels = (/"J", "F","M","A","M","J","J","A","S","O","N","D"/)
  ;ress@tmXBLabels = (/"-5", "-4","-3","-2","-1","0","+1","+2","+3","+4","+5","+6"/)
  ress@tiYAxisString = "Month"                  ; yaxis
  ress@tmYLLabelFontHeightF = 0.040

  ress@tiXAxisFontHeightF = 0.040
  ress@tiYAxisFontHeightF = 0.040
  
  ress@tmXTOn = True ; Turn off top tickmarks and labels
  ress@tmYROn = True ; Turn off left tickmarks and labels
  
  ress@tiXAxisFontHeightF = 0.040
  ress@tiYAxisFontHeightF = 0.040
  
  ress@vpHeightF 	  = 0.80        ; Changes the aspect ratio
  ress@vpWidthF 	  = 0.80
 
  ;---For higher quality large images
  ress@tiYAxisFontThicknessF  = 1.
  ress@tmXBMajorThicknessF    = 1.
  ress@tmXBMinorThicknessF    = 1.
  ress@tmYLMajorThicknessF    = 1.
  ress@tmYLMinorThicknessF    = 1.
  ress@tmYLLabelFontThicknessF= 1.
  ress@tmBorderThicknessF     = 1.
;-----------------------------------------------------------------------------
  resh                      = True         
  resh@gsnDraw              = False        ; don't draw yet
  resh@gsnFrame             = False        ; don't advance frame yet
  resh@gsnMaximize          = True  
  resh@gsnAddCyclic         = False
  
  resh@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels
  resh@cnMinLevelValF       = -.5            ; set min contour level
  resh@cnMaxLevelValF       =  .5             ; set max contour level
  resh@cnLevelSpacingF      =  0.05            ; set contour spacing   
  resh@cnFillOn             = True         ; turn on color fill
  resh@cnFillPalette        = "cmp_b2r"  
  resh@cnLevels             =  4.
 
  resh@cnFillOn             = True
  
  resh@cnSmoothingOn        = True
  resh@cnSmoothingDistanceF = 0.005
  resh@cnSmoothingTensionF  = 0.001
  
  resh@tfDoNDCOverlay       = True            ; necessary for correct overlay on map

  resh@cnLinesOn             = False            ; turn off contour lines
  resh@cnLineLabelsOn        = False            ; turn off contour labels
  
  resh@cnFillMode            = "AreaFill"
  resh@cnFillOpacityF        = 1.

  resh@cnInfoLabelOn        = False
  resh@lbLabelBarOn         = True       ; turn off individual lb's 
 
  resh@cnLineLabelInterval = 1           ; Default is 2
  
 ; resh@pmTickMarkDisplayMode = "Always"         ; turn on nicer tickmarks
  resh@tmXBLabelFontHeightF  = 0.040
  resh@tmYLLabelFontHeightF  = 0.040
 ; resh@tmYLLabelStride       = 2                ; label every other tickmark
  resh@tmXBLabelStride       = 2
  
  resh@tmXTOn                  = True 
  resh@tmYROn                  = True 
  ;resh@tmYRLalabelsbelsOn            = False
  resh@tmXTLabelsOn            = False
 
  resh@lbOrientation        = "vertical"
  resh@lbTitleFontHeightF  = 0.040 
  resh@lbLabelFontHeightF  = 0.040
  resh@lbLeftMarginF       = 0.10
  resh@lbRightMarginF      = 0.55
  ;resh@lbLabelOffsetF	  = 0.20
  resh@lbLabelFont   	  = "Helvetica-Bold"     ; label font
  resh@lbTitleString       = "SSH (m)"
  ;resh@lbLabelPosition	  = "Top"
  resh@lbTitleAngleF       = 270.
  resh@lbTitlePosition      = "Right"
  resh@lbTitleDirection     = "Across" 
  
  resh@lbLabelStride        = 10 
  
  resh@tmXBMode   = "Explicit"
  resh@tmXBValues = ispan(-50,10,10)
  resh@tmXBLabels = (/"50~S~o~N~W","40~S~o~N~W","30~S~o~N~W","20~S~o~N~W","10~S~o~N~W","0~S~o~N~","10~S~o~N~E","20~S~o~N~E"/)
  ;resh@tmXBMinorValues  = ispan(-50,20,5)
  resh@tmXBLabelFontHeightF = 0.030

  resh@tmYLMode = "Explicit"
  resh@tmYLValues = (/"0", "1","2","3","4","5","6","7","8","9","10","11"/)
  resh@tmYLLabels = (/"J", "F","M","A","M","J","J","A","S","O","N","D"/)
  ;resh@tmXBLabels = (/"-5", "-4","-3","-2","-1","0","+1","+2","+3","+4","+5","+6"/)
  resh@tiYAxisString = "Month"                  ; yaxis
  resh@tmYLLabelFontHeightF = 0.040
  
  resh@tiXAxisFontHeightF = 0.040
  resh@tiYAxisFontHeightF = 0.040
  
  resh@tmXTOn = True ; Turn off top tickmarks and labels
  resh@tmYROn = True ; Turn off left tickmarks and labels
  
  resh@tiXAxisFontHeightF = 0.040
  resh@tiYAxisFontHeightF = 0.040
  
  resh@vpHeightF 	       = 0.80        ; Changes the aspect ratio
  resh@vpWidthF 	       = 0.80
  
  ;---For higher quality large images
  resh@tiYAxisFontThicknessF  = 1.
  resh@tmXBMajorThicknessF    = 1.
  resh@tmXBMinorThicknessF    = 1.
  resh@tmYLMajorThicknessF    = 1.
  resh@tmYLMinorThicknessF    = 1.
  resh@tmYLLabelFontThicknessF= 1.
  resh@tmBorderThicknessF     = 1.
;-----------------------------------------------------------------------------
  ress1                      = True         
  ress1@gsnDraw              = False        ; don't draw yet
  ress1@gsnFrame             = False        ; don't advance frame yet
  ress1@gsnMaximize          = True  
  ress1@gsnAddCyclic         = False
  
  ress1@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels
  ress1@cnMinLevelValF       = -3.0            ; set min contour level
  ress1@cnMaxLevelValF       =  3.0             ; set max contour level
  ress1@cnLevelSpacingF      =  0.3            ; set contour spacing   
  ress1@cnFillOn             = True         ; turn on color fill
  ress1@cnFillPalette        = "cmp_b2r"  
;  ress1@cnLevels             =  4.
  
  ress1@cnFillOn             = True
  
  ress1@cnSmoothingOn        = True
  ress1@cnSmoothingDistanceF = 0.005
  ress1@cnSmoothingTensionF  = 0.001
  
  ress1@tfDoNDCOverlay        = True            ; necessary for correct overlay on map

  ress1@cnLinesOn             = False            ; turn off contour lines
  ress1@cnLineLabelsOn        = False            ; turn off contour labels
  
 ; ress1@cnFillMode            = "AreaFill"
 ; ress1@cnFillOpacityF        = 1.

  ress1@cnInfoLabelOn        = False
  ress1@lbLabelBarOn         = True       ; turn off individual lb's 
 
 ; ress1@cnLineLabelInterval = 1           ; Default is 2
  
 ; ress1@pmTickMarkDisplayMode = "Always"         ; turn on nicer tickmarks
  ress1@tmXBLabelFontHeightF  = 0.040
  ress1@tmYLLabelFontHeightF  = 0.040
 ; ress1@tmYLLabelStride       = 2                ; label every other tickmark
  ress1@tmXBLabelStride       = 2
 
  ress1@tmXTOn                  = True 
  ress1@tmYROn                  = True 
  ;ress1@tmYRLalabelsbelsOn            = False
  ress1@tmXTLabelsOn            = False
 
  ress1@lbOrientation        = "vertical"
  ress1@lbTitleFontHeightF  = 0.040 
  ress1@lbLabelFontHeightF  = 0.040
  ress1@lbLeftMarginF       = 0.10
  ress1@lbRightMarginF      = 0.55
  ;ress1@lbLabelOffsetF	  = 0.20
  ress1@lbLabelFont   	  = "Helvetica-Bold"     ; label font
  ress1@lbTitleString       = "Precipitation Anomalies (mm day~S~-1~N~)"
  ;ress1@lbLabelPosition	  = "Top"
  ress1@lbTitleAngleF       = 270.
  ress1@lbTitlePosition      = "Right"
  ress1@lbTitleDirection     = "Across"
  
  ress1@lbLabelStride        = 10 
  
  ress1@tmXBMode   = "Explicit"
  ress1@tmXBValues = ispan(-50,10,10)
  ress1@tmXBLabels = (/"50~S~o~N~W","40~S~o~N~W","30~S~o~N~W","20~S~o~N~W","10~S~o~N~W","0~S~o~N~","10~S~o~N~E","20~S~o~N~E"/)
  ;ress1@tmXBMinorValues  = ispan(-50,20,5)
  ress1@tmXBLabelFontHeightF = 0.040

  ress1@tmYLMode = "Explicit"
  ress1@tmYLValues = (/"0", "1","2","3","4","5","6","7","8","9","10","11"/)
  ress1@tmYLLabels = (/"J", "F","M","A","M","J","J","A","S","O","N","D"/)
  ;ress1@tmXBLabels = (/"-5", "-4","-3","-2","-1","0","+1","+2","+3","+4","+5","+6"/)
  ress1@tiYAxisString = "Month"                  ; yaxis
  ress1@tmYLLabelFontHeightF = 0.040
  
 
  ress1@tmXTOn = True ; Turn off top tickmarks and labels
  ress1@tmYROn = True ; Turn off left tickmarks and labels
  
  ress1@tiXAxisFontHeightF = 0.040
  ress1@tiYAxisFontHeightF = 0.040
 
  ress1@tiXAxisFontHeightF = 0.040
  ress1@tiYAxisFontHeightF = 0.040
  ress1@gsnLeftStringFontHeightF = 0.040
  
  ress1@vpHeightF 	  = 0.80        ; Changes the aspect ratio
  ress1@vpWidthF 	  = 0.80
 
  ;---For higher quality large images
  ress1@tiYAxisFontThicknessF  = 1.
  ress1@tmXBMajorThicknessF    = 1.
  ress1@tmXBMinorThicknessF    = 1.
  ress1@tmYLMajorThicknessF    = 1.
  ress1@tmYLMinorThicknessF    = 1.
  ress1@tmYLLabelFontThicknessF= 1.
  ress1@tmBorderThicknessF     = 1.
;-----------------------------------------------------------------------------
;===============================================================
  res1                      =   True
  res1@gsnDraw              = False        ; don't draw yet
  res1@gsnFrame             = False        ; don't advance frame yet
 ; res1@mpFillDrawOrder      = "preDraw"
  ;res1@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels
 ; res1@cnMinLevelValF       =  -1.0            ; set min contour level
 ; res1@cnMaxLevelValF       =  1.0             ; set max contour level
 ; res1@cnLevelSpacingF      =  .05            ; set contour spacing  
  res1@gsnAddCyclic         = False
  res1@cnLinesOn            = False       ; True is default
  res1@cnLineLabelsOn       = False        ; True is default
  res1@cnInfoLabelOn        = False
  res1@gsnRightString   = ""
  res1@gsnLeftString   = ""
  res1@cnLineLabelDensityF =  5
;-----------------------------------------------------------------------------
  opt = True
  opt@gsnShadeFillType = "pattern"        ; color is the default
  opt@gsnShadeLow = 17                     ; use fill pattern #17
  opt@gsnShadeHigh = 17                     ; use fill pattern #17
  opt@gsnShadeFillScaleF = 1.0
 
  
   delete(res@cnFillPalette)
   delete(ress@cnFillPalette)
   delete(ress1@cnFillPalette)
   delete(resh@cnFillPalette)
   res@cnFillPalette        	= "amwg_blueyellowred"       ;different colour bar
   ress@cnFillPalette        	= "amwg_blueyellowred"       ;different colour bar
   ress1@cnFillPalette        	= "amwg_blueyellowred"       ;different colour bar
   resh@cnFillPalette        	= "amwg_blueyellowred"       ;different colour bar
  ; res@gsnStringFontHeightF = 0.013
;=============================================================================
;===============================================================
   r3               = True                        ; plot mods desired
   r3@gsnDraw       = False                       ; do not draw
   r3@gsnFrame      = False                       ; do not advance frame
   r3@trYReverse            = True
 
   r3@xyLineThicknesses   = (/1.,1.,1.,1.,1.,1.,1.,1.,1./)       ; Define line thicknesses 
   r3@xyDashPatterns      = (/0.,0.,0.,0.,0.,0.,0.,0.,0./)           ; Dash patterns 
   r3@tmYLLabelFontColor = "foreground"   ; and right axes are same.
   r3@tiYAxisFontColor  = "foreground"
  
   plotA41  = gsn_csm_xy(wks,iso2&lon,iso2(tlev|0,lon|:),r3) ; create plot
;===============================================================
 
   ress@gsnLeftStringFontHeightF   = 0.050
   ress@gsnLeftString  		= "Canonical events"   
   plotA0 = gsn_csm_contour(wks,d_SSTP_C(month|:,{lon|-52:11}),ress)  ;           ; 
   plotS0 = gsn_csm_contour(wks,d_SSTdP_C(:,:),res1)       ;
   plotS0 = gsn_contour_shade(plotS0,p95,d_SSTdP_C@_FillValue,opt)
   overlay(plotA0,plotS0)
   plot(0) = plotA0
  
   ress@gsnLeftString  = "Non-canonical events"  
   plotA1 = gsn_csm_contour(wks,d_SSTP_N(month|:,{lon|-52:11}),ress)  ;           ; 
   plotS1 = gsn_csm_contour(wks,d_SSTdP_N(:,:),res1)       ;
   plotS1 = gsn_contour_shade(plotS1,p95,d_SSTdP_N@_FillValue,opt)
   overlay(plotA1,plotS1)
   plot(4) = plotA1
     
   ress1@gsnLeftStringFontHeightF   = 0.050
   ress1@gsnLeftString  		= "Canonical events"
   plotA2 = gsn_csm_contour(wks,rain_C(month|:,{lon|-52:13}),ress1)  ;  
   plotS2 = gsn_csm_contour(wks,raindP_C(:,:),res1)       ;
   plotS2 = gsn_contour_shade(plotS2,p95,raindP_C@_FillValue,opt)
   overlay(plotA2,plotS2)
   plot(1) = plotA2
   
   ress1@gsnLeftStringFontHeightF   = 0.050
   ress1@gsnLeftString  = "Non-canonical events"  
   plotA3 = gsn_csm_contour(wks,rain_N(month|:,{lon|-52:13}),ress1)  ;           ; 
   plotS3 = gsn_csm_contour(wks,raindP_N(:,:),res1)       ;
   plotS3 = gsn_contour_shade(plotS3,p95,raindP_N@_FillValue,opt)
   overlay(plotA3,plotS3)
   plot(5) = plotA3
   
   resh@gsnLeftStringFontHeightF   = 0.050
   resh@gsnLeftString  		= "Canonical events"
   plotA4 = gsn_csm_contour(wks,d_sshP_C(month|:,{lon|-52:11}),resh)  ;           ; 
   plotS4 = gsn_csm_contour(wks,d_sshdP_C(:,:),res1)       ;
   plotS4 = gsn_contour_shade(plotS4,p95,d_sshdP_C@_FillValue,opt)
   overlay(plotA4,plotS4)
   plot(2) = plotA4
   
   resh@gsnLeftStringFontHeightF   = 0.050
   resh@gsnLeftString  = "Non-canonical events"  
   plotA5 = gsn_csm_contour(wks,d_sshP_N(month|:,{lon|-52:11}),resh)  ;           ; 
   plotS5 = gsn_csm_contour(wks,d_sshdP_N(:,:),res1)       ;
   plotS5 = gsn_contour_shade(plotS5,p95,d_sshdP_N@_FillValue,opt)
   overlay(plotA5,plotS5)
   plot(6) = plotA5
   
   res@gsnLeftStringFontHeightF   = 0.050
   res@gsnLeftString  		= "Canonical events"
   plotA4 = gsn_csm_contour(wks,sub_C({level|0:300},{lon|-52:11}),res)  ;  
   plotS4 = gsn_csm_contour(wks,sub_SSTdP_C({level|0:300},{lon|-52:11}),res1)       ;
   plotS4 = gsn_contour_shade(plotS4,p95,sub_SSTdP_C@_FillValue,opt)
   overlay(plotA4,plotS4)
   overlay(plotA4,plotA41)
   plot(3) = plotA4
   
   res@gsnLeftStringFontHeightF   = 0.050
   res@gsnLeftString  = "Non-canonical events"  
   plotA51  = gsn_csm_xy(wks,iso2&lon,iso2(tlev|0,lon|:),r3) ; create plot
   plotA5 = gsn_csm_contour(wks,sub_N({level|0:300},{lon|-52:11}),res)  ;           ; 
   plotS5 = gsn_csm_contour(wks,sub_SSTdP_N({level|0:300},{lon|-52:11}),res1)       ;
   plotS5 = gsn_contour_shade(plotS5,p95,sub_SSTdP_N@_FillValue,opt)
   overlay(plotA5,plotS5)
   overlay(plotA5,plotA51)
   plot(7) = plotA5
;=============================================================================   
 delete(res@cnMinLevelValF)
 delete(res@cnMaxLevelValF)
 delete(res@cnLevelSpacingF)
;=====================================================================
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
  
  dumL0 = gsn_add_polyline(wks,plot(0), (/x_2,x_3/),(/y_2,y_3/) ,pres1)
  dumL1 = gsn_add_polyline(wks,plot(1), (/x_2,x_3/),(/y_2,y_3/) ,pres1)
  
  dumL2 = gsn_add_polyline(wks,plot(2), (/x_2,x_3/),(/y_2,y_3/) ,pres1)
  dumL3 = gsn_add_polyline(wks,plot(3), (/x_2,x_3/),(/y_2,y_3/) ,pres1)
  
  dumL4 = gsn_add_polyline(wks,plot(4), (/x_2,x_3/),(/y_2,y_3/) ,pres1)
  dumL5 = gsn_add_polyline(wks,plot(5), (/x_2,x_3/),(/y_2,y_3/) ,pres1)
  
  dumL6 = gsn_add_polyline(wks,plot(6), (/x_2,x_3/),(/y_2,y_3/) ,pres1)
 ; dumL7 = gsn_add_polyline(wks,plot(7), (/x_2,x_3/),(/y_2,y_3/) ,pres1)
;========================================================

  txres5               = True
  txres5@txAngleF      = 0.          ; Rotate string clockwise 
  txres5@txFontHeightF = 0.025   ; Use same font height as left axis
  txres5@txFont                = "Helvetica-Bold"
 ; txid = gsn_create_text(wks, "DT anomaly (~S~o~N~C)", txres)
   txid5 = gsn_create_text(wks, "a", txres5);
; Move text string to center/right edge of plot.;
  amres5                  = True
  amres5@amParallelPosF   =  -0.6 ; 0.5 is the right edge of the plot, so
                                 ; 0.6 is a little further to the right.
  amres5@amOrthogonalPosF =  -0.6   ; This is the center of the plot.
 ; tt_p0X = gsn_add_annotation(plot(0),txid5,amres5)
  
  txres5@txFontHeightF = 0.038   ; Use same font height as left axis
  amres5@amParallelPosF   =  -0.68 ; 0.5 is the right edge of the plot, so
  amres5@amOrthogonalPosF =  -0.55    ; This is the center of the plot.
  amres5@amJust           = "TopLeft"

 tt = (/"a","b","c","d","e","f","g","i"/)
create_tt = gsn_create_text(wks,tt,txres5);
; tt_p0 = new(nplots,graphic)                ; create graphic array ; only needed if paneling
tt_p0 = new(8,graphic)                ; create graphic array ; only needed if paneling
;---Attach the polyline
     do i = 0,7
	print(i)
        tt_p0(i) = gsn_add_annotation(plot(i),create_tt(i),amres5)
     end do
delete(i)
;==============================================================================================
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
  resP@gsnLeftString 			= ""
  ;resP@txString          		= "SST trends: 1984-2008"
  resP@txString          		= ""

;**************************************************************
   gsn_panel(wks,plot(0:),(/2,4/),resP)     ; now draw as one plot
;=====================================================================
;========================================================================
frame(wks)
end
