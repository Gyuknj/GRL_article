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
  printVarSummary(hadi)       ; variable overview
  printMinMax(hadi, True)
	
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
  printVarSummary(hadi)
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

  printVarSummary(hadisst_d)
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
;============================================================== ============================================  
  hadi_s = dtrend_msg_n(hadi_s&time,hadi_s,True,False,0)        ;removes linear trends
  printVarSummary(hadi_s)
;============================================================== ============================================ 
;------------------------------------------------------------------------------------------------------------------------------------
;------------------------------------------------------------------------------------------------------------------------------------
  ahadi_J = wgt_areaave_Wrap(hadi_s(time|:,month|5,{lon|-20:0}, {lat|-3:3}),1.0, 1.0, 0)        ;atl3 index
  printVarSummary(ahadi_J)
  printMinMax(ahadi_J,True)  
;============================================================== ============================================  
 ; ahadi_J = dim_standardize_n_Wrap(ahadi_J,1,0)      ; normalize 
 ; printVarSummary(ahadi_J)
 ; printMinMax(ahadi_J,True)  

;=================================================================================================================  
;=================================================================================================================================
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
;==========================================================================================
;==============================================================
  arain_M = wgt_areaave_Wrap(precip_s(time|:,month|4,{lon|-40:5}, {lat|-3:3}),1.0, 1.0, 0)        ;atl3 index
  printVarSummary(arain_M)
  printMinMax(arain_M,True)  
;============================================================== ============================================  
 ; arain_M = dim_standardize_n_Wrap(arain_M,1,0)      ; normalize 
 ; printVarSummary(arain_M)
 ; printMinMax(arain_M,True)  

;==============================================================
  ;Correlations
;==============================================================  
  r_sst_rain = escorc_n(ahadi_J,arain_M,0,0)
  print(r_sst_rain)
;===============================================================================================
  wks = gsn_open_wks("eps","Figure_S2")         ; send graphics to PNG file
  plot = new(1,graphic)                ; create graphic array ; only needed if paneling

   r0               = True                        ; plot mods desired
   r0@gsnDraw       = False                       ; do not draw
   r0@gsnFrame      = False                       ; do not advance frame
   r0@tiMainString  = ""                       ; title
   r0@tiMainFontHeightF  = 0.045              
   r0@tiXAxisString = ""                  ; yaxis
   r0@tiXAxisFontHeightF = 0.030
   r0@tiYAxisFontHeightF = 0.030
   r0@tiYAxisString = "CHIPRS"                  ; yaxis
   r0@tiXAxisString = "Precipitation Anomalies (mm day~S~-1~N~)"                  ; xaxis
   r0@tmXBLabelFontHeightF = 0.030
   r0@tmYLLabelFontHeightF = 0.030
   r0@trYMinF             =  -1.2                 ; manually set lower limit
   r0@trYMaxF             =   1.2           
   r0@trXMinF             =  -4.0                 ; manually set lower limit
   r0@trXMaxF             =   4.0               ;   "          upper
   r0@tmXTBorderOn            =  True
   r0@tmXTOn                  =  True
   r0@tmYLFormat	      = "@*+^sg"	; use "e" symbol, don't force 0's (see documentation)
   r0@tmXBFormat              = "@*+^sg"        ;
   r0@gsnCenterStringFontHeightF  = 0.020 
   r0@gsnCenterStringOrthogonalPosF = -0.20
   r0@gsnCenterStringFontColor 	    = "red"
   r0@tmXMajorGrid                = False          ; implement x grid 
   r0@tmXMajorGridThicknessF      = 1.0           ; 2.0 is default
   r0@tmXMajorGridLineDashPattern = 2             ; select short dash lines
   r0@tmXMajorGridLineColor	 = "black"
   r0@tmYMajorGrid                = False          ; implement y grid
   r0@tmYMajorGridThicknessF      = 1.0           ; 2.0 is default
   r0@tmYMajorGridLineDashPattern = 2             ; select short dash lines
   r0@tmYMajorGridLineColor	 = "black"
   r0@tiYAxisFontThicknessF  = 1.
   r0@tmXBMajorThicknessF    = 1.
   r0@tmXBMinorThicknessF    = 1.
   r0@tmYLMajorThicknessF    = 1.
   r0@tmYLMinorThicknessF    = 1.
   r0@tmYLLabelFontThicknessF= 1.
   r0@tmBorderThicknessF     = 1.
   r0@tmXBMajorOutwardLengthF = 0.0               ; draw tickmarks inward
   r0@tmXBMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
   r0@tmYLMajorOutwardLengthF = 0.0               ; draw tickmarks inward
   r0@tmYLMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
   r0@tmXTMajorOutwardLengthF = 0.0               ; draw tickmarks inward
   r0@tmXTMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
   r0@tmYRMajorOutwardLengthF = 0.0               ; draw tickmarks inward
   r0@tmYRMinorOutwardLengthF = 0.0               ; draw minor ticsk inward
   r0@xyMarkLineMode	      = "Markers"
   r0@xyMarkerColor	      = "blue"
   r0@xyMarker		      =  16
;=======================================================================================================================
  r0@gsnCenterStringFontHeightF = 0.025

; create individual plots
   r0@tiYAxisString = "SSTA (~S~o~N~C)"                  ; yaxis
   r0@tiMainString  = ""
   r0@gsnCenterString = "r ="+ sprintf("%5.2f", r_sst_rain)
        plot(0) = gsn_csm_xy (wks,arain_M(:),ahadi_J(:),r0)
;========================================================================
;===============================================================
;===============================================================
  x0 =  -4.0
  x1 =  4.0
  y0 =  -1.2
  y1 =  1.2
;---Resources for the polyline
  pres1                  = True           ; polyline mods desired
  pres1@gsLineThicknessF = 1.5            ; line thickness   
  pres1@gsLineColor      = "black"         ; color of lines
  pres1@gsLineDashPattern = 2
  pres1@tfPolyDrawOrder   = "postDraw"
  pres1@gsLineColor      = "black"         ; color of lines
  
   dumA0 = gsn_add_polyline(wks,plot(0), (/x0,x1/),(/y0,y1/) ,pres1)
;---------------------------------------------------------------------
;=====================================================================
;=====================================================================
;=====================================================================
  resT                                  = True         ; modify the panel plot
  resT@gsnDraw                          = True      ; don't draw yet
  resT@gsnFrame                         = False        ; don't advance frame yet
  resT@gsnPanelLabelBar                 = False        ; add common colorbar
  resT@gsnPanelYWhiteSpacePercent       = 0.0    ; default is 1.0
  resT@gsnPanelXWhiteSpacePercent       = 1.0    ; default is 1.0
  resT@gsnPanelTop                      = 0.99
  resT@gsnPanelBottom                   = 0.01
  resT@gsnPanelLeft                     = 0.01
  resT@gsnPanelRight                    = 0.99
  resT@gsnMaximize                      = True         ; large format
  resT@gsnPaperOrientation              = "portrait"   ; force portrait
  resT@gsnLeftString                    = ""
  ;resT@txString                        = "SST trends: 1984-2008"
  ;resT@txString                        = "JAS SST trends: 1979-2018"
;**************************************************************
;=====================================================================
  gsn_panel(wks,plot(0:),(/1,1/),resT)               ; create panel plots
;*******************************************
;===============================================================  

  frame(wks)
end  

