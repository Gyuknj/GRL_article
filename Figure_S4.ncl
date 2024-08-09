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
  
  p95 = 0.10
 
  diri = "/home/gyuk/Downloads/GYUK/GEOMAR/Final_work/"
  
  sst_data = "HadISST_sst.nc"
;==============================================================
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
  ;printVarSummary(hadi_s)


  hadi_s = dtrend_msg_n(hadi_s&time,hadi_s,True,False,0)        ;removes linear trends
  ;printVarSummary(hadi_s)
  ;printMinMax(hadi_s,True)

;---------------------------------------------------------------------------
  ahadi_s = wgt_areaave_Wrap(hadi_s(time|:,month|:,{lon|-20:0}, {lat|-3:3}),1.0, 1.0, 0)        ;atl3 index
  printVarSummary(ahadi_s)
  ;printMinMax(ahadi_s,True)
;----------------------------
  delete(hadi)
;  delete(hadi_s)

  jja_ahadi = dim_avg_n_Wrap(ahadi_s(time|:,month|5:7),1)		;JJA, ncl counts from zero
  print(jja_ahadi)

   jja_ahadi = dim_standardize_n(jja_ahadi, 0, 0)  ;Normalise deviations
   print(jja_ahadi)
;
;----------------------------------------------------------------------------
; Positive cases
;============================================================================
  n_atlP = dim_num_n(jja_ahadi.ge.thres,0)        ;Number of cases
  print(n_atlP)

  valsP = ind(jja_ahadi.ge.thres)
  print(valsP)
  
  ;print("=========================================")
  print(valsP+yrStrt)		;the years
;----------------------------------------------------------------------------
; Negative cases,
; note the "-" sign attached to the threshold here
;============================================================================
  n_atlN = dim_num_n(jja_ahadi.le.-thres,0)        ;Number of cases
  print(n_atlN)
;
  valsN = ind(jja_ahadi.le.-thres)
  ;print(valsN)
 ; print("=========================================")
  print(valsN+yrStrt)             ;the years
;==============================================================
;          Difference in Canonical Events  
;==============================================================
  dimZ = dimsizes(hadi_s)
  printVarSummary(dimZ)
      
        saodP_cases_maps  = new((/n_atlP,dimZ(1),dimZ(2),dimZ(3)/),float,hadi_s@_FillValue)
        printVarSummary(saodP_cases_maps)
    
        saodP_cases_maps = hadi_s(valsP,:,:,:)
        printVarSummary(saodP_cases_maps)        
;==============================================================       
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
;===================================================================
;==============================================================
;==============================================================

        saodN_cases_maps  = new((/n_atlN,dimZ(1),dimZ(2),dimZ(3)/),float,hadi_s@_FillValue)
        printVarSummary(saodN_cases_maps)
    
        saodN_cases_maps = hadi_s(valsN,:,:,:)
        printVarSummary(saodN_cases_maps)
        
        
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
;==============================================================

  pc1SST111 = (dim_avg_n_Wrap(sst_fp1,0) - dim_avg_n_Wrap(sst_fn1,0))/(2.)
  copy_VarCoords(hadi_s(0,:,:,:),pc1SST111)                             ; copy coordinate arrays 
  printVarSummary(pc1SST111)
  printMinMax(pc1SST111,True)
;===================================== 
  pc1SSTd111 = pc1SST111  ; Transfer metdata for plotting sig. values        

  pc1SSTd111 = ttest((dim_avg_n_Wrap(sst_fp1,0)),(dim_stddev_n_Wrap(sst_fp1,0))^2,5\
               ,(dim_avg_n_Wrap(sst_fn1,0)),(dim_stddev_n_Wrap(sst_fn1,0))^2,5,True, False)
  pc1SSTd111@long_name = "Probability: difference between means"
  printVarSummary(pc1SSTd111)

  pc1SSTd111 = where(abs(pc1SSTd111).le.p95,pc1SSTd111,pc1SSTd111@_FillValue) ; Sig.u winds, where UWND or VWND is significanlt
  ;print("==============================================================")  
  printVarSummary(pc1SSTd111)
  printMinMax(pc1SSTd111,True)
  
  d_hadiP61110 = (pc1SST111(month|0,lat|:, lon|:))  
  d_hadidP11110 =(pc1SSTd111(month|0,lat|:, lon|:))
  
  d_hadiP61111 = (pc1SST111(month|1,lat|:, lon|:))  
  d_hadidP11111 =(pc1SSTd111(month|1,lat|:, lon|:))
  
  d_hadiP61112 = (pc1SST111(month|2,lat|:, lon|:))  
  d_hadidP11112 = (pc1SSTd111(month|2,lat|:, lon|:))
  
  d_hadiP61113 = (pc1SST111(month|3,lat|:, lon|:))  
  d_hadidP11113 = (pc1SSTd111(month|3,lat|:, lon|:))
  
  d_hadiP61114 = (pc1SST111(month|4,lat|:, lon|:))  
  d_hadidP11114 = (pc1SSTd111(month|4,lat|:, lon|:))
  
  d_hadiP61115 = (pc1SST111(month|5,lat|:, lon|:))  
  d_hadidP11115 = (pc1SSTd111(month|5,lat|:, lon|:))
  
  d_hadiP61116 = (pc1SST111(month|6,lat|:, lon|:))  
  d_hadidP11116 =(pc1SSTd111(month|6,lat|:, lon|:))
  
  d_hadiP61117 = (pc1SST111(month|7,lat|:, lon|:))  
  d_hadidP11117 =(pc1SSTd111(month|7,lat|:, lon|:))
  
  d_hadiP61118 = (pc1SST111(month|8,lat|:, lon|:))  
  d_hadidP11118 = (pc1SSTd111(month|8,lat|:, lon|:))
  
  d_hadiP61119 = (pc1SST111(month|9,lat|:, lon|:))  
  d_hadidP11119 = (pc1SSTd111(month|9,lat|:, lon|:))
  
  d_hadiP611110 = (pc1SST111(month|10,lat|:, lon|:))  
  d_hadidP111110 = (pc1SSTd111(month|10,lat|:, lon|:))
  
  d_hadiP611111 = (pc1SST111(month|11,lat|:, lon|:))  
  d_hadidP111111 = (pc1SSTd111(month|11,lat|:, lon|:)) 
;=================================================================== 
;==============================================================
;==============================================================  
  wks = gsn_open_wks("eps","Figure_S4")         ; send graphics to PNG file
  plot = new(12,graphic)                ; create graphic array ; only needed if paneling
  ; plotS = new(nplots,graphic)                ; 
  

  res                      = True         
  res@gsnDraw              = False        ; don't draw yet
  res@gsnFrame             = False        ; don't advance frame yet

  res@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels
  res@cnMinLevelValF       = -1.            ; set min contour level
  res@cnMaxLevelValF       =  1.             ; set max contour level
  res@cnLevelSpacingF      =  0.1            ; set contour spacing  
 
  res@mpFillOn             = False        ; turn off map fill
  res@pmTickMarkDisplayMode = "Always"
  res@mpMinLatF            = latS         ; zoom in on map
  res@mpMaxLatF            = latN
  res@mpMinLonF            = lonL
  res@mpMaxLonF            = lonR
 
  ;res@tmXBLabelsOn          = False      ;no lon ticks
  ;res@tmYLMode   = "Explicit"
  ;res@tmYLValues = ispan(-10,10,10)
  ;res@tmYLLabels = (/"10~S~o~N~S", "0~S~o~N","10~S~o~N~N"/)
  ;res@tmYLMinorValues  = ispan(-15,20,5) 
  res@tmYLLabelFontHeightF = 0.015

 ; res@tmXBMode   = "Explicit"
  ;res@tmXBValues = ispan(-60,20,20)
  ;res@tmXBLabels = (/"60~S~o~N~W", "40~S~o~N~W","20~S~o~N~W","0~S~o~N", \
  ;                   "20~S~o~N~E"/)
  ;res@tmXBMinorValues  = ispan(-60,20,10)
  res@tmXBLabelFontHeightF = 0.015
      
;  res@mpGeophysicalLineColor  = "black" 
  res@mpFillDrawOrder      = "preDraw"

  res@gsnAddCyclic         = True
  res@cnFillOn             = True         ; turn on color fill
  res@cnLinesOn            = False        ; True is default
  res@cnLineLabelsOn       = False        ; True is default
  res@cnMissingValFillColor = "grey" ; comment this line out 
 
  res@mpOceanFillColor        = "white"
  res@mpGeophysicalLineColor  = "grey35" 
 
  res@cnFillPalette = "ncl_default"
  
  res@lbLabelBarOn         =  False       ; turn off individual lb's
  res@lbOrientation        = "vertical"
  res@lbTitleFontHeightF  = .013 
  res@lbLabelFontHeightF   = 0.013
  res@lbLeftMarginF         = 0.12
  res@lbRightMarginF        = 0.55
  res@lbLabelFont   		     = "Helvetica-Bold"     ; label font
 
  ;res@lbTitleString        = "SSTA (~S~o~N~C)"
  res@lbTitleAngleF        = 270.
  res@lbTitlePosition      = "Right"
  res@lbTitleDirection     = "Across"
  res@gsnLeftStringFontHeightF   = 0.025
  res@lbLabelStride        = 10

  ;res@vpHeightF 	  = 0.80        ; Changes the aspect ratio
  ;res@vpWidthF 		  = 0.40
 ; res@gsnCenterStringFontHeightF   = 0.025
 ; res@vpXF      = 0.30        ; change start locations
 ; res@vpYF      = 0.75        ; the plot
 
 
  res@tiYAxisFontThicknessF    = 1.
  res@tmXBMajorThicknessF      = 1.
  res@tmXBMinorThicknessF      = 1.
  res@tmYLMajorThicknessF      = 1.
  res@tmYLMinorThicknessF      = 1.
  res@tmYLLabelFontThicknessF  = 1.
  res@tmBorderThicknessF       = 1.  
;=============================================================================
;-----------------------------------------------------------------------------
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
  opt@gsnShadeLow = 17                     ; use fill pattern #4
  opt@gsnShadeHigh = 17                  ; use pattern #2   
  ;opt@gsnShadeMid = 2   
  
  delete(res@cnFillPalette)
  res@cnFillPalette        	= "amwg_blueyellowred"       ;different colour bar
  res@gsnCenterString  	= ""
      
   res@gsnLeftString  = "January"  
   plotA0 = gsn_csm_contour_map_ce(wks,d_hadiP61110(:,:),res)  ;           ; 
   plotS0 = gsn_csm_contour(wks,d_hadidP11110(:,:),res1)       ;
   plotS0 = gsn_contour_shade(plotS0,p95,d_hadidP11110@_FillValue,opt)
   overlay(plotA0,plotS0)
   plot(0) = plotA0
   
   res@gsnLeftString  = "February"  
   plotA1 = gsn_csm_contour_map_ce(wks,d_hadiP61111(:,:),res)  ;           ; 
   plotS1 = gsn_csm_contour(wks,d_hadidP11111(:,:),res1)       ;
   plotS1 = gsn_contour_shade(plotS1,p95,d_hadidP11111@_FillValue,opt)
   overlay(plotA1,plotS1)
   plot(1) = plotA1
   
   
   res@gsnLeftString  = "March"  
   plotA2 = gsn_csm_contour_map_ce(wks,d_hadiP61112(:,:),res)  ;           ; 
   plotS2 = gsn_csm_contour(wks,d_hadidP11112(:,:),res1)       ;
   plotS2 = gsn_contour_shade(plotS2,p95,d_hadidP11112@_FillValue,opt)
   overlay(plotA2,plotS2)
   plot(2) = plotA2
   
   
   res@gsnLeftString  = "April"  
   plotA3 = gsn_csm_contour_map_ce(wks,d_hadiP61113(:,:),res)  ;           ; 
   plotS3 = gsn_csm_contour(wks,d_hadidP11113(:,:),res1)       ;
   plotS3 = gsn_contour_shade(plotS3,p95,d_hadidP11113@_FillValue,opt)
   overlay(plotA3,plotS3)
   plot(3) = plotA3
   
   
   res@gsnLeftString  = "May"  
   plotA4 = gsn_csm_contour_map_ce(wks,d_hadiP61114(:,:),res)  ;           ; 
   plotS4 = gsn_csm_contour(wks,d_hadidP11114(:,:),res1)       ;
   plotS4 = gsn_contour_shade(plotS4,p95,d_hadidP11114@_FillValue,opt)
   overlay(plotA4,plotS4)
   plot(4) = plotA4
   
   
   res@gsnLeftString  = "June"  
   plotA5 = gsn_csm_contour_map_ce(wks,d_hadiP61115(:,:),res)  ;           ; 
   plotS5 = gsn_csm_contour(wks,d_hadidP11115(:,:),res1)       ;
   plotS5 = gsn_contour_shade(plotS5,p95,d_hadidP11115@_FillValue,opt)
   overlay(plotA5,plotS5)
   plot(5) = plotA5
   
   res@gsnLeftString  = "July"  
   plotA6 = gsn_csm_contour_map_ce(wks,d_hadiP61116(:,:),res)  ;           ; 
   plotS6 = gsn_csm_contour(wks,d_hadidP11116(:,:),res1)       ;
   plotS6 = gsn_contour_shade(plotS6,p95,d_hadidP11116@_FillValue,opt)
   overlay(plotA6,plotS6)
   plot(6) = plotA6
   
   res@gsnLeftString  = "August"  
   plotA7 = gsn_csm_contour_map_ce(wks,d_hadiP61117(:,:),res)  ;           ; 
   plotS7 = gsn_csm_contour(wks,d_hadidP11117(:,:),res1)       ;
   plotS7 = gsn_contour_shade(plotS7,p95,d_hadidP11117@_FillValue,opt)
   overlay(plotA7,plotS7)
   plot(7) = plotA7
   
   
   res@gsnLeftString  = "September"  
   plotA8 = gsn_csm_contour_map_ce(wks,d_hadiP61118(:,:),res)  ;           ; 
   plotS8 = gsn_csm_contour(wks,d_hadidP11118(:,:),res1)       ;
   plotS8 = gsn_contour_shade(plotS8,p95,d_hadidP11118@_FillValue,opt)
   overlay(plotA8,plotS8)
   plot(8) = plotA8
   
   
   res@gsnLeftString  = "October"  
   plotA9 = gsn_csm_contour_map_ce(wks,d_hadiP61119(:,:),res)  ;           ; 
   plotS9 = gsn_csm_contour(wks,d_hadidP11119(:,:),res1)       ;
   plotS9 = gsn_contour_shade(plotS9,p95,d_hadidP11119@_FillValue,opt)
   overlay(plotA9,plotS9)
   plot(9) = plotA9
   
   
   res@gsnLeftString  = "November"  
   plotA10 = gsn_csm_contour_map_ce(wks,d_hadiP611110(:,:),res)  ;           ; 
   plotS10 = gsn_csm_contour(wks,d_hadidP111110(:,:),res1)       ;
   plotS10 = gsn_contour_shade(plotS10,p95,d_hadidP111110@_FillValue,opt)
   overlay(plotA10,plotS10)
   plot(10) = plotA10
   
   
   res@gsnLeftString  = "December"  
   plotA11 = gsn_csm_contour_map_ce(wks,d_hadiP611111(:,:),res)  ;           ; 
   plotS11 = gsn_csm_contour(wks,d_hadidP111111(:,:),res1)       ;
   plotS11 = gsn_contour_shade(plotS11,p95,d_hadidP111111@_FillValue,opt)
   overlay(plotA11,plotS11)
   plot(11) = plotA11

 delete(res@cnMinLevelValF)
 delete(res@cnMaxLevelValF)
 delete(res@cnLevelSpacingF)
;===============================================================
;===============================================================
;===============================================================
;===============================================================
  resP                                  = True         ; modify the panel plot
  resP@gsnDraw                          = True      ; don't draw yet
  resP@gsnFrame                         = False        ; don't advance frame yet
  resP@gsnPanelLabelBar                 = True        ; add common colorbar
  resP@gsnPanelYWhiteSpacePercent       = 1.0    ; default is 1.0
  resP@gsnPanelXWhiteSpacePercent       = 1.0    ; default is 1.0
  resP@gsnPanelTop                      = 0.89
  resP@gsnPanelBottom                   = 0.59
  resP@gsnPanelLeft                     = 0.05
  resP@gsnPanelRight                    = 0.99
  resP@gsnMaximize                      = True         ; large format
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
  resP@gsnPanelFigureStrings            = (/"a","b","c","d","e","f","g","h","i","j","k","l"/) ; add strings to panel
  resP@amJust                           = "TopLeft"
  resP@gsnPanelFigureStringsFontHeightF =  0.007
  resP@gsnPanelFigureStringsFont        =  "Helvetica-Bold"
  
 ; resP@gsnPanelLeftString               = "Difference in Evolution of SSTA in Canonical Events"  
  resP@txFont                           = "Helvetica-Bold"
  resP@txString                        = "Evolution of SSTA in non-convective events"
  ;resP@txString                        = "JAS SST trends: 1979-2018"
;**************************************************************
   gsn_panel(wks,plot(0:),(/4,3/),resP)     ; now draw as one plot
;=====================================================================
;=====================================================================
;=====================================================================
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
  amres5@amParallelPosF   =  -0.48 ; 0.5 is the right edge of the plot, so
                                 ; 0.6 is a little further to the right.
  amres5@amOrthogonalPosF =  -0.48   ; This is the center of the plot.
 ; tt_p0X = gsn_add_annotation(plot(0),txid5,amres5)
  
  txres5@txFontHeightF = 0.035   ; Use same font height as left axis
  amres5@amParallelPosF   =  -0.57 ; 0.5 is the right edge of the plot, so
  amres5@amOrthogonalPosF =  -0.55    ; This is the center of the plot.
  amres5@amJust           = "TopLeft"

 tt = (/"a","b","c","d","e","f","g","h","i","j","k","l"/)
create_tt = gsn_create_text(wks,tt,txres5);
; tt_p0 = new(nplots,graphic)                ; create graphic array ; only needed if paneling
tt_p0 = new(1,graphic)                ; create graphic array ; only needed if paneling
;---Attach the polyline
;     do i = 0,4*nvar-1,1
;	print(i)
;        tt_p0(i) = gsn_add_annotation(plot(i),create_tt(i),amres5)
;     end do
;delete(i)
;
;========================================================================
  latEQ = (/ 5.,  5.,  -10.,   -10., 5./)
  lonEQ = (/-80., -50., -50., -80.,-80./)
;---Resources for the polyline
  pres0                  = True           ; polyline mods desired
  pres0@gsLineThicknessF = 2.            ; line thickness   
  pres0@gsLineColor      = "black"         ; color of lines
  pres0@gsLineDashPattern = 0
  pres0@tfPolyDrawOrder   = "PostDraw"
  

; dum0 = gsn_add_polyline(wks,plot(0),lonEQ,latEQ,pres0)

 latEQ1 = (/5.,  5.,  -10.,   -10., 5./)
 lonEQ1 = (/10., 30., 30., 10.,10./)
 ;dum1 = gsn_add_polyline(wks,plot(0),lonEQ1,latEQ1,pres0)
 
 latEQ2 = (/5.,  5.,  -10.,   -10., 5./)
 lonEQ2 = (/100., 150., 150., 100.,100./)
 ;dum2 = gsn_add_polyline(wks,plot(0),lonEQ2,latEQ2,pres0)

; dum0 = new(2*nvar,graphic)    
; dum1 = new(2*nvar,graphic)
; dum2 = new(2*nvar,graphic)

; do i = 0,2*nvar-1,1
;         dum0(i) =  gsn_add_polyline(wks,plot(i),lonEQ,latEQ,pres0)
;	 dum1(i) =  gsn_add_polyline(wks,plot(i),lonEQ1,latEQ1,pres0)
;	 dum2(i) =  gsn_add_polyline(wks,plot(i),lonEQ2,latEQ2,pres0)
; end do
; delete(i)
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
;=====================================================================
;=====================================================================
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
;=============================================================================
;=============================================================================
frame(wks)
end
