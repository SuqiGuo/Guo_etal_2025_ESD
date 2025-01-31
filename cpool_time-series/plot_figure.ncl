load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;....................
undef("fill_xy_plot_area")
procedure fill_xy_plot_area(wks,plot)
begin
  getvalues plot
    "trXMinF"   : xmin
    "trYMinF"   : ymin
    "trXMaxF"   : xmax
    "trYMaxF"   : ymax
  end getvalues

  xbox = (/xmin,xmax,xmax,xmin,xmin/)
  ybox = (/ymin,ymin,ymax,ymax,ymin/)

  gnres = True
  gnres@gsFillColor     = "white" ; "honeydew"
  gnres@tfPolyDrawOrder = "PreDraw"
  plot@$unique_string("box")$ = gsn_add_polygon(wks,plot,xbox,ybox,gnres)
end
;......................................

begin
myfile13 = addfile("/work/bm1147/b380949/web-monitoring/irri-crop/Eyr/cVeg/cVeg_irri-crop_mpiesm_signal-separated_revised.nc", "r")
ti4=cd_calendar(myfile13->time,-1)
myfile1 = addfile("/work/bm1147/b380949/paper1-figure-scripts/data/cropland_expansion_timeseries.nc", "r")
myfile2 = addfile("/work/bm1147/b380949/paper1-figure-scripts/data/Afforestation_timeseries.nc", "r")
myfile3 = addfile("/work/bm1147/b380949/paper1-figure-scripts/data/Irrigation_timeseries.nc", "r")
Gpool = myfile1->Gpool ;mpiesm
Gpool2 = myfile1->Gpool2 ;cesm
Gpool4 = myfile1->Gpool4 ;ecearth
printVarSummary(Gpool)
Gflux = myfile2->Gflux ;mpiesm
Gflux2 = myfile2->Gflux2 ;cesm
Gflux4 = myfile2->Gflux4 ;ecearth

Gpool3 = myfile3->Gflux ;mpiesm
Gpool32 = myfile3->Gflux2 ;cesm
Gpool43 = myfile3->Gflux4 ;ecearth



;.......plot...............
;yrfrac4  = yyyymm_to_yyyyfrac (ti4, 0)
yrfrac4  = ispan(0,159,1)
wks_type = "pdf"
wks_type@wkPaperWidth =  7.99   
wks_type@wkPaperHeight =8.2
;wks_type@wkPaperWidthF =  7.99   
;wks_type@wkPaperHeightF =8.2
wks = gsn_open_wks(wks_type,"cpool_final_16Oct")

lres=True
lres@gsnDraw         = False                                                          ;         do not draw the plot                
lres@gsnFrame        = False     
;lres@vpHeightF= 0.2                    ; change aspect ratio of plot
;lres@vpWidthF = 0.6     
lres@gsnPaperMargin = 0
lres@mpShapeMode  = "FreeAspect"
lres@vpHeightF= 0.3                    ; change aspect ratio of plot
lres@vpWidthF = 1              
lres@xyLineThicknessF = 2.8
lres@xyMonoDashPattern= True
lres@xyDashPattern=0
;lres@tmXBLabelDeltaF= 1.0

lres@tmXTOn        = False ; turn off top tickmarks
lres@tmYROn        = False ; turn off right tickmarks
lres@tmXBOn        = False
lres@tmXTBorderOn  = True ; turn off top axis
lres@tmYRBorderOn  = True ; turn off right axis
lres@tmBorderThicknessF= 1.0
     lres@tmXBLabelFontHeightF = 0.018
     lres@tmXTLabelFontHeightF = 0.018
     lres@tmYLLabelFontHeightF = 0.018
     lres@tiMainFontHeightF = 0.018
lres@tiXAxisFontHeightF= 0.019
lres@tiYAxisFontHeightF= 0.019
     lres@tiMainFont = "times-roman"    ;"Helvetica-Bold"
     lres@tmXBLabelFont="times-roman"
     lres@tmYLLabelFont="times-roman"
     lres@tiXAxisFont ="times-roman"
     lres@tiYAxisFont="times-roman"

lres@gsnLeftStringFontHeightF=0.018
lres@gsnLeftStringOrthogonalPosF=-0.18 ;0.025
lres@gsnLeftStringParallelPosF= 0.01 ;-0.15
lres@tiMainOffsetYF= 0.8 ;0.025  ;up and down
lres@tiMainOffsetXF=0.1  ;left and right
;lres@tiMainOffsetYF = -0.15
lres@tmXBMajorOutwardLengthF = 0.0
lres@tmXBMinorOutwardLengthF = 0.0
lres@tmYLMajorOutwardLengthF=0.0
lres@tmXBMajorLengthF=0.01
lres@tmXBMinorLengthF=0.007
lres@tmYLMajorLengthF=0.01

lres@gsnStringFont="times-Bold"
lres@tmYLMinorOn = False
;lres@tmGridDrawOrder= "PreDraw"
; lres@tmXMajorGrid          = True               ; Add white grid lines
; lres@tmXMajorGridLineColor = "grey77"
; lres@tmYMajorGrid          = True
; lres@tmYMajorGridLineColor = "grey77"

lres@trYMinF=-30
lres@trYMaxF=45
lres@trXMinF  =  0     ; 2015        ; min value on y-axis
lres@trXMaxF  =  159   ; 2177            ; max value on y-axis

lres@tmXBMode= "Manual"
lres@tmYLMode= "Manual"

lres@tmXBTickStartF = 0
lres@tmXBTickSpacingF = 30
lres@tmXBTickEndF = 160
lres@tmXBMinorPerMajor = 2 

lres@tmYLTickStartF = -30
lres@tmYLTickSpacingF = 15
lres@tmYLTickEndF = 45

lres@tiYAxisOffsetXF= -0.01


plot=new(3,graphic)
plot2=new(3,graphic)
plot3=new(3,graphic)

plotn=new(3,graphic)
lres@xyLineColors= (/"palegreen4","sandybrown","royalblue","black"/) 

lres@tmXBLabelsOn=False 
lres@tiYAxisString = " "                      ;Total Carbon in all Terrestrial Carbon Pools
lres@gsnLeftString="(a) "
;lres@tiMainString = "Cropland expansion (CROP-CTL)"
;lres@trYMinF= min(Gpool2(0:2,:))-2
;lres@trYMaxF=max(Gpool)+2
plot(0)=gsn_csm_xy(wks,yrfrac4,Gpool(:,:),lres)
print(Gpool(:,159))
print(dim_avg_n_Wrap(Gpool(:,130:159), 1))
print(Gpool2(:,159))
print(dim_avg_n_Wrap(Gpool2(:,130:159), 1))
print(Gpool4(:,159))
print(dim_avg_n_Wrap(Gpool4(:,130:159), 1))
;lres@trYMinF= min(Gpool2(0:2,:))-2
;lres@trYMaxF=max(Gpool2(0:2,:))+2
;lres@xyLineOpacityF=0.45
lres@xyDashPattern=12
lres@gsnStringFont="times-roman"
plot2(0)=gsn_csm_xy(wks,yrfrac4,Gpool2(:,:),lres)
lres@xyLineOpacityF=1.
lres@xyDashPattern=2
lres@gsnLeftString="      Cropland expansion (CROP-CTL)"
plot3(0)=gsn_csm_xy(wks,yrfrac4,Gpool4(:,:),lres)
lres@gsnStringFont="times-Bold"
overlay(plot3(0), plot2(0))
overlay(plot(0), plot3(0))
;..........afforestation
lres@xyDashPattern=0
lres@tmYLTickStartF = -30 ;-15
lres@tmYLTickSpacingF = 15
lres@tmYLTickEndF = 45 ;30
lres@gsnLeftString="(b) "
;lres@tiMainString = "Afforestation (FRST-CTL)"
;lres@tiYAxisString = "Global nonlocal carbon change  ~C~                               (GtC)"                      ;Total Carbon in all Terrestrial Carbon Pools
lres@tiYAxisString = "Global nonlocal carbon change (GtC)"
;lres@trYMinF= min(Gflux2)-2
;lres@trYMaxF= max(Gflux)+2
plot(1)=gsn_csm_xy(wks,yrfrac4,Gflux(:,:),lres)
print(Gflux(:,159))
print(dim_avg_n_Wrap(Gflux(:,130:159), 1))
print(Gflux2(:,159))
print(dim_avg_n_Wrap(Gflux2(:,130:159), 1))
print(Gflux4(:,159))
print(dim_avg_n_Wrap(Gflux4(:,130:159), 1))
;lres@trYMinF= min(Gflux2)-2
;lres@trYMaxF=max(Gflux2)+2
;lres@xyLineOpacityF=0.45
lres@xyDashPattern=12
plot2(1)=gsn_csm_xy(wks,yrfrac4,Gflux2(:,:),lres)
lres@xyLineOpacityF=1.
lres@xyDashPattern=2
lres@gsnStringFont="times-roman"
lres@gsnLeftString="      Afforestation (FRST-CTL)"
plot3(1)=gsn_csm_xy(wks,yrfrac4,Gflux4(:,:),lres)
lres@gsnStringFont="times-Bold"

overlay(plot(1), plot2(1))
overlay(plot(1), plot3(1))

;..irrigation.................
lres@tmXBOn= True
lres@xyDashPattern=0
lres@gsnLeftString="(c) "
;lres@tiMainString = "Irrigation expansion (IRR-CROP)"
lres@tiYAxisString = " "                      ;Total Carbon in all Terrestrial Carbon Pools
lres@tiXAxisString = "Simulation year"
lres@tmXBOn        = True
lres@tmXBLabelsOn=True
;lres@trYMinF= min(Gpool3)-2
;lres@trYMaxF=max(Gpool32)+2
plot(2)=gsn_csm_xy(wks,yrfrac4,Gpool3(:,:),lres)
print(Gpool3(:,159))
print(dim_avg_n_Wrap(Gpool3(:,130:159), 1))
print(Gpool32(:,159))
print(dim_avg_n_Wrap(Gpool32(:,130:159), 1))
print(Gpool43(:,159))
print(dim_avg_n_Wrap(Gpool43(:,130:159), 1))
;lres@trYMinF= min(Gpool32)-2
;lres@trYMaxF=max(Gpool32)+2
;lres@xyLineOpacityF=0.45
lres@xyDashPattern=12
plot2(2)=gsn_csm_xy(wks,yrfrac4,Gpool32(:,:),lres)
lres@xyLineOpacityF=1
lres@xyDashPattern=2
lres@gsnStringFont=""
lres@gsnStringFont="times-roman"
lres@gsnLeftString="      Irrigation expansion (IRR-CROP)"
plot3(2)=gsn_csm_xy(wks,yrfrac4,Gpool43(:,:),lres)
overlay(plot(2), plot2(2))
overlay(plot(2), plot3(2))
lres@gsnStringFont="times-Bold"

lgres                    = True
lgres@lgLineThicknessF   = lres@xyLineThicknessF  ; legend line thickness = contour line thickness
;lgres@lgItemType         = "Lines"        ; show lines only (default)
lgres@lgLabelFontHeightF = .07           ; set the legend label font thickness
lgres@lgLineLabelFont="times-roman"

lgres@vpWidthF           = 0.1           ; width of legend (NDC)
lgres@vpHeightF          = 0.05            ; height of legend (NDC)
lgres@lgPerimOn      = False       ; draw the box perimeter in orange
lgres@lgMonoDashIndex =True
amres = True 
amres@amParallelPosF   = 0.6   
 ;Add text string.
;*****************************************************   
; Manually create and attach legend
;***************************************************** 
  res_text                    = True                  ; text mods desired
  res_text@txFontHeightF      = 0.018                 ; change text size
  res_text@txJust             = "CenterLeft"          ; text justification
  res_text@txFont="times-roman"

  res_lines                   = True                  ; polyline mods desired
  ;res_lines@gsLineDashPattern = 0.                  
  res_lines@gsLineThicknessF  = 2.3                    ; line thicker
  res_lines@gsLineColor       = "grey69"                 ; line color
  ;add box 1
  x_b = (/4,4,81,81,4/)
  y_b = (/-24,-10,-10,-24,-24/)
  dum1b = gsn_add_polyline(wks,plot(2),x_b,y_b,res_lines)  

  x_a = (/87,87,156,156,87/)
  y_a = (/-24,-10,-10,-24,-24/)
  dum1a = gsn_add_polyline(wks,plot(2),x_a,y_a,res_lines)  
  
 ;add box 2
  ;xx = (/87,87/)
  ;yy = (/-22,-13/)
  ;dum1b = gsn_add_polyline(wks,plot(2),xx,yy,res_lines)  
  ;  xx = (/87,158/)
  ;yy = (/-13,-13/)
  ;dum1b = gsn_add_polyline(wks,plot(2),xx,yy,res_lines)  
  ;  xx = (/158,158/)
  ;yy = (/-22,-13/)
  ;dum1b = gsn_add_polyline(wks,plot(2),xx,yy,res_lines)  
  ;  xx = (/87,158/)
  ;yy = (/-22,-22/)
  ;dum1b = gsn_add_polyline(wks,plot(2),xx,yy,res_lines) 
res_lines@gsLineThicknessF  = lres@xyLineThicknessF                    ; line thicker
res_lines@gsLineColor       = "palegreen4"
  ;xx = (/4,15/)
  xx = (/5,16/)
  yy = (/-14,-14/)
  dum1 = gsn_add_polyline(wks,plot(2),xx,yy,res_lines)              ; add polyline
  dum2 = gsn_add_text(wks,plot(2),"Vegetation C",18,-14,res_text); add text
  
  res_lines@gsLineColor       = "sandybrown"                 ; line color
  ;xx = (/4,15/)
  xx = (/5,16/)
  yy = (/-20,-20/)
  dum3 = gsn_add_polyline(wks,plot(2),xx,yy,res_lines)              ; add polyline
  dum4 = gsn_add_text(wks,plot(2),"Soil C",18,-20,res_text); add text

  res_lines@gsLineColor       = "royalblue"                 ; line color
  ;xx = (/46,57/)
  xx = (/45,56/)
  yy = (/-14,-14/)
  dum5 = gsn_add_polyline(wks,plot(2),xx,yy,res_lines)              ; add polyline
  dum6 = gsn_add_text(wks,plot(2),"Litter C",58,-14,res_text); add text

  res_lines@gsLineColor       = "black"                 ; line color
  ;xx = (/46,57/)
  xx = (/45,56/)
  yy = (/-20,-20/)
  dum7 = gsn_add_polyline(wks,plot(2),xx,yy,res_lines)              ; add polyline
  dum8 = gsn_add_text(wks,plot(2),"Total land C",58,-20,res_text); add text
 ; "Litter C","Total land C"/),(/15,15,54.4,54.4/),(/-11.2,-22.3,-11.2,-22.3/),txres)

  res_lines@gsLineColor       = "grey55"                 ; line color
  ;xx = (/86,97/)
  xx = (/89,100/)
  yy = (/-14,-14/)
  dum9 = gsn_add_polyline(wks,plot(2),xx,yy,res_lines)              ; add polyline
  dum10 = gsn_add_text(wks,plot(2),"MPI-ESM",102,-14,res_text); add text

  res_lines@gsLineColor       = "grey55"                 ; line color
  res_lines@gsLineDashPattern = 12
  ;xx = (/86,97/)
  xx = (/89,100/)
  yy = (/-20,-20/)
  dum11 = gsn_add_polyline(wks,plot(2),xx,yy,res_lines)              ; add polyline
  dum12 = gsn_add_text(wks,plot(2),"CESM",102,-20,res_text); add text

  res_lines@gsLineColor       = "grey55"                 ; line color
  res_lines@gsLineDashPattern = 2
  res_lines@gsLineOpacityF= 1.
  xx = (/126,137/)
  yy = (/-14,-14/)
  dum13 = gsn_add_polyline(wks,plot(2),xx,yy,res_lines)              ; add polyline
  dum14 = gsn_add_text(wks,plot(2),"EC-Earth",139,-14,res_text); add text

lnres=True
lnres@gsnDraw         = False                                                          ;         do not draw the plot                
lnres@gsnFrame        = False     
lnres@vpHeightF= 0.4                    ; change aspect ratio of plot
lnres@vpWidthF = 0.8                  
lnres@tiMainString = ""
lnres@xyLineThicknessF = 2.8
lnres@xyDashPatterns=0
lnres@tiYAxisString = " "
lnres@gsnLeftStringFontHeightF = 0.03
lnres@gsnPanelYWhiteSpacePercent = 0.7
;lnres@gsnPanelYWhiteSpacePercent = 0.9
lnres@gsnPanelBottom=0.5
midi=(/0,0/)
 lnres@xyLineColors = (/"gray60","gray60"/);
 do i=0,2
  ;fill_xy_plot_area(wks,plot(i))
 ;plotn(i)=gsn_add_polyline(wks, plot(i),(/0,159/), midi, lnres)
end do

 pres                  = True         ; Set panel resources.
 ;pres@gsnMaximize      = True         ; Maximize plots in panel.
 pres@gsnPaperMargin = 0
  pres@gsnBoxMargin  = 0

pres@gsnPanelYWhiteSpacePercent = 0.01
pres@gsnPanelXWhiteSpacePercent = 0
pres@gsnPanelBottom=0.03
pres@gsnPanelRight = 1
pres@gsnPanelLeft = 0.03
;drawNDCGrid(wks)
;psres = True                                                               
 ;maximize_output(wks,psres)  ; calls draw and frame for you 

;
; If you try to use gsn_panel to panel different sized plots, you
; will get unexpected results. You can turn on gsnPanelDebug to
; see what viewport values are used, and then use these values
; to set the viewport resources for each individual plot.
 ;pres               = True
 ;pres@gsnPanelDebug = True    ; Use to print values for position, width,

                                ; and height.
gsn_panel(wks,plot(:),(/3,1/),pres)

end