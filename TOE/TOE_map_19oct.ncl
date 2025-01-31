load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin
;......................................................................
myfileCT = addfile("/work/bm1147/b380949/web-monitoring/CTL_MPI/cVeg_ctl.nc", "r")
myfileCT1 = addfile("/work/bm1147/b380949/web-monitoring/CTL_MPI/cSoil_ctl.nc", "r")

myfile1 = addfile("/home/b/b380949/TOE_1std.nc", "r")
myfile2 = addfile("/home/b/b380949/TOE_1std_CESM.nc", "r")
myfile3 = addfile("/home/b/b380949/TOE_1std_ECEARTH.nc", "r")


TOE_cveg_crop=myfile1->TOE_cveg_crop(15,:,:)
TOE_csoil_crop=myfile1->TOE_csoil_crop(15,:,:)
TOE_cveg_frst=myfile1->TOE_cveg_frst(15,:,:)
TOE_csoil_frst=myfile1->TOE_csoil_frst(15,:,:) 
TOE_cveg_irri=myfile1->TOE_cveg_irri(15,:,:)
TOE_csoil_irri=myfile1->TOE_csoil_irri(15,:,:)

TOE_cveg_crop_cesm=myfile2->TOE_cveg_crop(15,:,:)
TOE_csoil_crop_cesm=myfile2->TOE_csoil_crop(15,:,:)
TOE_cveg_frst_cesm=myfile2->TOE_cveg_frst(15,:,:)
TOE_csoil_frst_cesm=myfile2->TOE_csoil_frst(15,:,:) 
TOE_cveg_irri_cesm=myfile2->TOE_cveg_irri(15,:,:)
TOE_csoil_irri_cesm=myfile2->TOE_csoil_irri(15,:,:)

TOE_cveg_crop_ecearth=myfile3->TOE_cveg_crop(15,:,:)
TOE_csoil_crop_ecearth=myfile3->TOE_csoil_crop(15,:,:)
TOE_cveg_frst_ecearth=myfile3->TOE_cveg_frst(15,:,:)
TOE_csoil_frst_ecearth=myfile3->TOE_csoil_frst(15,:,:) 
TOE_cveg_irri_ecearth=myfile3->TOE_cveg_irri(15,:,:)
TOE_csoil_irri_ecearth=myfile3->TOE_csoil_irri(15,:,:)

wks_type = "pdf"
wks_type@wkPaperWidthF =  9.0
wks_type@wkPaperHeightF = 9.0
wks = gsn_open_wks(wks_type,"TOE_map_19oct")
gsn_define_colormap(wks,"cmocean_turbid")
res = True
res@cnFillOn  = True
res@cnLinesOn = False
res@gsnPaperMargin = 0
res@wkBackgroundColor="white"
res@cnFillMode = "RasterFill"
res@gsnAddCyclic =True
res@vpWidthF   = 0.5
res@vpHeightF  =0.2
res@tmXTOn                = False 
res@tmXTLabelsOn          = False 
res@tmYROn                = False
res@tmYRLabelsOn          = False
res@cnLineThicknessF = 1.0
res@cnLineLabelsOn = False
res@tmYLMinorOn           = True
res@lbLabelBarOn         = False ; Turn off individual labelbar;
;res@lbBoxEndCapStyle    = "RectangleEnds";"TriangleBothEnds"
res@lbTitleOffsetF      =-0.03
res@lbTitleJust         ="CenterCenter"
res@lbTitlePosition     = "Bottom"
res@mpMinLatF            = -60
res@gsnDraw = False
res@gsnFrame = False
res@mpGeophysicalLineColor="gray70"
res@mpGeophysicalLineThicknessF=0.6
res@tmXBBorderOn=False
res@tmXTBorderOn=False
res@tmYLBorderOn=False
res@tmYRBorderOn=False
res@tmXBOn=False
res@tmXTOn=False
res@tmYLOn=False
res@tmYROn=False
res@cnInfoLabelOn = False
res@cnLevelSelectionMode = "ExplicitLevels"
res@cnLabelBarEndStyle= "ExcludeOuterBoxes"
res@cnLevels  =  (/0,10,20,30,40,70,100,130,160/);17
res@cnFillColors =   (/0,2,41,65,97,145,178,225,255/)

lres=True
lres@gsnDraw         = False                  ;         do not draw the plot                
lres@gsnFrame        = False     
lres@gsnPaperMargin = 0
lres@vpHeightF= res@vpHeightF                    ; change aspect ratio of plot
lres@vpWidthF = 0.42*res@vpWidthF                  
lres@xyLineThicknessF = 1.5
lres@tmXTOn        = False ; turn off top tickmarks
lres@tmYLOn        = False ; turn off right tickmarks
lres@tmBorderThicknessF= 0.8
lres@tmXBLabelFontHeightF = 0.014
;lres@tmXBLabelAngleF = 15
lres@tmYRLabelFontHeightF = 0.014
lres@tmXBLabelFont="times-roman"
lres@tmXBLabelDeltaF = -0.7
lres@gsnLeftStringParallelPosF= 0.02
lres@gsnLeftStringOrthogonalPosF=-0.14
lres@gsnLeftStringFontHeightF=0.018
lres@tmXBMajorThicknessF = 0.5
lres@tmYRMajorThicknessF = 0.5

lres@tmXBMajorOutwardLengthF = 0.0
lres@tmYRMajorOutwardLengthF=0.0
lres@tmXBMajorLengthF=0.005
lres@tmYRMajorLengthF=0.005

lres@tmYRLabelsOn =True
lres@tmYUseLeft =False
lres@tmYRLabelFont="times-roman"
lres@tmYRMode              = "Explicit"
lres@tmYRValues            = ispan(-60,80,20)
lres@tmYRLabels  = (/"60~S~o~N~S","40~S~o~N~S","20~S~o~N~S","0~S~o~N","20~S~o~N~N",\  
                     "40~S~o~N~N","60~S~o~N~N","80~S~o~N~N"/) 
printVarSummary(TOE_cveg_irri({-60:90},:))
printVarSummary(TOE_cveg_crop_cesm({-60:90},:))
  y = TOE_cveg_crop&lat   ;y(16::)  -60to90 80points
y_cesm = TOE_cveg_crop_cesm&lat
y_ecearth = TOE_cveg_crop_ecearth&lat

crop_zonal=dim_avg_n_Wrap(TOE_cveg_crop({-60:90},:),1)
crop_zonal_cesm=dim_avg_n_Wrap(TOE_cveg_crop_cesm({-60:90},:),1)
crop_zonal_ecearth=dim_avg_n_Wrap(TOE_cveg_crop_ecearth({-60:90},:),1)

frst_zonal=dim_avg_n_Wrap(TOE_cveg_frst({-60:90},:),1)
frst_zonal_cesm=dim_avg_n_Wrap(TOE_cveg_frst_cesm({-60:90},:),1)
frst_zonal_ecearth=dim_avg_n_Wrap(TOE_cveg_frst_ecearth({-60:90},:),1)

irri_zonal=dim_avg_n_Wrap(TOE_cveg_irri({-60:90},:),1)
irri_zonal_cesm=dim_avg_n_Wrap(TOE_cveg_irri_cesm({-60:90},:),1)
irri_zonal_ecearth=dim_avg_n_Wrap(TOE_cveg_irri_ecearth({-60:90},:),1)

;...plot.............................
plot1=new(12,graphic)
plot2=new(12,graphic)
plot3=new(12,graphic)

res@tiXAxisFont="times-roman"
res@tiXAxisFontHeightF=0.01
res@tiYAxisFont="times-roman"
res@tiMainFont="times-roman"
res@tiMainFontHeightF=0.019  ; 0.018
res@tiMainOffsetYF = -0.008
res@tiYAxisFontHeightF=0.019
res@gsnLeftStringParallelPosF= 0.04
res@gsnLeftStringOrthogonalPosF=-0.11
 res@gsnRightString    = " "
res@tmXBLabelFontHeightF =0.014
res@tmYLLabelFontHeightF =0.014
res@gsnStringFont ="times-Bold"
res@gsnLeftStringFontHeightF=0.019  ;0.018
res@gsnRightStringFontHeightF=0.011
res@gsnCenterStringFontHeightF=0.018
res@tmXBLabelFont         = "times-roman"
res@tmYLLabelFont         = "times-roman" 

    lres@xyLineColor= (/0.84,0.52,0.22,1./); "darkgoldenrod3"
   res@cnMonoFillColor=False
     res@tiYAxisString  = "CROP-CTL"
     res@tiMainString ="MPI-ESM"
    res@gsnLeftString    = "(a)"
    plot3(0)=gsn_csm_contour_map(wks, TOE_cveg_crop(:,:), res)
    res@tiYAxisString  = " "
    res@tiMainString ="CESM"
    res@gsnLeftString    = "(b)"
    plot3(1)=gsn_csm_contour_map(wks,TOE_cveg_crop_cesm, res)

    res@tiMainString ="EC-Earth"
    res@gsnLeftString    = "(c)"
    plot3(2)=gsn_csm_contour_map(wks,TOE_cveg_crop_ecearth, res)
        res@gsnCenterString = ""
    res@tiMainString =""

    lres@gsnStringFont ="times-Bold"
    lres@gsnLeftStringFontHeightF=0.019
    lres@gsnLeftString    = "(d)"
    lres@trXMinF= -20
    lres@trXMaxF=160
    lres@tmXBMode              = "Explicit"
    lres@tmXBValues            =(/0,30,60,90,120,150/)   
    lres@tmXBLabels  = (/"0","30","60","90","120","150"/)  

    plot3(3)=gsn_csm_xy(wks,crop_zonal_ecearth(:),y_ecearth(32::),lres)

res@tiYAxisString  = "FRST-CTL"
    res@gsnLeftString    = "(d)"
    plot3(4)=gsn_csm_contour_map(wks, TOE_cveg_frst(:,:), res)
        res@tiYAxisString  = " "
    res@gsnLeftString    = "(e)"
    plot3(5)=gsn_csm_contour_map(wks, TOE_cveg_frst_cesm(:,:), res)

          res@gsnLeftString    = "(f)"

    plot3(6)=gsn_csm_contour_map(wks, TOE_cveg_frst_ecearth(:,:), res)

lres@gsnLeftString    = "(h)"
        lres@tiXAxisString=" "
    plot3(7)=gsn_csm_xy(wks,frst_zonal_ecearth(:),y_ecearth(32::),lres)
   res@tiYAxisString  = "IRR-CROP"
    res@gsnLeftString    = "(g)"
    plot3(8)=gsn_csm_contour_map(wks, TOE_cveg_irri(:,:), res)
    res@tiYAxisString  = " "
    res@gsnLeftString    = "(h)"
    plot3(9)=gsn_csm_contour_map(wks, TOE_cveg_irri_cesm(:,:), res)

    res@gsnLeftString    = "(i)"
    plot3(10)=gsn_csm_contour_map(wks, TOE_cveg_irri_ecearth(:,:), res)

    lres@gsnLeftString    = "(l)"
    lres@tiXAxisString= "";"kg/m~S~2~N~"
    lres@tiXAxisFont="times-roman"
    lres@tiXAxisFontHeightF=0.018

    plot3(11)=gsn_csm_xy(wks,irri_zonal_ecearth(:),y_ecearth(32::),lres)
    lres@tiXAxisString=" "
;..............................................
    lres@xyLineColor= (/0.0,0.0,0.0,1.0/)  ; "black"
    lres@gsnLeftString    = " "
    lres@tiXAxisString=" "
    plot1(3)=gsn_csm_xy(wks,crop_zonal_cesm(:),y_cesm(32::),lres)
    plot1(7)=gsn_csm_xy(wks,frst_zonal_cesm(:),y_cesm(32::),lres)
    plot1(11)=gsn_csm_xy(wks,irri_zonal_cesm(:),y_cesm(32::),lres)
       overlay(plot3(3),plot1(3))
       overlay(plot3(7),plot1(7))
       overlay(plot3(11),plot1(11))
        lres@xyLineColor= (/0.39,0.66,0.26,1./) ;"green3"
    plot2(3)=gsn_csm_xy(wks,crop_zonal(:),y(16::),lres)
    plot2(7)=gsn_csm_xy(wks,frst_zonal(:),y(16::),lres)
    plot2(11)=gsn_csm_xy(wks,irri_zonal(:),y(16::),lres)
       overlay(plot3(3),plot2(3))
       overlay(plot3(7),plot2(7))
       overlay(plot3(11),plot2(11))
lgres                    = True
lgres@lgLineThicknessF   = lres@xyLineThicknessF  ; legend line thickness = contour line thickness
;lgres@lgItemType         = "Lines"        ; show lines only (default)
lgres@lgLabelFontHeightF = 0.07          ; set the legend label font thickness
lgres@lgLineLabelFont="times-roman"
lgres@lgLabelJust = "CenterLeft"
lgres@lgLabelFont = "times-roman"
lgres@vpWidthF           = 0.1           ; width of legend (NDC)
lgres@vpHeightF          = 0.05            ; height of legend (NDC)
lgres@lgPerimOn      = False       ; draw the box perimeter in orange
lgres@lgMonoDashIndex =True
amres = True 
amres@amParallelPosF   = 0.6    
strk=(/" "," "," "/)
strk1=(/"  MPI-ESM","  CESM"," EC-Earth"/)

colors=(/"royalblue","orangered4","gold3"/)
lgres@lgLineColors       = (/(/0.39,0.66,0.26,1./),(/0.0,0.0,0.0,1./),(/0.84,0.52,0.22,1./)/)
gsn_legend_ndc(wks,3,(/"  MPI-ESM","  CESM","  EC-Earth"/),0.73,0.33,lgres)


  plres = True
 ;plres@gsFillOpacityF         =  0.5
plres@gsLineThicknessF = 0.7
 plres@gsLineOpacityF         =  3

plot7=new(7,graphic)
plx=new((/7,5/),"integer")
ply=new((/7,5/),"integer")
;North_America{25:49},{255:300} ;Amazon{-18:2},{283:307} ;Congo{-9:9},{13:29} 
;North Eurasia{53:68},{30:90};East Asia{40:57},{107:139};Southeast Asia {15:30},{105:117};Australia {-23:-12},{115:150}
plx = (/(/260,260,300,300,260/),(/283,283,308,308,283/),(/13,13,29,29,13/),(/32,32,92,92,32/),(/106,106,140,140,106/),(/101,101,119,119,101/),(/115,115,150,150,115/)/)
ply = (/(/28,54,54,28,28/),(/-16,3,3,-16,-16/),(/-9,9,9,-9,-9/),(/56,67,67,56,56/),(/46,62,62,46,46/),(/22,36,36,22,22/),(/-23,-12,-12,-23,-23/)/)

do i= 0,6
    plres@gsLineColor            = "red"
     ;   plot7(i)=gsn_add_polyline(wks, plot3(0,0), plx(i,:), ply(i,:), plres)
end do

pres                  = True         ; Set panel resources.
pres@gsnPanelBottom=0. ;0.03
pres@gsnPanelRight = 1.; 1
pres@gsnPanelLeft = 0. ;0.03
pres@gsnPaperMargin = 0.0
pres@gsnBoxMargin  = 0.0
pres@gsnPanelLabelBar = True         ; Turn on panel labelbar.
pres@pmLabelBarWidthF    = 0.6      ; default is shorter
pres@pmLabelBarHeightF   = 0.04       ; default is taller
pres@lbLabelFontHeightF  = 0.01       ; make labels smaller ( default=0.02 )
pres@lbPerimOn             = False
pres@pmLabelBarOrthogonalPosF = 0.00 ;0.02
pres@pmLabelBarParallelPosF   =  -0.1  ;0.02
pres@lbBoxLineColor        = "white"
pres@lbBoxEndCapStyle    = "RectangleEnds" ;"TriangleBothEnds"

pres@lbTitleFontHeightF  =0.01
pres@lbTitleFont         ="times-roman"
pres@lbTitleString       ="Time of emergence (years)  "
pres@lbTitleOffsetF      = 0.2 ;-0.2
pres@lbTitleJust         ="CenterCenter"
pres@lbTitlePosition     = "Bottom"
pres@gsnPanelYWhiteSpacePercent = 0.005
pres@gsnPanelXWhiteSpacePercent =0.0
pres@tiYAxisString=""
pres@txFontHeightF  =0.018
pres@txFont ="times-roman"
pres@gsnPanelXF = (/0.0185,0.255,0.492,0.729,0.0185,0.255,0.492,0.729,0.0185,0.255,0.492,0.729/) ;0.4
pres@gsnPanelYF = (/0.664,0.664,0.664,0.664,0.554,0.554,0.554,0.554,0.444,0.444,0.444,0.444/) ;0.4
pres@gsnPanelDebug = True
gsn_panel(wks,plot3(:),(/3,4/),pres)
end
