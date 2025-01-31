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
wks_type@wkPaperWidthF =  8.0
wks_type@wkPaperHeightF = 8.0
wks = gsn_open_wks(wks_type,"TOE_latitudinal_cSoil")
gsn_define_colormap(wks,"cmocean_turbid")
lres=True
lres@gsnDraw         = False                  ;         do not draw the plot                
lres@gsnFrame        = False     
lres@gsnPaperMargin = 0
lres@vpHeightF= 0.6                  ; change aspect ratio of plot
lres@vpWidthF = 0.6                 
lres@xyLineThicknessF = 2;1.5
lres@tmXTOn        = False ; turn off top tickmarks
lres@tmYLOn        = True ; turn off Left tickmarks
lres@tmYROn        = False ; turn off right tickmarks
lres@tmYLMinorOn =False
lres@tmBorderThicknessF=1.0; 0.8
lres@tmXBLabelFontHeightF = 0.024; 0.014
lres@tmYLLabelFontHeightF =0.024; 0.014
lres@tmXBLabelFont="times-roman"
lres@tmXBLabelDeltaF = -0.7
lres@gsnLeftStringParallelPosF= 0.02
lres@gsnLeftStringOrthogonalPosF=-0.14
lres@gsnLeftStringFontHeightF=0.018
lres@tmXBMajorThicknessF = 1.
lres@tmYLMajorThicknessF = 1.

lres@tmXBMajorOutwardLengthF = 0.0
lres@tmYLMajorOutwardLengthF= 0.0
lres@tmXBMajorLengthF=0.01
lres@tmYLMajorLengthF=0.01

printVarSummary(TOE_csoil_irri({-60:90},:))
printVarSummary(TOE_csoil_crop_cesm({-60:90},:))
  y = TOE_csoil_crop&lat   ;y(16::)  -60to90 80points
y_cesm = TOE_csoil_crop_cesm&lat
y_ecearth = TOE_csoil_crop_ecearth&lat

crop_zonal=dim_avg_n_Wrap(TOE_csoil_crop({-60:90},:),1)
crop_zonal_cesm=dim_avg_n_Wrap(TOE_csoil_crop_cesm({-60:90},:),1)
crop_zonal_ecearth=dim_avg_n_Wrap(TOE_csoil_crop_ecearth({-60:90},:),1)

frst_zonal=dim_avg_n_Wrap(TOE_csoil_frst({-60:90},:),1)
frst_zonal_cesm=dim_avg_n_Wrap(TOE_csoil_frst_cesm({-60:90},:),1)
frst_zonal_ecearth=dim_avg_n_Wrap(TOE_csoil_frst_ecearth({-60:90},:),1)

irri_zonal=dim_avg_n_Wrap(TOE_csoil_irri({-60:90},:),1)
irri_zonal_cesm=dim_avg_n_Wrap(TOE_csoil_irri_cesm({-60:90},:),1)
irri_zonal_ecearth=dim_avg_n_Wrap(TOE_csoil_irri_ecearth({-60:90},:),1)
;...plot.............................
plot1=new(12,graphic)
plot2=new(12,graphic)
plot3=new(12,graphic)
     lres@tiMainFont = "times-roman"
lres@gsnLeftStringFontHeightF=0.026
lres@gsnLeftStringOrthogonalPosF=0.009
lres@gsnLeftStringParallelPosF= 0.0
lres@gsnStringFont="times-Bold"

      lres@xyLineColor= (/0.75,0.20,0.1,1./) ;red 
lres@tmYLLabelDeltaF = 0
lres@tmYLLabelFont="times-roman"
lres@tmYLMode              = "Explicit"
lres@tmYLValues            = ispan(-60,80,20)
lres@tmYLLabels  = (/"60~S~o~N~S","40~S~o~N~S","20~S~o~N~S","0~S~o~N","20~S~o~N~N",\  
                     "40~S~o~N~N","60~S~o~N~N","80~S~o~N~N"/) 
    lres@tiMainFontHeightF = 0.026
     lres@tiMainOffsetYF = -0.024
     lres@tiMainOffsetXF = 0.01
     lres@tiMainString ="Cropland expansion (CROP-CTL)"
    lres@gsnStringFont ="times-Bold"
    lres@gsnLeftString    = "(a)"
    lres@trXMinF= 0
    lres@trXMaxF=160
    lres@tmXBMode              = "Explicit"
    lres@tmXBValues            =(/0,30,60,90,120,150/)   
    lres@tmXBLabels  = (/"0","30","60","90","120","150"/)  
  
    plot3(0)=gsn_csm_xy(wks,crop_zonal_ecearth(:),y_ecearth(32::),lres)

lres@tmYLLabels  = (/"","","","","",\  
                     "","",""/)
lres@tiMainString  = "Afforestation (FRST-CTL)"
lres@tiMainOffsetXF = -0.015
lres@gsnLeftString    = "(b)"
lres@tiXAxisFontHeightF=0.027
lres@tiXAxisFont ="times-roman"
lres@tiXAxisString="Latitudinal mean time of emergence (years)"
lres@tiXAxisOffsetYF= -0.01
    plot3(1)=gsn_csm_xy(wks,frst_zonal_ecearth(:),y_ecearth(32::),lres)

       lres@tiMainOffsetXF = 0.01
   lres@tiMainString  = "Irrigation expansion (IRR-CROP)"
    lres@gsnLeftString    = "(c)"
    lres@tiXAxisString= "";"kg/m~S~2~N~"
    lres@tiXAxisFont="times-roman"
    lres@tiXAxisFontHeightF=0.018
    plot3(2)=gsn_csm_xy(wks,irri_zonal_ecearth(:),y_ecearth(32::),lres)
    lres@tiXAxisString=" "
;..............................................
    lres@xyLineColor= (/0.0,0.0,0.0,1.0/)  ; "black"
    lres@gsnLeftString    = " "
    lres@tiXAxisString=" "
    plot1(0)=gsn_csm_xy(wks,crop_zonal_cesm(:),y_cesm(32::),lres)
    plot1(1)=gsn_csm_xy(wks,frst_zonal_cesm(:),y_cesm(32::),lres)
    plot1(2)=gsn_csm_xy(wks,irri_zonal_cesm(:),y_cesm(32::),lres)
     do i=0,2
       overlay(plot3(i),plot1(i))
     end do
        lres@xyLineColor= (/0.11,0.44,0.93,1./) ;blue 

    plot2(0)=gsn_csm_xy(wks,crop_zonal(:),y(16::),lres)
    plot2(1)=gsn_csm_xy(wks,frst_zonal(:),y(16::),lres)
    plot2(2)=gsn_csm_xy(wks,irri_zonal(:),y(16::),lres)
        do i=0,2
       overlay(plot3(i),plot2(i))
     end do

lgres                    = True
lgres@lgLineThicknessF   = lres@xyLineThicknessF  ; legend line thickness = contour line thickness
lgres@lgLabelFontHeightF = 0.08          ; set the legend label font thickness
lgres@lgLineLabelFont="times-roman"
lgres@lgLabelJust = "CenterLeft"
lgres@lgLabelFont = "times-roman"
lgres@vpWidthF           = 0.13           ; width of legend (NDC)
lgres@vpHeightF          = 0.07            ; height of legend (NDC)
lgres@lgPerimOn      = False       ; draw the box perimeter in orange
lgres@lgMonoDashIndex =True
amres = True 
amres@amParallelPosF   = 0.6    
strk=(/" "," "," "/)
strk1=(/"  MPI-ESM","  CESM"," EC-Earth"/)

colors=(/"royalblue","orangered4","gold3"/)
lgres@lgLineColors       = (/(/0.75,0.20,0.1,1./),(/0.0,0.0,0.0,1./),(/0.11,0.44,0.93,1./)/)
gsn_legend_ndc(wks,3,(/"  EC-Earth","  CESM","  MPI-ESM"/),0.825,0.45,lgres)

     pnlres                = True
     pnlres@gsnMaximize    = True   ; Maximize in frame
     pnlres@gsnPaperMargin =0 
  pnlres@gsnBoxMargin  = 0
pnlres@gsnPanelRight = 1.0
pnlres@gsnPanelLeft =0.01
pnlres@gsnPanelXF = (/0.07,0.375,0.68/)
     pnlres@gsnPanelYWhiteSpacePercent = 0.00
     pnlres@gsnPanelXWhiteSpacePercent = 0.00
    pnlres@gsnPanelDebug = True

     pnlres@gsnPanelFigureStringsFontHeightF=0.02
     pnlres@amJust="TopLeft"
     pnlres@gsnPanelFigureStringsPerimOn      = False  

gsn_panel(wks,plot3(0:2),(/1,3/),pnlres)
end
