load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin
;..When plotting different scenarios, change myfile1,2 to nyfile5,6 or myfile13,14.........................
;...When plotting different scenarios, change myfile12,22 to myfile52,62 or myfile132,142................
landsea_ecearth=addfile("/work/bm1147/b380949/EC-earth/dataset/anlysisneed/landmask_ecearth_regridded.nc","r")
landmask_ecearth=landsea_ecearth->landmask(0,{-60:80},:)

myfile1 = addfile("/home/b/b380949/mulreg/attri_crop_cVeg.nc", "r")
myfile2 = addfile("/home/b/b380949/mulreg/attri_crop_cSoil.nc", "r")

myfile5 = addfile("/home/b/b380949/mulreg/attri_frst_cveg.nc", "r")
myfile6 = addfile("/home/b/b380949/mulreg/attri_frst_csoil.nc", "r")
myfile13 = addfile("/home/b/b380949/mulreg/attri_irri_cveg.nc", "r")
myfile14 = addfile("/home/b/b380949/mulreg/attri_irri_csoil.nc", "r")

myfile12 = addfile("/home/b/b380949/mulreg/attri_crop_cVeg_CESM.nc", "r")
myfile22 = addfile("/home/b/b380949/mulreg/attri_crop_cSoil_CESM.nc", "r")

myfile52 = addfile("/home/b/b380949/mulreg/attri_frst_cVeg_CESM.nc", "r")
myfile62 = addfile("/home/b/b380949/mulreg/attri_frst_cSoil_CESM.nc", "r")
myfile132 = addfile("/home/b/b380949/mulreg/attri_irri_cVeg_CESM.nc", "r")
myfile142 = addfile("/home/b/b380949/mulreg/attri_irri_cSoil_CESM.nc", "r")

myfile123= addfile("/home/b/b380949/mulreg/attri_crop_cVeg_ecearth.nc", "r")
myfile223= addfile("/home/b/b380949/mulreg/attri_crop_cSoil_ecearth.nc", "r")

myfile523= addfile("/home/b/b380949/mulreg/attri_frst_cVeg_ecearth.nc", "r")
myfile623= addfile("/home/b/b380949/mulreg/attri_frst_cSoil_ecearth.nc", "r")
myfile1323= addfile("/home/b/b380949/mulreg/attri_irri_cVeg_ecearth.nc", "r")
myfile1423= addfile("/home/b/b380949/mulreg/attri_irri_cSoil_ecearth.nc", "r")

;.........deforestation..... ..................
Gpool=new((/2,75,192/),double) 
Gpool(0,:,:) = myfile5->b(1,:,:) ;cVeg_mpi_t2m
Gpool(1,:,:) =myfile6->b(1,:,:)   ;csoil_mpi_t2m
;..............afforestation...................
Gpool2=new((/2,75,192/),double) 
Gpool2(0,:,:) = myfile5->b(2,:,:)  ;cVeg_mpi_surface soil moisture
Gpool2(1,:,:) =myfile6->b(2,:,:)   ;csoil_mpi_soil moisture
;.....................irrigation......................................
Gpool4=new((/2,75,192/),double) 
Gpool4(0,:,:) = myfile5->correlation   ;correlation
Gpool4(1,:,:) =myfile6->correlation
;................white for diff=0........................
Gpool_revised=mask(Gpool, Gpool.eq.0, False)*100   ;crop(5,75,192)5:veg,soil,product,litter,land
copy_VarMeta(Gpool, Gpool_revised)
Gpool2_revised=mask(Gpool2, Gpool2.eq.0, False)*100  ;frst
copy_VarMeta(Gpool2, Gpool2_revised)
Gpool4_revised=mask(Gpool4, Gpool4.eq.0, False)  ;irri
copy_VarMeta(Gpool4, Gpool4_revised)
;...............CESM........t2m..
pool22=new((/2,149,288/),double) 
pool22(0,:,:) = myfile52->b(1,:,:)
pool22(1,:,:) =myfile62->b(1,:,:)
;.............mrsos..........
flux2=new((/2,149,288/),double)  
flux2(0,:,:) = myfile52->b(2,:,:)
flux2(1,:,:) =myfile62->b(2,:,:)
;................cor...................
pool23=new((/2,149,288/),double)  
pool23(0,:,:)= myfile52->correlation
pool23(1,:,:)=myfile62->correlation
;...................
pool2_revised=mask(pool22, pool22.eq.0, False)*100;crop;4,149,288;4:veg,soil,litter,land
copy_VarMeta(pool22, pool2_revised)

pool23_revised=mask(pool23, pool23.eq.0, False) ;irri
copy_VarMeta(pool23, pool23_revised)

flux2_revised=mask(flux2, flux2.eq.0, False)*100;frst
copy_VarMeta(flux2, flux2_revised)
;.........ecearth.....t2m..
bt2m_eceartho=new((/2,149,288/),double) 
bt2m_eceartho(0,:,:) = myfile523->b(1,:,:)
bt2m_eceartho(1,:,:) =myfile623->b(1,:,:)
;.............mrsos..........
bmrsos_eceartho=new((/2,149,288/),double)  
bmrsos_eceartho(0,:,:) = myfile523->b(2,:,:)
bmrsos_eceartho(1,:,:) =myfile623->b(2,:,:)
;................cor...................
cor_eceartho=new((/2,149,288/),double)  
cor_eceartho(0,:,:)= myfile523->correlation
cor_eceartho(1,:,:)=myfile623->correlation
;.......deal with missing value...........
bt2m_eceartho=mask(bt2m_eceartho,isnan_ieee(bt2m_eceartho),False)
bmrsos_eceartho=mask(bmrsos_eceartho,isnan_ieee(bmrsos_eceartho),False)
cor_eceartho=mask(cor_eceartho,isnan_ieee(cor_eceartho),False)
print(any(isnan_ieee(bt2m_eceartho)))
print(any(isnan_ieee(bmrsos_eceartho)))
print(any(isnan_ieee(cor_eceartho)))
bt2m_ecearth=mask(bt2m_eceartho, landmask_ecearth.eq.1, True)
copy_VarMeta(bt2m_eceartho, bt2m_ecearth)
bmrsos_ecearth=mask(bmrsos_eceartho, landmask_ecearth.eq.1, True)
copy_VarMeta(bmrsos_eceartho, bmrsos_ecearth)
cor_ecearth=mask(cor_eceartho, landmask_ecearth.eq.1, True)
copy_VarMeta(cor_eceartho, cor_ecearth)
;..................
bt2m_ecearth_revised=mask(bt2m_ecearth, bt2m_ecearth.eq.0, False)*100;crop;4,149,288;4:veg,soil,litter,land
copy_VarMeta(bt2m_ecearth, bt2m_ecearth_revised)

cor_ecearth_revised=mask(cor_ecearth, cor_ecearth.eq.0, False) ;irri
copy_VarMeta(cor_ecearth, cor_ecearth_revised)

bmrsos_ecearth_revised=mask(bmrsos_ecearth, bmrsos_ecearth.eq.0, False)*100;frst
copy_VarMeta(bmrsos_ecearth, bmrsos_ecearth_revised)
;......plot.........................................
wks_type = "pdf"
wks_type@wkPaperWidthF =  9.0
wks_type@wkPaperHeightF = 9.0
wks = gsn_open_wks(wks_type,"attri_cVeg_frst_21oct")
gsn_define_colormap(wks,"BrYlGr")
res = True
res@cnFillOn  = True
res@cnLinesOn = False
res@wkBackgroundColor="white"
res@cnFillMode = "RasterFill"
res@gsnAddCyclic =True
res@vpWidthF   = 0.5
res@vpHeightF  =0.2
res@tmXTOn                = False 
res@tmXTLabelsOn          = False 
res@tmYROn                = False
res@tmYRLabelsOn          = False
res@cnLevelSelectionMode = "ExplicitLevels"
res@cnLineThicknessF = 0
res@cnLineLabelsOn = False
res@tmYLMinorOn           = True
res@lbLabelBarOn         = True ; Turn off individual labelbar;
res@pmLabelBarHeightF   = 0.04 
res@lbLabelFontHeightF  = 0.015 
res@lbBoxEndCapStyle    = "TriangleBothEnds"
res@lbTitleOffsetF      =  0.1
res@lbTitleJust         ="CenterCenter"
res@lbTitlePosition     = "Bottom"
res@lbLabelOffsetF = 0.1 ;0.02
res@lbTitleFontHeightF  =0.018
res@pmLabelBarOrthogonalPosF= -0.05
res@mpMinLatF            = -58
res@mpMaxLatF            = 78
res@mpPerimOn       =False
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
res@gsnLeftStringParallelPosF= 0.01
res@gsnLeftStringOrthogonalPosF= -0.06
lres=True
lres@gsnDraw         = False                  ;         do not draw the plot                
lres@gsnFrame        = False     
lres@vpHeightF= res@vpHeightF                    ; change aspect ratio of plot
lres@vpWidthF = 0.4*res@vpWidthF                  
lres@xyLineThicknessF = 1.5
lres@tmXTOn        = False ; turn off top tickmarks
lres@tmYLOn        = False ; turn off right tickmarks
lres@tmXBLabelFontHeightF = 0.010
lres@tmYRLabelFontHeightF = 0.010
lres@tmXBLabelFont="times-roman"
lres@gsnLeftStringFontHeightF=0.018
lres@gsnLeftStringParallelPosF=  0.01
lres@gsnLeftStringOrthogonalPosF= -0.15
lres@gsnStringFont ="times-Bold"
lres@tmBorderThicknessF= 0.8
lres@tmYRLabelsOn =True
lres@tmYUseLeft =False
lres@tmYRLabelFont="times-roman"
lres@tmYRMode              = "Explicit"
lres@tmYRValues            = ispan(-60,80,20)
lres@tmYRLabels  = (/"60~S~o~N~S","40~S~o~N~S","20~S~o~N~S","0~S~o~N","20~S~o~N~N",\  
                     "40~S~o~N~N","60~S~o~N~N","80~S~o~N~N"/) 
lres@tiXAxisFont="times-roman"
lres@tiXAxisFontHeightF=0.012
lres@tmXBLabelFontHeightF = 0.014
lres@tmYRLabelFontHeightF = 0.014
lres@tmXBLabelFont="times-roman"
lres@tmXBLabelDeltaF = -0.7
lres@tmXBMajorThicknessF = 0.5
lres@tmYRMajorThicknessF = 0.5
lres@tmXBMajorOutwardLengthF = 0.0
lres@tmYRMajorOutwardLengthF=0.0
lres@tmXBMajorLengthF=0.005
lres@tmYRMajorLengthF=0.005

y = Gpool&lat   ;y(16::)  -60to90 80points
y_cesm = pool22&lat
printVarSummary(y)
printVarSummary(y_cesm)
crop_zonal=dim_avg_n_Wrap(Gpool, 2)
crop_zonal_cesm=dim_avg_n_Wrap(pool22, 2)
crop_zonal_ecearth=dim_avg_n_Wrap(bt2m_ecearth, 2)

printVarSummary(crop_zonal)
printVarSummary(crop_zonal_cesm)
printVarSummary(crop_zonal_ecearth)

frst_zonal=dim_avg_n_Wrap(Gpool2, 2)
frst_zonal_cesm=dim_avg_n_Wrap(flux2, 2)
frst_zonal_ecearth=dim_avg_n_Wrap(bmrsos_ecearth, 2)

irri_zonal=dim_avg_n_Wrap(Gpool4, 2)
irri_zonal_cesm=dim_avg_n_Wrap(pool23, 2)
irri_zonal_ecearth=dim_avg_n_Wrap(cor_ecearth, 2)
;...plot.............................
plot1=new((/12,3/),graphic)
plot2=new((/12,3/),graphic)
plot3=new((/12,3/),graphic)
res@tiXAxisFont="times-roman"
res@tiXAxisFontHeightF=0.01
res@tiYAxisFont="times-roman"
res@tiYAxisFontHeightF=0.012
res@gsnRightString    = " "
res@gsnStringFontHeightF = 0.012
res@tmXBLabelFontHeightF =0.01
res@tmYLLabelFontHeightF =0.01
res@gsnStringFont ="times-Bold"
res@tiMainFont="times-roman"
res@tiMainFontHeightF=0.019  ; 0.018
res@tiMainOffsetYF = -0.008
res@gsnLeftStringFontHeightF=0.019
res@gsnRightStringFontHeightF=0.011
res@gsnCenterStringFontHeightF=0.018
res@tmXBLabelFont         = "times-roman"
res@tmYLLabelFont         = "times-roman"   
res@cnMonoFillColor=False
res@cnMissingValFillColor="white" 
    do j=0,1
    res@tiMainString = "MPI-ESM"
    res@gsnLeftString    = "(a)"
res@cnLevels      =  (/-12.,-9.,-6.,-3.,0,3.,6.,9.,12./)
res@cnFillColors =   (/2,3,4,5,6,8,9,10,11,12/)
    plot3(0,j)=gsn_csm_contour_map(wks, Gpool_revised(j,:,:), res)
    res@tiYAxisString  = " "
    res@tiMainString = "CESM"
    res@gsnLeftString    = "(b)"
    res@lbTitleString  = "Sensitivity of vegetation C to surface temperature (10~S~-2~N~ kg/m~S~2~N~ per K)"

res@cnLevels  =  (/-3.2,-2.4,-1.6,-0.8,0,0.8,1.6,2.4,3.2/)
res@cnFillColors =   (/2,3,4,5,6,8,9,10,11,12/)
    plot3(1,j)=gsn_csm_contour_map(wks, pool2_revised(j,:,:), res)

    res@tiYAxisString  = " "
    res@tiMainString = "EC-Earth"
    res@gsnLeftString    = "(c)"
    res@lbTitleString  = " "
res@cnLevels      =  (/-12.,-9.,-6.,-3.,0,3.,6.,9.,12./)
res@cnFillColors =   (/2,3,4,5,6,8,9,10,11,12/)
printMinMax(bt2m_ecearth_revised(j,:,:), 1)
    plot3(2,j)=gsn_csm_contour_map(wks,bt2m_ecearth_revised(j,:,:), res)

    lres@xyLineColor= (/0.75,0.20,0.1,1./); "darkgoldenrod3"
    lres@gsnLeftString    = "(d)"
    lres@tiXAxisString=" "
    lres@trXMinF= -0.06
    lres@trXMaxF=0.08
    lres@tmXBMode              = "Explicit"
    lres@tmXBValues            = (/-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08/)
    lres@tmXBLabels  = (/" ","-4"," ",\  
                     "0"," ","4"," ","8"/)
    plot3(3,j)=gsn_csm_xy(wks,crop_zonal_ecearth(j,:),y_cesm,lres)
    res@tiMainString = " "
    res@gsnLeftString    = "(e)"
    res@cnLevels      =  (/-12.,-9.,-6.,-3.,0,3.,6.,9.,12./)
res@cnFillColors =   (/2,3,4,5,6,8,9,10,11,12/)

    plot3(4,j)=gsn_csm_contour_map(wks, Gpool2_revised(j,:,:), res)
    res@gsnCenterString = " "
    res@gsnLeftString    = "(f)"
    res@lbTitleString  = "Sensitivity of vegetation C to surface soil moisture (10~S~-2~N~ kg/m~S~2~N~ per kg/m~S~2~N~ water)"
    res@cnLevels  =  (/-3.2,-2.4,-1.6,-0.8,0,0.8,1.6,2.4,3.2/)
res@cnFillColors =   (/2,3,4,5,6,8,9,10,11,12/)

    plot3(5,j)=gsn_csm_contour_map(wks, flux2_revised(j,:,:), res)

    res@gsnCenterString = " "
    res@gsnLeftString    = "(g)"
    res@lbTitleString  = " "
    res@cnLevels      =  (/-48,-36,-24,-12,0,12,24,36,48/)
res@cnFillColors =   (/2,3,4,5,6,8,9,10,11,12/)

printMinMax(bmrsos_ecearth_revised(j,:,:), 1)
    plot3(6,j)=gsn_csm_contour_map(wks, bmrsos_ecearth_revised(j,:,:), res)

    lres@gsnLeftString    = "(h)"
    lres@tiXAxisString=" "


    lres@trXMinF= -0.16
    lres@trXMaxF=0.48
    lres@tmXBMode              = "Explicit"
    lres@tmXBValues            = (/-0.08,0,0.08,0.16,0.24,0.32,0.40,0.48/)
    lres@tmXBLabels  = (/" ","0"," ","16"," ","32"," ","48"/)  ;  lres@tiXAxisString=" "

    plot3(7,j)=gsn_csm_xy(wks,frst_zonal_ecearth(j,:),y_cesm,lres)

res@gsnCenterString = " "
res@lbLabelBarOn         = False ; Turn off individual labelbar;
res@gsnLeftString    = "(i)"
res@cnLevels      =  (/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9/)
res@cnFillColors =   (/0,7,14,15,16,17,18,19,20,21/)
    plot3(8,j)=gsn_csm_contour_map(wks, Gpool4_revised(j,:,:), res)
    res@tiYAxisString  = " "
    res@gsnCenterString = " "
    res@gsnLeftString    = "(j)"
    res@lbLabelBarOn         = True ; Turn off individual labelbar;
    res@pmLabelBarHeightF   = 0.05
    res@pmLabelBarWidthF = 0.8
res@lbTitleOffsetF      =0.1
res@lbLabelOffsetF = 0.1 ;0.02
res@lbTitleFontHeightF  =0.018
res@pmLabelBarOrthogonalPosF= 0.01
    res@lbTitleString= "R-squared"
    plot3(9,j)=gsn_csm_contour_map(wks, pool23_revised(j,:,:), res)

    res@tiYAxisString  = " "
    res@gsnCenterString = " "
    res@gsnLeftString    = "(k)"
    res@lbLabelBarOn         = False ; Turn off individual labelbar;
    res@lbTitleString= " "
    res@tiXAxisString= " "
    plot3(10,j)=gsn_csm_contour_map(wks, cor_ecearth_revised(j,:,:), res)

    lres@gsnLeftString    = "(l)"
    lres@trXMinF= 0
    lres@trXMaxF=0.8
    lres@tmXBMode              = "Explicit"
    lres@tmXBValues            = (/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8/)
    lres@tmXBLabels  = (/" ","0.2"," ","0.4"," ","0.6"," ","0.8"/)   
    lres@tiXAxisString=" "
    lres@lbTitleString= " "
    lres@tiXAxisFont="times-roman"
    lres@tiXAxisFontHeightF=0.012
    plot3(11,j)=gsn_csm_xy(wks,irri_zonal_ecearth(j,:),y_cesm,lres)

    lres@tiXAxisString=" "
   ;................................
    lres@xyLineColor= (/0.0,0.0,0.0,1.0/)  ; "black" cesm
    res@gsnLeftString    = " "
    res@gsnCenterString = " "
    res@gsnRightString=" "  
    res@cnMonoFillColor=True
    res@cnFillColor="Transparent"
    res@cnMissingValFillColor="white"
    res@lbLabelBarOn         = False ; Turn off individual labelbar;
    lres@gsnLeftString    = " "
    plot1(0,j)=gsn_csm_contour(wks, Gpool(j,:,:), res)
    plot1(1,j)=gsn_csm_contour(wks, pool22(j,:,:), res)
    plot1(2,j)=gsn_csm_contour(wks, bt2m_ecearth(j,:,:), res)
    plot1(3,j)=gsn_csm_xy(wks,crop_zonal_cesm(j,:),y_cesm,lres)

    plot1(4,j)=gsn_csm_contour(wks, Gpool2(j,:,:), res)
    plot1(5,j)=gsn_csm_contour(wks, flux2(j,:,:), res)
    plot1(6,j)=gsn_csm_contour(wks, bmrsos_ecearth(j,:,:), res)
    plot1(7,j)=gsn_csm_xy(wks,frst_zonal_cesm(j,:),y_cesm,lres)

    plot1(8,j)=gsn_csm_contour(wks, Gpool4(j,:,:), res)
    plot1(9,j)=gsn_csm_contour(wks, pool23(j,:,:), res)
    plot1(10,j)=gsn_csm_contour(wks, cor_ecearth(j,:,:), res)
    plot1(11,j)=gsn_csm_xy(wks,irri_zonal_cesm(j,:),y_cesm,lres)

    lres@xyLineColor= (/0.11,0.44,0.93,1./) ;green3 ;MPI-ESM
    plot2(0,j)=gsn_csm_contour(wks, Gpool(j,:,:), res)
    plot2(1,j)=gsn_csm_contour(wks, pool22(j,:,:), res)
    plot2(2,j)=gsn_csm_contour(wks, bt2m_ecearth(j,:,:), res)
    plot2(3,j)=gsn_csm_xy(wks,crop_zonal(j,:),y,lres)

    plot2(4,j)=gsn_csm_contour(wks, Gpool2(j,:,:), res)
    plot2(5,j)=gsn_csm_contour(wks, flux2(j,:,:), res)
    plot2(6,j)=gsn_csm_contour(wks, bmrsos_ecearth(j,:,:), res)
    plot2(7,j)=gsn_csm_xy(wks,frst_zonal(j,:),y,lres)

    plot2(8,j)=gsn_csm_contour(wks, Gpool4(j,:,:), res)
    plot2(9,j)=gsn_csm_contour(wks, pool23(j,:,:), res)
    plot2(10,j)=gsn_csm_contour(wks, cor_ecearth(j,:,:), res)
    plot2(11,j)=gsn_csm_xy(wks,irri_zonal(j,:),y,lres)

end do 
do i=0,11
    do j=0,1
       overlay(plot3(i,j),plot1(i,j))
       overlay(plot3(i,j),plot2(i,j))
end do
end do

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
lgres@lgLineColors       = (/(/0.75,0.20,0.1,1./),(/0.0,0.0,0.0,1./),(/0.11,0.44,0.93,1./)/)
gsn_legend_ndc(wks,3,(/"  EC-Earth","  CESM","  MPI-ESM"/),0.77,0.302,lgres)


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
      ;  plot7(i)=gsn_add_polyline(wks, plot3(0,0), plx(i,:), ply(i,:), plres)
end do
;end
pres                  = True         ; Set panel resources.
pres@gsnMaximize      = True         ; Maximize plots in panel.
pres@gsnPanelLabelBar = False       ; Turn on panel labelbar.
pres@pmLabelBarOrthogonalPosF=0.027
pres@pmLabelBarParallelPosF=-0.15
pres@pmLabelBarWidthF    = 0.6      ; default is shorter
pres@lbPerimOn             = False
pres@lbBoxLineColor        = "white"
pres@lbBoxEndCapStyle    = "TriangleBothEnds"
pres@lbTitleFontHeightF  =0.008
pres@lbTitleFont         ="times-roman"
pres@lbTitleString       =" "
pres@lbTitleOffsetF      =-0.1
pres@lbTitleJust         ="CenterCenter"
pres@lbTitlePosition     = "Bottom"
pres@tiYAxisString=""
pres@txFontHeightF  =0.005
pres@gsnPanelXF = (/0.0185,0.265,0.511,0.759,0.0185,0.265,0.511,0.759,0.0185,0.265,0.511,0.759/) ;0.4
gsn_panel(wks,plot3(:,0),(/3,4/),pres)
end
