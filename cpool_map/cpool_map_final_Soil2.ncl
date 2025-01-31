load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin
 wks2 = gsn_open_wks("pdf","testBrYLGr")
  gsn_define_colormap(wks2,"BrYlGr")
  gsn_draw_colormap(wks2)
;......................................................................
myfile1 = addfile("/work/bm1147/b380949/web-monitoring/crop-ctl/Eyr/cVeg/cVeg_crop-ctl_mpiesm_signal-separated_revised.nc", "r")
myfile2 = addfile("/work/bm1147/b380949/web-monitoring/crop-ctl/Eyr/cSoil/cSoil_crop-ctl_mpiesm_signal-separated_revised.nc", "r")
myfile3 = addfile("/work/bm1147/b380949/web-monitoring/crop-ctl/Eyr/cProduct/cProduct_crop-ctl_mpiesm_signal-separated_revised.nc", "r")
myfile4= addfile("/work/bm1147/b380949/web-monitoring/crop-ctl/Eyr/cLitter/cLitter_crop-ctl_mpiesm_signal-separated_revised.nc", "r")

myfile5 = addfile("/work/bm1147/b380949/web-monitoring/frst-ctl/Eyr/cVeg/cVeg_frst-ctl_mpiesm_signal-separated_revised.nc", "r")
myfile6 = addfile("/work/bm1147/b380949/web-monitoring/frst-ctl/Eyr/cSoil/cSoil_frst-ctl_mpiesm_signal-separated_revised.nc", "r")
myfile7 = addfile("/work/bm1147/b380949/web-monitoring/frst-ctl/Eyr/cProduct/cProduct_frst-ctl_mpiesm_signal-separated_revised.nc", "r")
myfile8 = addfile("/work/bm1147/b380949/web-monitoring/frst-ctl/Eyr/cLitter/cLitter_frst-ctl_mpiesm_signal-separated_revised.nc", "r")

myfile13 = addfile("/work/bm1147/b380949/web-monitoring/irri-crop/Eyr/cVeg/cVeg_irri-crop_mpiesm_signal-separated_revised.nc", "r")
myfile14 = addfile("/work/bm1147/b380949/web-monitoring/irri-crop/Eyr/cSoil/cSoil_irri-crop_mpiesm_signal-separated_revised.nc", "r")
myfile15 = addfile("/work/bm1147/b380949/web-monitoring/irri-crop/Eyr/cProduct/cProduct_irri-crop_mpiesm_signal-separated_revised.nc", "r")
myfile16 = addfile("/work/bm1147/b380949/web-monitoring/irri-crop/Eyr/cLitter/cLitter_irri-crop_mpiesm_signal-separated_revised.nc", "r")

ti=cd_calendar(myfile1->time,-1)
ti2=cd_calendar(myfile5->time,-1)
ti4=cd_calendar(myfile13->time,-1)
;........mpiesm...........
;.........deforestation..... ..................
pool=new((/5,162,96,192/),double) 
pool(0,:,:,:) = myfile1->cVeg_nonlocal(:,::-1,:)
pool&time=ti
pool(1,:,:,:) =myfile2->cSoil_nonlocal(:,::-1,:)
pool(2,:,:,:) =myfile3->cProduct_nonlocal(:,::-1,:)
pool(3,:,:,:) =myfile4->cLitter_nonlocal(:,::-1,:)
pool(4,:,:,:)=dim_sum_n_Wrap(pool(0:3,:,:,:),0)
;..............afforestation...................
pool2=new((/5,161,96,192/),double) 
pool2(0,:,:,:) = myfile5->cVeg_nonlocal(:,::-1,:)
pool2&time=ti2
pool2(1,:,:,:) =myfile6->cSoil_nonlocal(:,::-1,:)
pool2(2,:,:,:) =myfile7->cProduct_nonlocal(:,::-1,:)
pool2(3,:,:,:) =myfile8->cLitter_nonlocal(:,::-1,:)
pool2(4,:,:,:)=dim_sum_n_Wrap(pool2(0:3,:,:,:),0)
;.....................irrigation......................................
pool4=new((/5,160,96,192/),double) 
pool4(0,:,:,:) = myfile13->cVeg_nonlocal(:,::-1,:)
pool4&time=ti4
pool4(1,:,:,:) =myfile14->cSoil_nonlocal(:,::-1,:)
pool4(2,:,:,:) =myfile15->cProduct_nonlocal(:,::-1,:)
pool4(3,:,:,:) =myfile16->cLitter_nonlocal(:,::-1,:)
pool4(4,:,:,:)=dim_sum_n_Wrap(pool4(0:3,:,:,:),0)

;.............last 30 year ave................................
Gpool=new((/5,96,192/),double) 
Gpool2=new((/5,96,192/),double) 
Gpool4=new((/5,96,192/),double) 
do i=0,3
Gpool(i,:,:)=dim_avg_n_Wrap(pool(i,130:159,:,:),0)  ;kg/m^2 ；last30years
Gpool2(i,:,:)=dim_avg_n_Wrap(pool2(i,130:159,:,:),0)  ;kg/m^2
Gpool4(i,:,:)=dim_avg_n_Wrap(pool4(i,130:159,:,:),0)  ;kg/m^2
end do
printVarSummary(Gpool)
Gpool(4,:,:)=dim_sum_n_Wrap(Gpool(0:3,:,:),0)
Gpool2(4,:,:)=dim_sum_n_Wrap(Gpool2(0:3,:,:),0)
Gpool4(4,:,:)=dim_sum_n_Wrap(Gpool4(0:3,:,:),0)
;................white for diff=0........................
Gpool_revised=mask(Gpool, Gpool.eq.0, False)   ;crop(5,96,192)5:veg,soil,product,litter,land
copy_VarMeta(Gpool, Gpool_revised)
Gpool2_revised=mask(Gpool2, Gpool2.eq.0, False)  ;frst
copy_VarMeta(Gpool2, Gpool2_revised)
Gpool4_revised=mask(Gpool4, Gpool4.eq.0, False)  ;irri
copy_VarMeta(Gpool4, Gpool4_revised)
;....CESM.................
myfile12 = addfile("/work/bm1147/b380949/web-monitoring/CESM/crop-ctl/Lmon/TOTVEGC/TOTVEGC_crop-ctl_cesm_signal-separated.nc", "r")
myfile22 = addfile("/work/bm1147/b380949/web-monitoring/CESM/crop-ctl/Lmon/TOTSOMC/TOTSOMC_crop-ctl_cesm_signal-separated.nc", "r")
myfile42 = addfile("/work/bm1147/b380949/web-monitoring/CESM/crop-ctl/Lmon/TOTLITC/TOTLITC_crop-ctl_cesm_signal-separated.nc", "r")
myfile32 = addfile("/work/bm1147/b380949/web-monitoring/CESM/crop-ctl/Lmon/TOTECOSYSC/TOTECOSYSC_crop-ctl_cesm_signal-separated.nc", "r")

myfile52 = addfile("/work/bm1147/b380949/web-monitoring/CESM/frst-ctl/Lmon/TOTVEGC/TOTVEGC_frst-ctl_cesm_signal-separated.nc", "r")
myfile62 = addfile("/work/bm1147/b380949/web-monitoring/CESM/frst-ctl/Lmon/TOTSOMC/TOTSOMC_frst-ctl_cesm_signal-separated.nc", "r")
myfile82 = addfile("/work/bm1147/b380949/web-monitoring/CESM/frst-ctl/Lmon/TOTLITC/TOTLITC_frst-ctl_cesm_signal-separated.nc", "r")
myfile72 = addfile("/work/bm1147/b380949/web-monitoring/CESM/frst-ctl/Lmon/TOTECOSYSC/TOTECOSYSC_frst-ctl_cesm_signal-separated.nc", "r")

myfile132 = addfile("/work/bm1147/b380949/web-monitoring/CESM/irri-crop/Lmon/TOTVEGC/TOTVEGC_irri-crop_cesm_signal-separated.nc", "r")
myfile142 = addfile("/work/bm1147/b380949/web-monitoring/CESM/irri-crop/Lmon/TOTSOMC/TOTSOMC_irri-crop_cesm_signal-separated.nc", "r")
myfile152 = addfile("/work/bm1147/b380949/web-monitoring/CESM/irri-crop/Lmon/TOTLITC/TOTLITC_irri-crop_cesm_signal-separated.nc", "r")
myfile162 = addfile("/work/bm1147/b380949/web-monitoring/CESM/irri-crop/Lmon/TOTECOSYSC/TOTECOSYSC_irri-crop_cesm_signal-separated.nc", "r")

myfile=addfile("/home/b/b380949/gridarea/landmask_glacier_cesm.nc", "r")
slm = myfile->landmask
;.........deforestation..... ..................
pool2n=new((/4,192,288/),double) ;crop
pool2n(0,:,:) = dim_avg_n_Wrap(myfile12->TOTVEGC_nonlocal(1560:1919,:,:),0)
pool2n(1,:,:) =dim_avg_n_Wrap(myfile22->TOTSOMC_nonlocal(1560:1919,:,:),0)
pool2n(2,:,:) =dim_avg_n_Wrap(myfile42->TOTLITC_nonlocal(1560:1919,:,:),0)
pool2n(3,:,:) =dim_avg_n_Wrap(myfile32->TOTECOSYSC_nonlocal(1560:1919,:,:),0)

;.............afforestation..........
flux2n=new((/4,192,288/),double)  ;frst
flux2n(0,:,:) = dim_avg_n_Wrap(myfile52->TOTVEGC_nonlocal(1560:1919,:,:),0)
flux2n(1,:,:) =dim_avg_n_Wrap(myfile62->TOTSOMC_nonlocal(1560:1919,:,:),0)
flux2n(2,:,:) =dim_avg_n_Wrap(myfile82->TOTLITC_nonlocal(1560:1919,:,:),0)
flux2n(3,:,:) =dim_avg_n_Wrap(myfile72->TOTECOSYSC_nonlocal(1560:1919,:,:),0)

;print(flx2(0,0,159,:,:))
;................irri...................
pool23n=new((/4,192,288/),double)  ;irrigation
pool23n(0,:,:)= dim_avg_n_Wrap(myfile132->TOTVEGC_nonlocal(1560:1919,:,:),0)
pool23n(1,:,:)=dim_avg_n_Wrap(myfile142->TOTSOMC_nonlocal(1560:1919,:,:),0)
pool23n(2,:,:)= dim_avg_n_Wrap(myfile152->TOTLITC_nonlocal(1560:1919,:,:),0)
pool23n(3,:,:)= dim_avg_n_Wrap(myfile162->TOTECOSYSC_nonlocal(1560:1919,:,:),0)
pool2m=pool2n*0.001
flux2m=flux2n*0.001
pool23m=pool23n*0.001
copy_VarMeta(pool2n, pool2m)
copy_VarMeta(flux2n, flux2m)
copy_VarMeta(pool23n, pool23m)
;.......Deal with missing values............
pool22=mask(pool2m,  slm.eq.0, False)
copy_VarMeta(pool2n, pool22)

pool23=mask(pool23m, slm.eq.0, False)
copy_VarMeta(pool23n, pool23)

flux2=mask(flux2m, slm.eq.0, False)
copy_VarMeta(flux2n, flux2)
;...................
pool2_revised=mask(pool22, pool22.eq.0, False);crop;4,192,288;4:veg,soil,litter,land
copy_VarMeta(pool2n, pool2_revised)

pool23_revised=mask(pool23, pool23.eq.0, False) ;irri
copy_VarMeta(pool23n, pool23_revised)

flux2_revised=mask(flux2, flux2.eq.0, False);frst
copy_VarMeta(flux2n, flux2_revised)
;............ec-earth.........................
landsea_ecearth=addfile("/work/bm1147/b380949/EC-earth/dataset/anlysisneed/landmask_ecearth_regridded.nc","r")
landmask_ecearth=landsea_ecearth->landmask(0,:,:)
myfile13n = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/crop-ctl/Eyr/cVeg/interped_cVeg_crop-ctl_ecearth_signal-separated.nc", "r")
myfile23 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/crop-ctl/Eyr/cSoil/interped_cSoil_crop-ctl_ecearth_signal-separated.nc", "r")
myfile43 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/crop-ctl/Eyr/cLitter/interped_cLitter_crop-ctl_ecearth_signal-separated.nc", "r")
myfile33 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/crop-ctl/Eyr/cLand/interped_cLand_crop-ctl_ecearth_signal-separated.nc", "r")

myfile53 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/frst-ctl/Eyr/cVeg/interped_cVeg_frst-ctl_ecearth_signal-separated.nc", "r")
myfile63 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/frst-ctl/Eyr/cSoil/interped_cSoil_frst-ctl_ecearth_signal-separated.nc", "r")
myfile83 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/frst-ctl/Eyr/cLitter/interped_cLitter_frst-ctl_ecearth_signal-separated.nc", "r")
myfile73 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/frst-ctl/Eyr/cLand/interped_cLand_frst-ctl_ecearth_signal-separated.nc", "r")

myfile133 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/irr-crop/Eyr/cVeg/interped_cVeg_irr-crop_ecearth_signal-separated.nc", "r")
myfile143 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/irr-crop/Eyr/cSoil/interped_cSoil_irr-crop_ecearth_signal-separated.nc", "r")
myfile153 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/irr-crop/Eyr/cLitter/interped_cLitter_irr-crop_ecearth_signal-separated.nc", "r")
myfile163 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/irr-crop/Eyr/cLand/interped_cLand_irr-crop_ecearth_signal-separated.nc", "r")

pool4eo=new((/4,192,288/),double) 
pool4eo(0,:,:) = dim_avg_n_Wrap(myfile13n->cVeg_nonlocal(130:159,:,:),0)
pool4eo(1,:,:) =dim_avg_n_Wrap(myfile23->cSoil_nonlocal(130:159,:,:),0)
pool4eo(2,:,:) =dim_avg_n_Wrap(myfile43->cLitter_nonlocal(130:159,:,:),0)
pool4eo(3,:,:) =dim_avg_n_Wrap(myfile33->cLand_nonlocal(130:159,:,:),0)
;printVarSummary(pool4eo)
;.............afforestation..........
flux4eo=new((/4,192,288/),double) 
flux4eo(0,:,:) = dim_avg_n_Wrap(myfile53->cVeg_nonlocal(130:159,:,:),0)
flux4eo(1,:,:) =dim_avg_n_Wrap(myfile63->cSoil_nonlocal(130:159,:,:),0)
flux4eo(2,:,:) =dim_avg_n_Wrap(myfile83->cLitter_nonlocal(130:159,:,:),0)
flux4eo(3,:,:) =dim_avg_n_Wrap(myfile73->cLand_nonlocal(130:159,:,:),0)
;print(flux4(0,0,159,:,:))
;................irri...................
pool43eo=new((/4,192,288/),double) 
pool43eo(0,:,:)= dim_avg_n_Wrap(myfile133->cVeg_nonlocal(130:159,:,:),0)
pool43eo(1,:,:)=dim_avg_n_Wrap(myfile143->cSoil_nonlocal(130:159,:,:),0)
pool43eo(2,:,:)= dim_avg_n_Wrap(myfile153->cLitter_nonlocal(130:159,:,:),0)
pool43eo(3,:,:) =dim_avg_n_Wrap(myfile163->cLand_nonlocal(130:159,:,:),0)
;.......Deal with missing values.................
pool4eo=mask(pool4eo,isnan_ieee(pool4eo),False)
flux4eo=mask(flux4eo,isnan_ieee(flux4eo),False)
pool43eo=mask(pool43eo,isnan_ieee(pool43eo),False)
print(any(isnan_ieee(pool4eo)))
print(any(isnan_ieee(flux4eo)))
print(any(isnan_ieee(pool43eo)))
pool4e=mask(pool4eo, landmask_ecearth.eq.1, True)
copy_VarMeta(pool4eo, pool4e)
flux4e=mask(flux4eo, landmask_ecearth.eq.1, True)
copy_VarMeta(flux4eo, flux4e)
pool43e=mask(pool43eo, landmask_ecearth.eq.1, True)
copy_VarMeta(pool43eo, pool43e)
;.................................
pool4e_revised=mask(pool4e, pool4e.eq.0, False);crop;4,192,288;4:veg,soil,litter,land
copy_VarMeta(pool4e, pool4e_revised)

flux4e_revised=mask(flux4e, flux4e.eq.0, False) ;irri
copy_VarMeta(flux4e, flux4e_revised)

pool43e_revised=mask(pool43e, pool43e.eq.0, False);frst
copy_VarMeta(pool43e, pool43e_revised)
;..............plot.................................................
wks_type = "pdf"
wks_type@wkPaperWidthF =  9.0
wks_type@wkPaperHeightF = 9.0
wks_type@wkPDFResolution = 3000
wks = gsn_open_wks(wks_type,"cpool_map_final_cSoil_16Oct")
gsn_define_colormap(wks,"BrYlGr")
;gsn_reverse_colormap(wks)
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
res@cnLevelSelectionMode = "ExplicitLevels"
res@cnLevels  =  (/-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8/)
res@cnFillColors =   (/2,3,4,5,6,7,8,9,10,11,12/)
res@cnLineThicknessF = 1.0
res@cnLineLabelsOn = False
res@tmYLMinorOn           = True
res@lbLabelBarOn         = False ; Turn off individual labelbar;
res@lbBoxEndCapStyle    = "TriangleBothEnds"
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
lres@tmXBLabelAngleF = 15
lres@tmYRLabelFontHeightF = 0.014
lres@tmXBLabelFont="times-roman"
lres@tmXBLabelDeltaF = -0.7
lres@gsnLeftStringParallelPosF= 0.04
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
printVarSummary(Gpool4(0,{-60:90},:))
printVarSummary(pool22(0,{-60:90},:))
printVarSummary(pool4e(0,{-60:90},:))
y = Gpool&lat   ;y(16::)  -60to90 80points
y_cesm = pool22&lat
y_ecearth=pool4e&lat

crop_zonal=dim_avg_n_Wrap(Gpool(:,{-60:90},:), 2)
crop_zonal_cesm=dim_avg_n_Wrap(pool22(:,{-60:90},:), 2)
crop_zonal_ecearth=dim_avg_n_Wrap(pool4e(:,{-60:90},:), 2)

frst_zonal=dim_avg_n_Wrap(Gpool2(:,{-60:90},:), 2)
frst_zonal_cesm=dim_avg_n_Wrap(flux2(:,{-60:90},:), 2)
frst_zonal_ecearth=dim_avg_n_Wrap(flux4e(:,{-60:90},:), 2)

irri_zonal=dim_avg_n_Wrap(Gpool4(:,{-60:90},:), 2)
irri_zonal_cesm=dim_avg_n_Wrap(pool23(:,{-60:90},:), 2)
irri_zonal_ecearth=dim_avg_n_Wrap(pool43e(:,{-60:90},:), 2)
;...plot.............................
plot1=new((/12,3/),graphic)
plot3=new((/12,3/),graphic)
plot2=new((/12,3/),graphic)

str=(/"cVeg","cSoil","cProduct","cLitter","cLand"/)
str2=(/"_Total","_Local","_Nonlocal"/)
res@tiXAxisFont="times-roman"
res@tiXAxisFontHeightF=0.01
res@tiYAxisFont="times-roman"
res@tiMainFont="times-roman"
res@tiMainFontHeightF= 0.019
res@tiMainOffsetYF = -0.008
res@tiYAxisFontHeightF=0.02
res@gsnLeftStringParallelPosF= 0.04
res@gsnLeftStringOrthogonalPosF=-0.11
 res@gsnRightString    = " "
res@tmXBLabelFontHeightF =0.014
res@tmYLLabelFontHeightF =0.014
res@gsnStringFont ="times-Bold"
res@gsnLeftStringFontHeightF=0.019
res@gsnRightStringFontHeightF=0.011
res@gsnCenterStringFontHeightF=0.018
res@tmXBLabelFont         = "times-roman"
res@tmYLLabelFont         = "times-roman"    
    do j=0,2
  lres@xyLineColor= (/0.75,0.20,0.1,1./) ;red 
   res@cnMonoFillColor=False
    res@cnMissingValFillColor="white"
     res@tiYAxisString  = "CROP-CTL"
     res@tiMainString ="MPI-ESM"
    res@gsnLeftString    = "(a)"
    plot3(0,j)=gsn_csm_contour_map(wks, Gpool_revised(j,:,:), res)
    res@tiYAxisString  = " "
    res@tiMainString ="CESM"
    res@gsnLeftString    = "(b)"
    plot3(1,j)=gsn_csm_contour_map(wks, pool2_revised(j,:,:), res)
    res@tiMainString ="EC-Earth"
    res@gsnLeftString    = "(c)"
    plot3(2,j)=gsn_csm_contour_map(wks, pool4e_revised(j,:,:), res)
    res@gsnCenterString = ""
    res@tiMainString =""

    lres@gsnStringFont ="times-Bold"
    lres@gsnLeftStringFontHeightF=0.019
    lres@gsnLeftString    = "(d)"
    lres@trXMinF= -0.4
    lres@trXMaxF=max(crop_zonal_ecearth(1,:))
    lres@tmXBMode              = "Explicit"
    lres@tmXBValues            = (/-0.4,-0.2,0,0.2,0.4/)
    lres@tmXBLabels  = (/"-0.4","-0.2","0","0.2","0.4"/)

plest                   = True                      ; polyline mods desired
  plest@gsLineColor       = "grey"                    ; color of lines
  plest@gsLineThicknessF  = 1.2                       ; thickness of lines
  plest@gsLineOpacityF = 0.8
  dum=new((/3,2/),"graphic")
  dum(0,0)=gsn_add_polyline(wks,plot3(3,0),(/0,0/),(/-60,90/),plest)
  dum(1,0)=gsn_add_polyline(wks,plot3(7,0),(/0,0/),(/-60,90/),plest)
  dum(2,0)=gsn_add_polyline(wks,plot3(11,0),(/0,0/),(/-60,90/),plest)
  dum(0,1)=gsn_add_polyline(wks,plot3(3,1),(/0,0/),(/-60,90/),plest)
  dum(1,1)=gsn_add_polyline(wks,plot3(7,1),(/0,0/),(/-60,90/),plest)
  dum(2,1)=gsn_add_polyline(wks,plot3(11,1),(/0,0/),(/-60,90/),plest)
    plot3(3,j)=gsn_csm_xy(wks,crop_zonal_ecearth(j,:),y_ecearth(32::),lres) ;ecearth
   res@tiYAxisString  = "FRST-CTL"
    res@gsnLeftString    = "(e)"
    plot3(4,j)=gsn_csm_contour_map(wks, Gpool2_revised(j,:,:), res)
        res@tiYAxisString  = " "
    res@gsnLeftString    = "(f)"
    plot3(5,j)=gsn_csm_contour_map(wks, flux2_revised(j,:,:), res)
    res@gsnLeftString    = "(g)"
    plot3(6,j)=gsn_csm_contour_map(wks, flux4e_revised(j,:,:), res)

    lres@gsnLeftString    = "(h)"
    plot3(7,j)=gsn_csm_xy(wks,frst_zonal_ecearth(j,:),y_ecearth(32::),lres) ;ecearth

   res@tiYAxisString  = "IRR-CROP"

    res@gsnLeftString    = "(i)"
    plot3(8,j)=gsn_csm_contour_map(wks, Gpool4_revised(j,:,:), res)
        res@tiYAxisString  = " "
    res@gsnLeftString    = "(j)"
    plot3(9,j)=gsn_csm_contour_map(wks, pool23_revised(j,:,:), res)

    res@gsnLeftString    = "(k)"
    plot3(10,j)=gsn_csm_contour_map(wks, pool43e_revised(j,:,:), res)
    lres@gsnLeftString    = "(l)"

    lres@tiXAxisString= "";"kg/m~S~2~N~"
    lres@tiXAxisFont="times-roman"
    lres@tiXAxisFontHeightF=0.018
    plot3(11,j)=gsn_csm_xy(wks,irri_zonal_ecearth(j,:),y_ecearth(32::),lres) ;ecearth
    lres@tiXAxisString=" "

   ;................................
    lres@xyLineColor= (/0.0,0.0,0.0,1.0/)  ; "black"
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
    plot1(2,j)=gsn_csm_contour(wks, pool4e(j,:,:), res)


    plot2(0,j)=gsn_csm_contour(wks, Gpool(j,:,:), res)
    plot2(1,j)=gsn_csm_contour(wks, pool22(j,:,:), res)
    plot2(2,j)=gsn_csm_contour(wks, pool4e(j,:,:), res)

    plot1(3,j)=gsn_csm_xy(wks,crop_zonal_cesm(j,:),y_cesm(32::),lres)
        lres@xyLineColor= (/0.11,0.44,0.93,1./) ;blue 
    printMinMax(crop_zonal_ecearth(0,:), 0)
    plot2(3,j)=gsn_csm_xy(wks,crop_zonal(j,:),y(16::),lres) ;MPI-ESM
    lres@xyLineColor= (/0.0,0.0,0.0,1./)
    plot1(4,j)=gsn_csm_contour(wks, Gpool2(j,:,:), res)
    plot1(5,j)=gsn_csm_contour(wks, flux2(j,:,:), res)
    plot1(6,j)=gsn_csm_contour(wks, flux4e(j,:,:), res)

    plot2(4,j)=gsn_csm_contour(wks, Gpool2(j,:,:), res)
    plot2(5,j)=gsn_csm_contour(wks, flux2(j,:,:), res)
    plot2(6,j)=gsn_csm_contour(wks, flux4e(j,:,:), res)
    plot1(7,j)=gsn_csm_xy(wks,frst_zonal_cesm(j,:),y_cesm(32::),lres)
        lres@xyLineColor= (/0.11,0.44,0.93,1./) ;blue
        printMinMax(frst_zonal_ecearth(0,:),0)
    plot2(7,j)=gsn_csm_xy(wks,frst_zonal(j,:),y(16::),lres) ;MPI-ESM
   
    plot1(8,j)=gsn_csm_contour(wks, Gpool4(j,:,:), res)
    plot1(9,j)=gsn_csm_contour(wks, pool23(j,:,:), res)
    plot1(10,j)=gsn_csm_contour(wks, pool43e(j,:,:), res)

    plot2(8,j)=gsn_csm_contour(wks, Gpool4(j,:,:), res)
    plot2(9,j)=gsn_csm_contour(wks, pool23(j,:,:), res)
    plot2(10,j)=gsn_csm_contour(wks, pool43e(j,:,:), res)
     lres@xyLineColor= (/0.0,0.0,0.0,1./)
    plot1(11,j)=gsn_csm_xy(wks,irri_zonal_cesm(j,:),y_cesm(32::),lres)
        lres@xyLineColor= (/0.11,0.44,0.93,1./) ;blue
    printMinMax(irri_zonal_ecearth(0,:), 0)
    plot2(11,j)=gsn_csm_xy(wks,irri_zonal(j,:),y(16::),lres) ;MPI-ESM
end do 

do i=0,11
    do j=0,2
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
strk=(/" "," "," "/)
strk1=(/"  MPI-ESM","  CESM"," EC-Earth"/)

colors=(/"royalblue","orangered4","gold3"/)
lgres@lgLineColors       = (/(/0.75,0.20,0.1,1./),(/0.0,0.0,0.0,1./),(/0.11,0.44,0.93,1./)/)

gsn_legend_ndc(wks,3,(/"  EC-Earth","  CESM","  MPI-ESM"/),0.75,0.33,lgres)

  plres = True
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
    plres@gsLineColor            = (/0.75,0.34,0.37,0.8/)
    plres@gsLineThicknessF = 1.5
        plot7(i)=gsn_add_polyline(wks, plot3(0,1), plx(i,:), ply(i,:), plres)
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
pres@lbBoxEndCapStyle    = "TriangleBothEnds"
pres@lbTitleFontHeightF  =0.01
pres@lbTitleFont         ="times-roman"
pres@lbTitleString       ="Nonlocal soil carbon change (kg/m~S~2~N~)  "
pres@lbTitleOffsetF      = 0.2 ;-0.2
pres@lbTitleJust         ="CenterCenter"
pres@lbTitlePosition     = "Bottom"
pres@gsnPanelYWhiteSpacePercent = 0.005
pres@gsnPanelXWhiteSpacePercent =0.0
pres@tiYAxisString=""
pres@txFontHeightF  =0.018
pres@txFont ="times-roman"
pres@gsnPanelXF = (/0.0385,0.275,0.512,0.749,0.0385,0.275,0.512,0.749,0.0385,0.275,0.512,0.749/) ;0.4
pres@gsnPanelYF = (/0.664,0.664,0.664,0.664,0.554,0.554,0.554,0.554,0.444,0.444,0.444,0.444/) ;0.4
pres@gsnPanelDebug = True
gsn_panel(wks,plot3(:,1),(/3,4/),pres)
end
