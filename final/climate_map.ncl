load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin
;...........Climate map uses 150-year mean corresponding to the ensemble mean published in De Hertog et al., 2023.
;.............from time =11 to time =160
;..............cesm...crop-ctl.............................
myfile = addfile("/pf/b/b380949/gridarea/gridarea_sea0.nc","r")
my=addfile("/pf/b/b380949/gridarea/AREA_CESM.nc","r")
myfile1 = addfile("/work/bm1147/b380949/web-monitoring/CESM/crop-ctl/Lmon/SOILWATER_10CM/SOILWATER_10CM_crop-ctl_cesm_signal-separated.nc", "r")
;;kg m-2 mrsos
;(time,depth,lat,lon) kg m-2 mrsol
myfile2 = addfile("/work/bm1147/b380949/web-monitoring/CESM/crop-ctl/Lmon/TREFHT/TREFHT_crop-ctl_cesm_signal-separated.nc", "r")
tas_crop_cesm=new((/192,288/),double) ;3varia, 3 simulation
mrsos_crop_cesm=new((/192,288/),double) ;
tas_crop_cesm = dim_avg_n_Wrap(myfile2->TREFHT_nonlocal,0)
mrsos_crop_cesm = dim_avg_n_Wrap(myfile1->SOILWATER_10CM_nonlocal(120:1919,:,:),0)
;......ecearth.....crop-ctl....150 time steps.......from time= 11 to time= 160.......................
myfile3 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/crop-ctl/Eyr/mrsos/interped_mrsos_crop-ctl_ecearth_signal-separated_new.nc", "r")
;ecearth mrsos unit:m^3/m^3-->*10kg/m^2i.e.1m^3/m^3=10kg/m^2
;kg m-2 mrsos    96,192
;(time,depth,lat,lon) kg m-2 mrsol
myfile4 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/crop-ctl/Eyr/tas/interped_tas_crop-ctl_ecearth_signal-separated_new.nc", "r")
tas_crop_ecearth=new((/192,288/),double) ;3varia, 3 simulation
mrsosn_crop_ecearth=new((/192,288/),double) ;
tas_crop_ecearth= dim_avg_n_Wrap(myfile4->tas_nonlocal,0) 
mrsosn_crop_ecearth = dim_avg_n_Wrap(myfile3->mrsos_nonlocal,0)
mrsos_crop_ecearth=mrsosn_crop_ecearth*10
copy_VarMeta(mrsosn_crop_ecearth, mrsos_crop_ecearth)
;.......mpiesm.......crop-ctl..............................
myfile5 = addfile("/work/bm1147/b380949/web-monitoring/crop-ctl/Lmon/mrsos/mrsos_crop-ctl_mpiesm_signal-separated.nc", "r")
;;kg m-2 mrsos;
myfile6 = addfile("/work/bm1147/b380949/web-monitoring/crop-ctl/Amon/tas/tas_crop-ctl_mpiesm_signal-separated.nc", "r")
tas_crop_mpiesm=new((/96,192/),double) ;3varia, 3 simulation
mrsos_crop_mpiesm =new((/96,192/),double) ;
tas_crop_mpiesm = dim_avg_n_Wrap(myfile6->tas_nonlocal(120:1919,:,:),0)
mrsos_crop_mpiesm = dim_avg_n_Wrap(myfile5->mrsos_nonlocal(120:1919,:,:),0)
;.......cesm.............frst-ctl...........................
myfile7 = addfile("/work/bm1147/b380949/web-monitoring/CESM/frst-ctl/Lmon/SOILWATER_10CM/SOILWATER_10CM_frst-ctl_cesm_signal-separated.nc", "r")
myfile8 = addfile("/work/bm1147/b380949/web-monitoring/CESM/frst-ctl/Lmon/TREFHT/TREFHT_frst-ctl_cesm_signal-separated.nc", "r")
tas_frst_cesm=new((/192,288/),double) ;3varia, 3 simulation
mrsos_frst_cesm=new((/192,288/),double) ;
tas_frst_cesm= dim_avg_n_Wrap(myfile8->TREFHT_nonlocal,0)
mrsos_frst_cesm = dim_avg_n_Wrap(myfile7->SOILWATER_10CM_nonlocal(120:1919,:,:),0)
;........ecearth..........frst-ctl.................
myfile9 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/frst-ctl/Eyr/mrsos/interped_mrsos_frst-ctl_ecearth_signal-separated_new.nc", "r")
;ecearth mrsos unit:m^3/m^3-->*10kg/m^2i.e.1m^3/m^3=10kg/m^2
;kg m-2 mrsos
;(time,depth,lat,lon) kg m-2 mrsol
myfile10 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/frst-ctl/Eyr/tas/interped_tas_frst-ctl_ecearth_signal-separated_new.nc", "r")
tas_frst_ecearth=new((/192,288/),double) ;3varia, 3 simulation
mrsosn_frst_ecearth=new((/192,288/),double) ;
tas_frst_ecearth= dim_avg_n_Wrap(myfile10->tas_nonlocal,0)
mrsosn_frst_ecearth = dim_avg_n_Wrap(myfile9->mrsos_nonlocal,0)
mrsos_frst_ecearth = mrsosn_frst_ecearth*10
copy_VarMeta(mrsosn_frst_ecearth, mrsos_frst_ecearth)
;.........mpiesm...........frst-ctl...............
myfile11 = addfile("/work/bm1147/b380949/web-monitoring/frst-ctl/Lmon/mrsos/mrsos_frst-ctl_mpiesm_signal-separated.nc", "r")
;;kg m-2 mrsos;
myfile12 = addfile("/work/bm1147/b380949/web-monitoring/frst-ctl/Amon/tas/tas_frst-ctl_mpiesm_signal-separated.nc", "r")
tas_frst_mpiesm=new((/96,192/),double) ;3varia, 3 simulation
mrsos_frst_mpiesm=new((/96,192/),double) ;
tas_frst_mpiesm= dim_avg_n_Wrap(myfile12->tas_nonlocal(120:1919,:,:),0)
mrsos_frst_mpiesm = dim_avg_n_Wrap(myfile11->mrsos_nonlocal(120:1919,:,:),0)
;........ecearth..........irr-crop.................
myfile13 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/irr-crop/Eyr/mrsos/interped_mrsos_irr-crop_ecearth_signal-separated_new.nc", "r")
;ecearth mrsos unit:m^3/m^3-->*10kg/m^2i.e.1m^3/m^3=10kg/m^2
;kg m-2 mrsos
;(time,depth,lat,lon) kg m-2 mrsol
myfile14 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/irr-crop/Eyr/tas/interped_tas_irr-crop_ecearth_signal-separated_new.nc", "r")
tas_irri_ecearth=new((/192,288/),double) ;3varia, 3 simulation
mrsosn_irri_ecearth=new((/192,288/),double) ;
tas_irri_ecearth= dim_avg_n_Wrap(myfile14->tas_nonlocal,0)
mrsosn_irri_ecearth = dim_avg_n_Wrap(myfile13->mrsos_nonlocal,0)
mrsos_irri_ecearth=mrsosn_irri_ecearth*10
copy_VarMeta(mrsosn_irri_ecearth, mrsos_irri_ecearth)
;........cesm............irr-crop.................
myfile15 = addfile("/work/bm1147/b380949/web-monitoring/CESM/irri-crop/Lmon/SOILWATER_10CM/SOILWATER_10CM_irri-crop_cesm_signal-separated.nc", "r")
myfile16 = addfile("/work/bm1147/b380949/web-monitoring/CESM/irri-crop/Lmon/TREFHT/TREFHT_irri-crop_cesm_signal-separated.nc", "r")
tas_irri_cesm=new((/192,288/),double) ;3varia, 3 simulation
mrsos_irri_cesm=new((/192,288/),double) ;
tas_irri_cesm= dim_avg_n_Wrap(myfile16->TREFHT_nonlocal,0)
mrsos_irri_cesm = dim_avg_n_Wrap(myfile15->SOILWATER_10CM_nonlocal(120:1919,:,:),0)
;........mpiesm.....irr-crop.......................
myfile17 = addfile("/work/bm1147/b380949/web-monitoring/irri-crop/Lmon/mrsos/mrsos_irri-crop_mpiesm_signal-separated.nc", "r")
;;kg m-2 mrsos;
myfile18 = addfile("/work/bm1147/b380949/web-monitoring/irri-crop/Amon/tas/tas_irri-crop_mpiesm_signal-separated.nc", "r")
tas_irri_mpiesm=new((/96,192/),double) ;3varia, 3 simulation
mrsos_irri_mpiesm=new((/96,192/),double) ;
tas_irri_mpiesm= dim_avg_n_Wrap(myfile18->tas_nonlocal(120:1919,:,:),0)
mrsos_irri_mpiesm = dim_avg_n_Wrap(myfile17->mrsos_nonlocal(120:1919,:,:),0)
;..............plot.................................................
tas_crop_cesm_revised=mask(tas_crop_cesm, tas_crop_cesm.eq.0, False);crop;4,192,288;4:veg,soil,litter,land
copy_VarMeta(tas_crop_cesm, tas_crop_cesm_revised)
mrsos_crop_cesm_revised=mask(mrsos_crop_cesm, mrsos_crop_cesm.eq.0, False);crop;4,192,288;4:veg,soil,litter,land
copy_VarMeta(mrsos_crop_cesm, mrsos_crop_cesm_revised)

tas_frst_cesm_revised=mask(tas_frst_cesm, tas_frst_cesm.eq.0, False);frst;4,192,288;4:veg,soil,litter,land
copy_VarMeta(tas_frst_cesm, tas_frst_cesm_revised)
mrsos_frst_cesm_revised=mask(mrsos_frst_cesm, mrsos_frst_cesm.eq.0, False);frst;4,192,288;4:veg,soil,litter,land
copy_VarMeta(mrsos_frst_cesm, mrsos_frst_cesm_revised)

tas_irri_cesm_revised=mask(tas_irri_cesm, tas_irri_cesm.eq.0, False);irri;4,192,288;4:veg,soil,litter,land
copy_VarMeta(tas_irri_cesm, tas_irri_cesm_revised)
mrsos_irri_cesm_revised=mask(mrsos_irri_cesm, mrsos_irri_cesm.eq.0, False);irri;4,192,288;4:veg,soil,litter,land
copy_VarMeta(mrsos_irri_cesm, mrsos_irri_cesm_revised)

tas_crop_mpiesm_revised=mask(tas_crop_mpiesm, tas_crop_mpiesm.eq.0, False);crop;4,192,288;4:veg,soil,litter,land
copy_VarMeta(tas_crop_mpiesm, tas_crop_mpiesm_revised)
mrsos_crop_mpiesm_revised=mask(mrsos_crop_mpiesm, mrsos_crop_mpiesm.eq.0, False);crop;4,192,288;4:veg,soil,litter,land
copy_VarMeta(mrsos_crop_mpiesm, mrsos_crop_mpiesm_revised)

tas_frst_mpiesm_revised=mask(tas_frst_mpiesm, tas_frst_mpiesm.eq.0, False);frst;4,192,288;4:veg,soil,litter,land
copy_VarMeta(tas_frst_mpiesm, tas_frst_mpiesm_revised)
mrsos_frst_mpiesm_revised=mask(mrsos_frst_mpiesm, mrsos_frst_mpiesm.eq.0, False);frst;4,192,288;4:veg,soil,litter,land
copy_VarMeta(mrsos_frst_mpiesm, mrsos_frst_mpiesm_revised)

tas_irri_mpiesm_revised=mask(tas_irri_mpiesm, tas_irri_mpiesm.eq.0, False);irri;4,192,288;4:veg,soil,litter,land
copy_VarMeta(tas_irri_mpiesm, tas_irri_mpiesm_revised)
mrsos_irri_mpiesm_revised=mask(mrsos_irri_mpiesm, mrsos_irri_mpiesm.eq.0, False);irri;4,192,288;4:veg,soil,litter,land
copy_VarMeta(mrsos_irri_mpiesm, mrsos_irri_mpiesm_revised)

tas_crop_ecearth_revised=mask(tas_crop_ecearth, tas_crop_ecearth.eq.0, False);crop;4,192,288;4:veg,soil,litter,land
copy_VarMeta(tas_crop_ecearth, tas_crop_ecearth_revised)
mrsos_crop_ecearth_revised=mask(mrsos_crop_ecearth, mrsos_crop_ecearth.eq.0, False);crop;4,192,288;4:veg,soil,litter,land
copy_VarMeta(mrsos_crop_ecearth, mrsos_crop_ecearth_revised)

tas_frst_ecearth_revised=mask(tas_frst_ecearth, tas_frst_ecearth.eq.0, False);frst;4,192,288;4:veg,soil,litter,land
copy_VarMeta(tas_frst_ecearth, tas_frst_ecearth_revised)
mrsos_frst_ecearth_revised=mask(mrsos_frst_ecearth, mrsos_frst_ecearth.eq.0, False);frst;4,192,288;4:veg,soil,litter,land
copy_VarMeta(mrsos_frst_ecearth, mrsos_frst_ecearth_revised)

tas_irri_ecearth_revised=mask(tas_irri_ecearth, tas_irri_ecearth.eq.0, False);irri;4,192,288;4:veg,soil,litter,land
copy_VarMeta(tas_irri_ecearth, tas_irri_ecearth_revised)
mrsos_irri_ecearth_revised=mask(mrsos_irri_ecearth, mrsos_irri_ecearth.eq.0, False);irri;4,192,288;4:veg,soil,litter,land
copy_VarMeta(mrsos_irri_ecearth, mrsos_irri_ecearth_revised)
;..........................................................................
wks_type = "pdf"
wks_type@wkPaperWidthF =  9.0
wks_type@wkPaperHeightF = 9.0
wks = gsn_open_wks(wks_type,"climate_map_revised")
gsn_define_colormap(wks,"temp_19lev")
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
res@cnLevels  =  (/-2,-1.5,-1.0,-0.5,0,0.5,1.0,1.5,2.0/)
res@cnFillColors = (/2,4,6,8,10,12,14,16,18,20/)
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
lres@gsnLeftStringParallelPosF= 0.0      ;0.04
lres@gsnLeftStringOrthogonalPosF= 0.0 ;-0.14
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
printVarSummary(tas_crop_mpiesm({-60:90},:))
printVarSummary(tas_crop_cesm({-60:90},:))
printVarSummary(tas_crop_ecearth({-60:90},:))
y = tas_crop_mpiesm&lat   ;y(16::)  -60to90 80points
y_cesm = tas_crop_cesm&lat
y_ecearth = tas_crop_ecearth&lat

crop_zonal_tas=dim_avg_n_Wrap(tas_crop_mpiesm({-60:90},:), 1)
crop_zonal_cesm_tas=dim_avg_n_Wrap(tas_crop_cesm({-60:90},:), 1)
crop_zonal_ecearth_tas=dim_avg_n_Wrap(tas_crop_ecearth({-60:90},:), 1)
printVarSummary(crop_zonal_ecearth_tas)
;print(crop_zonal_ecearth_tas)
frst_zonal_tas=dim_avg_n_Wrap(tas_frst_mpiesm({-60:90},:), 1)
frst_zonal_cesm_tas=dim_avg_n_Wrap(tas_frst_cesm({-60:90},:), 1)
frst_zonal_ecearth_tas=dim_avg_n_Wrap(tas_frst_ecearth({-60:90},:), 1)

irri_zonal_tas=dim_avg_n_Wrap(tas_irri_mpiesm({-60:90},:), 1)
irri_zonal_cesm_tas=dim_avg_n_Wrap(tas_irri_cesm({-60:90},:), 1)
irri_zonal_ecearth_tas=dim_avg_n_Wrap(tas_irri_ecearth({-60:90},:), 1)

crop_zonal_mrsos=dim_avg_n_Wrap(mrsos_crop_mpiesm({-60:90},:), 1)
crop_zonal_cesm_mrsos=dim_avg_n_Wrap(mrsos_crop_cesm({-60:90},:), 1)
crop_zonal_ecearth_mrsos=dim_avg_n_Wrap(mrsos_crop_ecearth({-60:90},:), 1)

frst_zonal_mrsos=dim_avg_n_Wrap(mrsos_frst_mpiesm({-60:90},:), 1)
frst_zonal_cesm_mrsos=dim_avg_n_Wrap(mrsos_frst_cesm({-60:90},:), 1)
frst_zonal_ecearth_mrsos=dim_avg_n_Wrap(mrsos_frst_ecearth({-60:90},:), 1)

irri_zonal_mrsos=dim_avg_n_Wrap(mrsos_irri_mpiesm({-60:90},:), 1)
irri_zonal_cesm_mrsos=dim_avg_n_Wrap(mrsos_irri_cesm({-60:90},:), 1)
irri_zonal_ecearth_mrsos=dim_avg_n_Wrap(mrsos_irri_ecearth({-60:90},:), 1)
;...plot.............................
plot1=new(12,graphic)
plot3=new(12,graphic)
plot2=new(12,graphic)
res@tiXAxisFont="times-roman"
res@tiXAxisFontHeightF=0.01
res@tiYAxisFont="times-roman"
res@tiMainFont="times-roman"
res@tiMainFontHeightF=0.019  ; 0.018
res@tiMainOffsetYF =  -0.017        ;-0.008

res@tiYAxisFontHeightF=0.019
res@gsnLeftStringParallelPosF= 0.0 ;0.04
res@gsnLeftStringOrthogonalPosF= 0.0    ;  -0.11
 res@gsnRightString    = " "
res@tmXBLabelFontHeightF =0.014
res@tmYLLabelFontHeightF =0.014
res@gsnStringFont ="times-Bold"
res@gsnLeftStringFontHeightF=0.019  ;0.018
res@gsnRightStringFontHeightF=0.011
res@gsnCenterStringFontHeightF=0.018
res@tmXBLabelFont         = "times-roman"
res@tmYLLabelFont         = "times-roman"    
  lres@xyLineColor= (/0.75,0.20,0.1,1./) ;red 
   res@cnMonoFillColor=False
    res@cnMissingValFillColor="white"
     res@tiYAxisString  = "CROP-CTL"
     res@tiMainString ="MPI-ESM"
    res@gsnLeftString    = "(a)"
    plot3(0)=gsn_csm_contour_map(wks, tas_crop_mpiesm_revised(:,:), res)
    res@tiYAxisString  = " "
    res@tiMainString ="CESM"
    res@gsnLeftString    = "(b)"
    plot3(1)=gsn_csm_contour_map(wks, tas_crop_cesm_revised(:,:), res)
    res@tiMainString ="EC-Earth"
    res@gsnLeftString    = "(c)"
    plot3(2)=gsn_csm_contour_map(wks, tas_crop_ecearth_revised(:,:), res)
    res@gsnCenterString = ""
    res@tiMainString =""

    lres@gsnStringFont ="times-Bold"
    lres@gsnLeftStringFontHeightF=0.019
    lres@gsnLeftString    = "(d)"
    lres@trXMinF= -1.6
    lres@trXMaxF=1.1
    lres@tmXBMode              = "Automatic"
    lres@tmXBMode              = "Explicit"
    lres@tmXBValues            = (/-1.5,-1.0,-0.5,0,0.5,1.0,1.0,1.0/)
    lres@tmXBLabels  = (/"","-1","-0.5",\  
                     "0","0.5","1.0","",""/)
    plot3(3)=gsn_csm_xy(wks,crop_zonal_ecearth_tas(:),y_ecearth(32::),lres) ;ecearth
   res@tiYAxisString  = "FRST-CTL"
    res@gsnLeftString    = "(e)"
    plot3(4)=gsn_csm_contour_map(wks, tas_frst_mpiesm_revised(:,:), res)
        res@tiYAxisString  = " "
    res@gsnLeftString    = "(f)"
    plot3(5)=gsn_csm_contour_map(wks, tas_frst_cesm_revised(:,:), res)
    res@gsnLeftString    = "(g)"
    plot3(6)=gsn_csm_contour_map(wks, tas_frst_ecearth_revised(:,:), res)

    lres@gsnLeftString    = "(h)"
    plot3(7)=gsn_csm_xy(wks,frst_zonal_ecearth_tas(:),y_ecearth(32::),lres) ;ecearth
   res@tiYAxisString  = "IRR-CROP"
    res@gsnLeftString    = "(i)"
    plot3(8)=gsn_csm_contour_map(wks, tas_irri_mpiesm_revised(:,:), res)
        res@tiYAxisString  = " "
    res@gsnLeftString    = "(j)"
    plot3(9)=gsn_csm_contour_map(wks, tas_irri_cesm_revised(:,:), res)
    res@gsnLeftString    = "(k)"
    plot3(10)=gsn_csm_contour_map(wks, tas_irri_ecearth_revised(:,:), res)
    lres@gsnLeftString    = "(l)"
    lres@tiXAxisString= "";"kg/m~S~2~N~"
    lres@tiXAxisFont="times-roman"
    lres@tiXAxisFontHeightF=0.018
    plot3(11)=gsn_csm_xy(wks,irri_zonal_ecearth_tas(:),y_ecearth(32::),lres) ;ecearth
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

    plot1(0)=gsn_csm_contour(wks, tas_crop_mpiesm(:,:), res)
    plot1(1)=gsn_csm_contour(wks, tas_crop_cesm(:,:), res)
    plot1(2)=gsn_csm_contour(wks, tas_crop_ecearth(:,:), res)

    plot2(0)=gsn_csm_contour(wks, tas_crop_mpiesm(:,:), res)
    plot2(1)=gsn_csm_contour(wks, tas_crop_cesm(:,:), res)
    plot2(2)=gsn_csm_contour(wks, tas_crop_ecearth(:,:), res)

    plot1(3)=gsn_csm_xy(wks,crop_zonal_cesm_tas(:),y_cesm(32::),lres)
      lres@xyLineColor= (/0.11,0.44,0.93,1./) ;blue
    plot2(3)=gsn_csm_xy(wks,crop_zonal_tas(:),y(16::),lres) ;MPI-ESM
    lres@xyLineColor= (/0.0,0.0,0.0,1./)
    plot1(4)=gsn_csm_contour(wks, tas_frst_mpiesm(:,:), res)
    plot1(5)=gsn_csm_contour(wks, tas_frst_cesm(:,:), res)
    plot1(6)=gsn_csm_contour(wks, tas_frst_ecearth(:,:), res)

    plot2(4)=gsn_csm_contour(wks, tas_frst_mpiesm(:,:), res)
    plot2(5)=gsn_csm_contour(wks, tas_frst_cesm(:,:), res)
    plot2(6)=gsn_csm_contour(wks, tas_frst_ecearth(:,:), res)
    plot1(7)=gsn_csm_xy(wks,frst_zonal_cesm_tas(:),y_cesm(32::),lres)
    lres@xyLineColor= (/0.11,0.44,0.93,1./) ;blue


    plot2(7)=gsn_csm_xy(wks,frst_zonal_tas(:),y(16::),lres) ;MPI-ESM
   
    plot1(8)=gsn_csm_contour(wks, tas_irri_mpiesm(:,:), res)
    plot1(9)=gsn_csm_contour(wks, tas_irri_cesm(:,:), res)
    plot1(10)=gsn_csm_contour(wks, tas_irri_ecearth(:,:), res)

    plot2(8)=gsn_csm_contour(wks, tas_irri_mpiesm(:,:), res)
    plot2(9)=gsn_csm_contour(wks, tas_irri_cesm(:,:), res)
    plot2(10)=gsn_csm_contour(wks, tas_irri_ecearth(:,:), res)
     lres@xyLineColor= (/0.0,0.0,0.0,1./)
    plot1(11)=gsn_csm_xy(wks,irri_zonal_cesm_tas(:),y_cesm(32::),lres)
                        lres@xyLineColor= (/0.11,0.44,0.93,1./) ;blue
    plot2(11)=gsn_csm_xy(wks,irri_zonal_tas(:),y(16::),lres) ;MPI-ESM
       overlay(plot3(3),plot1(3))
       overlay(plot3(3),plot2(3))

       overlay(plot3(7),plot1(7))
       overlay(plot3(7),plot2(7))

       overlay(plot3(11),plot1(11))
       overlay(plot3(11),plot2(11))

lgres                    = True
lgres@lgLineThicknessF   = lres@xyLineThicknessF  ; legend line thickness = contour line thickness
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
lgres@lgLineColors       = (/(/0.75,0.20,0.1,1./),(/0.0,0.0,0.0,1./),(/0.11,0.44,0.93,1./)/)
gsn_legend_ndc(wks,3,(/"  EC-Earth","  CESM","  MPI-ESM"/),0.73,0.33,lgres)
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
        plot7(i)=gsn_add_polyline(wks, plot3(0), plx(i,:), ply(i,:), plres)
end do

plest                   = True                      ; polyline mods desired
  plest@gsLineColor       = "grey"                    ; color of lines
  plest@gsLineThicknessF  = 1.2                       ; thickness of lines
  plest@gsLineOpacityF = 0.8
  dum=new(3,graphic)
  dum(0)=gsn_add_polyline(wks,plot3(3),(/0,0/),(/-60,90/),plest)
  dum(1)=gsn_add_polyline(wks,plot3(7),(/0,0/),(/-60,90/),plest)
  dum(2)=gsn_add_polyline(wks,plot3(11),(/0,0/),(/-60,90/),plest)
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
pres@lbTitleString       ="Near surface temperature (K)  "
pres@lbTitleOffsetF      = 0.2 ;-0.2
pres@lbTitleJust         ="CenterCenter"
pres@lbTitlePosition     = "Bottom"
pres@gsnPanelYWhiteSpacePercent = 0.005
pres@gsnPanelXWhiteSpacePercent =0.0
pres@tiYAxisString=""
pres@txFontHeightF  =0.018
pres@txFont ="times-roman"
pres@gsnPanelXF = (/0.0185,0.255,0.492,0.729,0.0185,0.255,0.492,0.729,0.0185,0.255,0.492,0.729/) ;0.4
pres@gsnPanelYF = (/0.664,0.664,0.664,0.664,0.549,0.549,0.549,0.549,0.434,0.434,0.434,0.434/) ;0.4
pres@gsnPanelDebug = True
gsn_panel(wks,plot3(:),(/3,4/),pres)
end

