load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin
;...........................................
mpicon1 = addfile("/work/bm1147/b380949/web-monitoring/crop-ctl/Eyr/cVeg/cVeg_ctl.nc", "r")
mpicon2 = addfile("/work/bm1147/b380949/web-monitoring/crop-ctl/Eyr/cSoil/cSoil_ctl.nc", "r")
mpicon3 = addfile("/work/bm1147/b380949/web-monitoring/crop-ctl/Eyr/cProduct/cProduct_ctl.nc", "r")
mpicon4= addfile("/work/bm1147/b380949/web-monitoring/crop-ctl/Eyr/cLitter/cLitter_ctl.nc", "r")

myfile1 = addfile("/work/bm1147/b380949/web-monitoring/crop-ctl/Eyr/cVeg/cVeg_crop.nc", "r")
myfile2 = addfile("/work/bm1147/b380949/web-monitoring/crop-ctl/Eyr/cSoil/cSoil_crop.nc", "r")
myfile3 = addfile("/work/bm1147/b380949/web-monitoring/crop-ctl/Eyr/cProduct/cProduct_crop.nc", "r")
myfile4= addfile("/work/bm1147/b380949/web-monitoring/crop-ctl/Eyr/cLitter/cLitter_crop.nc", "r")

myfile5 = addfile("/work/bm1147/b380949/web-monitoring/frst-ctl/Eyr/cVeg/cVeg_frst.nc", "r")
myfile6 = addfile("/work/bm1147/b380949/web-monitoring/frst-ctl/Eyr/cSoil/cSoil_frst.nc", "r")
myfile7 = addfile("/work/bm1147/b380949/web-monitoring/frst-ctl/Eyr/cProduct/cProduct_frst.nc", "r")
myfile8 = addfile("/work/bm1147/b380949/web-monitoring/frst-ctl/Eyr/cLitter/cLitter_frst.nc", "r")

myfile13 = addfile("/work/bm1147/b380949/web-monitoring/irri-crop/Eyr/cVeg/cVeg_irri.nc", "r")
myfile14 = addfile("/work/bm1147/b380949/web-monitoring/irri-crop/Eyr/cSoil/cSoil_irri.nc", "r")
myfile15 = addfile("/work/bm1147/b380949/web-monitoring/irri-crop/Eyr/cProduct/cProduct_irri.nc", "r")
myfile16 = addfile("/work/bm1147/b380949/web-monitoring/irri-crop/Eyr/cLitter/cLitter_irri.nc", "r")

chessboard_mpi=addfile("/home/b/b380949/chessboard/lamaclima_experiments_chessboard_pattern_mpiesm.nc", "r")
myfile = addfile("/home/b/b380949/gridarea/gridarea_sea0.nc","r")
my=addfile("/home/b/b380949/gridarea/AREA_CESM.nc","r")

cb_mpi=chessboard_mpi->chessboard_pattern(::-1,:)

gridarea = myfile->slm
gridarea1= my->area
gridarea2=todouble(gridarea1)
copy_VarMeta(gridarea1, gridarea2)
gridarea0=where(ismissing(gridarea2), 0, gridarea2)
copy_VarMeta(gridarea2, gridarea0)
gridarea2n=gridarea0*1000000 ;km^2->m^2
copy_VarMeta(gridarea2, gridarea2n)
;...........Read MPIESM data and calculate the average for the last 30 years
;..........ctl......................
mpicon=new((/5,96,192/),double) ;ctl 2015-2176
mpicon(0,:,:) = dim_avg_n_Wrap(mpicon1->cVeg(130:159,:,:),0)
mpicon(1,:,:) =dim_avg_n_Wrap(mpicon2->cSoil(130:159,:,:),0)
mpicon(2,:,:) =dim_avg_n_Wrap(mpicon3->cProduct(130:159,:,:),0)
mpicon(3,:,:) =dim_avg_n_Wrap(mpicon4->cLitter(130:159,:,:),0)
mpicon(4,:,:)=dim_sum_n_Wrap(mpicon(0:3,:,:),0)
;..............crop expansion..........
pool=new((/5,96,192/),double) ;crop 2015-2176
pool(0,:,:) = dim_avg_n_Wrap(myfile1->cVeg(130:159,:,:),0)
pool(1,:,:) =dim_avg_n_Wrap(myfile2->cSoil(130:159,:,:),0)
pool(2,:,:) =dim_avg_n_Wrap(myfile3->cProduct(130:159,:,:),0)
pool(3,:,:) =dim_avg_n_Wrap(myfile4->cLitter(130:159,:,:),0)
pool(4,:,:)=dim_sum_n_Wrap(pool(0:3,:,:),0)
;..............afforestation...................
pool2=new((/5,96,192/),double) ;aff:2015-2175
pool2(0,:,:) = dim_avg_n_Wrap(myfile5->cVeg(130:159,:,:),0)
pool2(1,:,:) =dim_avg_n_Wrap(myfile6->cSoil(130:159,:,:),0)
pool2(2,:,:) =dim_avg_n_Wrap(myfile7->cProduct(130:159,:,:),0)
pool2(3,:,:) =dim_avg_n_Wrap(myfile8->cLitter(130:159,:,:),0)
pool2(4,:,:)=dim_sum_n_Wrap(pool2(0:3,:,:),0)
;.....................irrigation......................................
pool4=new((/5,96,192/),double) ;irri:2015-2174
pool4(0,:,:) = dim_avg_n_Wrap(myfile13->cVeg(130:159,:,:),0)
pool4(1,:,:) =dim_avg_n_Wrap(myfile14->cSoil(130:159,:,:),0)
pool4(2,:,:) =dim_avg_n_Wrap(myfile15->cProduct(130:159,:,:),0)
pool4(3,:,:) =dim_avg_n_Wrap(myfile16->cLitter(130:159,:,:),0)
pool4(4,:,:)=dim_sum_n_Wrap(pool4(0:3,:,:),0)
;.......calculat the difference at first anf then mask with chessboard ........ cVeg
diffveg=new((/3,96,192/),double)
diffveg(0,:,:)=pool(1,:,:)-mpicon(1,:,:);crop-ctl
diffveg(1,:,:)=pool2(1,:,:)-mpicon(1,:,:);frst-ctl
diffveg(2,:,:)=pool4(1,:,:)-pool(1,:,:);irri-crop
copy_VarMeta(pool(1,:,:), diffveg(0,:,:))
;.......................................
diffveg_mask=mask(diffveg,cb_mpi, 0)
copy_VarMeta(pool(1,:,:), diffveg_mask(0,:,:))
;.......global vs regional..............
;.......global...........
total_g=wgt_areasum2(diffveg(:,{-90:90},{0:360}),gridarea({-90:90},{0:360}),0)*0.000000000001  ;kg/m^2->GtC
nonlocal_g=wgt_areasum2(diffveg_mask(:,{-90:90},{0:360}),gridarea({-90:90},{0:360}),0)*0.000000000001  ;kg/m^2->GtC since gridarerevised unit2 km

total_g!0="scenario"
nonlocal_g!0="scenario"
;.region:North_America{28:54},{260:300} ;Amazon{-16:3},{283:308} ;Congo{-9:9},{13:29} 
;North Eurasia{56:67},{32:92};East Asia{46:62},{106:140};Southeast Asia {22:36},{101:119}
;Australia {-23:-12},{115:150}
total_NA=wgt_areasum2(diffveg(:,{28:54},{260:300}),gridarea({28:54},{260:300}),0)*0.000000000001
total_AM=wgt_areasum2(diffveg(:,{-16:3},{283:308}),gridarea({-16:3},{283:308}),0)*0.000000000001
total_CG=wgt_areasum2(diffveg(:,{-9:9},{13:29}),gridarea({-9:9},{13:29}),0)*0.000000000001
total_NE=wgt_areasum2(diffveg(:,{56:67},{32:92}),gridarea({56:67},{32:92}),0)*0.000000000001
total_EA=wgt_areasum2(diffveg(:,{46:62},{106:140}),gridarea({46:62},{106:140}),0)*0.000000000001
total_SEA=wgt_areasum2(diffveg(:,{22:36},{101:119}),gridarea({22:36},{101:119}),0)*0.000000000001
total_AU=wgt_areasum2(diffveg(:,{-23:-12},{115:150}),gridarea({-23:-12},{115:150}),0)*0.000000000001

nonlocal_NA=wgt_areasum2(diffveg_mask(:,{28:54},{260:300}),gridarea({28:54},{260:300}),0)*0.000000000001
nonlocal_AM=wgt_areasum2(diffveg_mask(:,{-16:3},{283:308}),gridarea({-16:3},{283:308}),0)*0.000000000001
nonlocal_CG=wgt_areasum2(diffveg_mask(:,{-9:9},{13:29}),gridarea({-9:9},{13:29}),0)*0.000000000001
nonlocal_NE=wgt_areasum2(diffveg_mask(:,{56:67},{32:92}),gridarea({56:67},{32:92}),0)*0.000000000001
nonlocal_EA=wgt_areasum2(diffveg_mask(:,{46:62},{106:140}),gridarea({46:62},{106:140}),0)*0.000000000001
nonlocal_SEA=wgt_areasum2(diffveg_mask(:,{22:36},{101:119}),gridarea({22:36},{101:119}),0)*0.000000000001
nonlocal_AU=wgt_areasum2(diffveg_mask(:,{-23:-12},{115:150}),gridarea({-23:-12},{115:150}),0)*0.000000000001

ratio_mpi=new((/8,3/),double)
ratio_mpi(0,:)=nonlocal_g/total_g
ratio_mpi(1,:)=nonlocal_NA/total_NA
ratio_mpi(2,:)=nonlocal_AM/total_AM
ratio_mpi(3,:)=nonlocal_CG/total_CG
ratio_mpi(4,:)=nonlocal_NE/total_NE
ratio_mpi(5,:)=nonlocal_EA/total_EA
ratio_mpi(6,:)=nonlocal_SEA/total_SEA
ratio_mpi(7,:)=nonlocal_AU/total_AU

;calculat the difference at first anf then mask with chessboard 
;after that calculate temperal and spacial mean.
;then divided.
;....Read CESM data and calculate the average for the last 30 years..........
control12 = addfile("/work/bm1147/b380949/web-monitoring/CESM/crop-ctl/Lmon/TOTVEGC/TOTVEGC_ctl.nc", "r")
control22 = addfile("/work/bm1147/b380949/web-monitoring/CESM/crop-ctl/Lmon/TOTSOMC/TOTSOMC_ctl.nc", "r")
control42 = addfile("/work/bm1147/b380949/web-monitoring/CESM/crop-ctl/Lmon/TOTLITC/TOTLITC_ctl.nc", "r")
control32 = addfile("/work/bm1147/b380949/web-monitoring/CESM/crop-ctl/Lmon/TOTECOSYSC/TOTECOSYSC_ctl.nc", "r")

myfile12 = addfile("/work/bm1147/b380949/web-monitoring/CESM/crop-ctl/Lmon/TOTVEGC/TOTVEGC_crop.nc", "r")
myfile22 = addfile("/work/bm1147/b380949/web-monitoring/CESM/crop-ctl/Lmon/TOTSOMC/TOTSOMC_crop.nc", "r")
myfile42 = addfile("/work/bm1147/b380949/web-monitoring/CESM/crop-ctl/Lmon/TOTLITC/TOTLITC_crop.nc", "r")
myfile32 = addfile("/work/bm1147/b380949/web-monitoring/CESM/crop-ctl/Lmon/TOTECOSYSC/TOTECOSYSC_crop.nc", "r")

myfile52 = addfile("/work/bm1147/b380949/web-monitoring/CESM/frst-ctl/Lmon/TOTVEGC/TOTVEGC_frst.nc", "r")
myfile62 = addfile("/work/bm1147/b380949/web-monitoring/CESM/frst-ctl/Lmon/TOTSOMC/TOTSOMC_frst.nc", "r")
myfile82 = addfile("/work/bm1147/b380949/web-monitoring/CESM/frst-ctl/Lmon/TOTLITC/TOTLITC_frst.nc", "r")
myfile72 = addfile("/work/bm1147/b380949/web-monitoring/CESM/frst-ctl/Lmon/TOTECOSYSC/TOTECOSYSC_frst.nc", "r")

myfile132 = addfile("/work/bm1147/b380949/web-monitoring/CESM/irri-crop/Lmon/TOTVEGC/TOTVEGC_irri.nc", "r")
myfile142 = addfile("/work/bm1147/b380949/web-monitoring/CESM/irri-crop/Lmon/TOTSOMC/TOTSOMC_irri.nc", "r")
myfile152 = addfile("/work/bm1147/b380949/web-monitoring/CESM/irri-crop/Lmon/TOTLITC/TOTLITC_irri.nc", "r")
myfile162 = addfile("/work/bm1147/b380949/web-monitoring/CESM/irri-crop/Lmon/TOTECOSYSC/TOTECOSYSC_irri.nc", "r")
myfile=addfile("/home/b/b380949/gridarea/landmask_glacier_cesm.nc", "r")
slm = myfile->landmask

chessboard_cesm=addfile("/home/b/b380949/chessboard/lamaclima_experiments_chessboard_pattern_cesm.nc", "r")
cb_cesm=chessboard_cesm->chessboard_pattern
;.....ctl................
cesmcon=new((/4,192,288/),double) ;crop
cesmcon(0,:,:) = dim_avg_n_Wrap(control12->TOTVEGC(1560:1919,:,:),0)
cesmcon(1,:,:) =dim_avg_n_Wrap(control22->TOTSOMC(1560:1919,:,:),0)
cesmcon(2,:,:) =dim_avg_n_Wrap(control42->TOTLITC(1560:1919,:,:),0)
cesmcon(3,:,:) =dim_avg_n_Wrap(control32->TOTECOSYSC(1560:1919,:,:),0)
;......crop.........
pool2n=new((/4,192,288/),double) ;crop
pool2n(0,:,:) = dim_avg_n_Wrap(myfile12->TOTVEGC(1560:1919,:,:),0)
pool2n(1,:,:) =dim_avg_n_Wrap(myfile22->TOTSOMC(1560:1919,:,:),0)
pool2n(2,:,:) =dim_avg_n_Wrap(myfile42->TOTLITC(1560:1919,:,:),0)
pool2n(3,:,:) =dim_avg_n_Wrap(myfile32->TOTECOSYSC(1560:1919,:,:),0)
;.............afforestation..........
flux2n=new((/4,192,288/),double)  ;frst
flux2n(0,:,:) = dim_avg_n_Wrap(myfile52->TOTVEGC(1560:1919,:,:),0)
flux2n(1,:,:) =dim_avg_n_Wrap(myfile62->TOTSOMC(1560:1919,:,:),0)
flux2n(2,:,:) =dim_avg_n_Wrap(myfile82->TOTLITC(1560:1919,:,:),0)
flux2n(3,:,:) =dim_avg_n_Wrap(myfile72->TOTECOSYSC(1560:1919,:,:),0)
;................irri...................
pool23n=new((/4,192,288/),double)  ;irrigation
pool23n(0,:,:)= dim_avg_n_Wrap(myfile132->TOTVEGC(1560:1919,:,:),0)
pool23n(1,:,:)=dim_avg_n_Wrap(myfile142->TOTSOMC(1560:1919,:,:),0)
pool23n(2,:,:)= dim_avg_n_Wrap(myfile152->TOTLITC(1560:1919,:,:),0)
pool23n(3,:,:)= dim_avg_n_Wrap(myfile162->TOTECOSYSC(1560:1919,:,:),0)
;.......g-->kg
cesmconm=cesmcon*0.001
pool2m=pool2n*0.001
flux2m=flux2n*0.001
pool23m=pool23n*0.001
copy_VarMeta(cesmcon, cesmconm)
copy_VarMeta(pool2n, pool2m)
copy_VarMeta(flux2n, flux2m)
copy_VarMeta(pool23n, pool23m)
;.......deal with missing value...........
cesmcon2=mask(cesmconm,  slm.eq.0, False)
copy_VarMeta(cesmcon, cesmcon2)
pool22=mask(pool2m,  slm.eq.0, False);crop
copy_VarMeta(pool2n, pool22)
pool23=mask(pool23m, slm.eq.0, False);irri
copy_VarMeta(pool23n, pool23)
flux2=mask(flux2m, slm.eq.0, False);frst
copy_VarMeta(flux2n, flux2)
;.......calculat the difference at first anf then mask with chessboard ........ cVeg
diffveg_cesm=new((/3,192,288/),double)
diffveg_cesm(0,:,:)=pool22(1,:,:)-cesmcon2(1,:,:);crop-ctl
diffveg_cesm(1,:,:)=flux2(1,:,:)-cesmcon2(1,:,:);frst-ctl
diffveg_cesm(2,:,:)=pool23(1,:,:)-pool22(1,:,:);irri-crop
copy_VarMeta(pool22(1,:,:), diffveg_cesm(0,:,:))
;.......................................
diffveg_cesm_mask=mask(diffveg_cesm,cb_cesm, 0)
copy_VarMeta(pool22(1,:,:), diffveg_cesm_mask(0,:,:))
;.......global vs regional..............
;.......global...........
total_cesm_g=wgt_areasum2(diffveg_cesm,gridarea2n,0)*0.000000000001  ;kg/m^2->GtC
nonlocal_cesm_g=wgt_areasum2(diffveg_cesm_mask,gridarea2n,0)*0.000000000001  ;kg/m^2->GtC since gridarerevised unit2 km
total_cesm_g!0="scenario"
nonlocal_cesm_g!0="scenario"
;.region:North_America{28:54},{260:300} ;Amazon{-16:3},{283:308} ;Congo{-9:9},{13:29} 
;North Eurasia{56:67},{32:92};East Asia{46:62},{106:140};Southeast Asia {22:36},{101:119}
;Australia {-23:-12},{115:150}
total_cesm_NA=wgt_areasum2(diffveg_cesm(:,{28:54},{260:300}),gridarea2n({28:54},{260:300}),0)*0.000000000001
total_cesm_AM=wgt_areasum2(diffveg_cesm(:,{-16:3},{283:308}),gridarea2n({-16:3},{283:308}),0)*0.000000000001
total_cesm_CG=wgt_areasum2(diffveg_cesm(:,{-9:9},{13:29}),gridarea2n({-9:9},{13:29}),0)*0.000000000001
total_cesm_NE=wgt_areasum2(diffveg_cesm(:,{56:67},{32:92}),gridarea2n({56:67},{32:92}),0)*0.000000000001
total_cesm_EA=wgt_areasum2(diffveg_cesm(:,{46:62},{106:140}),gridarea2n({46:62},{106:140}),0)*0.000000000001
total_cesm_SEA=wgt_areasum2(diffveg_cesm(:,{22:36},{101:119}),gridarea2n({22:36},{101:119}),0)*0.000000000001
total_cesm_AU=wgt_areasum2(diffveg_cesm(:,{-23:-12},{115:150}),gridarea2n({-23:-12},{115:150}),0)*0.000000000001

nonlocal_cesm_NA=wgt_areasum2(diffveg_cesm_mask(:,{28:54},{260:300}),gridarea2n({28:54},{260:300}),0)*0.000000000001
nonlocal_cesm_AM=wgt_areasum2(diffveg_cesm_mask(:,{-16:3},{283:308}),gridarea2n({-16:3},{283:308}),0)*0.000000000001
nonlocal_cesm_CG=wgt_areasum2(diffveg_cesm_mask(:,{-9:9},{13:29}),gridarea2n({-9:9},{13:29}),0)*0.000000000001
nonlocal_cesm_NE=wgt_areasum2(diffveg_cesm_mask(:,{56:67},{32:92}),gridarea2n({56:67},{32:92}),0)*0.000000000001
nonlocal_cesm_EA=wgt_areasum2(diffveg_cesm_mask(:,{46:62},{106:140}),gridarea2n({46:62},{106:140}),0)*0.000000000001
nonlocal_cesm_SEA=wgt_areasum2(diffveg_cesm_mask(:,{22:36},{101:119}),gridarea2n({22:36},{101:119}),0)*0.000000000001
nonlocal_cesm_AU=wgt_areasum2(diffveg_cesm_mask(:,{-23:-12},{115:150}),gridarea2n({-23:-12},{115:150}),0)*0.000000000001

ratio_cesm=new((/8,3/),double);8 region;3:crop-ctl;frst-ctl;irri-crop
ratio_cesm(0,:)=nonlocal_cesm_g/total_cesm_g
ratio_cesm(1,:)=nonlocal_cesm_NA/total_cesm_NA
ratio_cesm(2,:)=nonlocal_cesm_AM/total_cesm_AM
ratio_cesm(3,:)=nonlocal_cesm_CG/total_cesm_CG
ratio_cesm(4,:)=nonlocal_cesm_NE/total_cesm_NE
ratio_cesm(5,:)=nonlocal_cesm_EA/total_cesm_EA
ratio_cesm(6,:)=nonlocal_cesm_SEA/total_cesm_SEA
ratio_cesm(7,:)=nonlocal_cesm_AU/total_cesm_AU
;......Read EC-Earth data and calculate the average for the last 30 years...............
landsea_ecearth=addfile("/work/bm1147/b380949/EC-earth/dataset/anlysisneed/landmask_ecearth_regridded.nc","r")
landmask_ecearth=landsea_ecearth->landmask(0,:,:)
myfile13n = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/crop-ctl/Eyr/cVeg/cVeg_crop-ctl_climate_beforeinter.nc", "r")
myfile23 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/crop-ctl/Eyr/cSoil/cSoil_crop-ctl_climate_beforeinter.nc", "r")
myfile43 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/crop-ctl/Eyr/cLitter/cLitter_crop-ctl_climate_beforeinter.nc", "r")
myfile33 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/crop-ctl/Eyr/cLand/cLand_crop-ctl_climate_beforeinter.nc", "r")

myfile53 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/frst-ctl/Eyr/cVeg/cVeg_frst-ctl_climate_beforeinter.nc", "r")
myfile63 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/frst-ctl/Eyr/cSoil/cSoil_frst-ctl_climate_beforeinter.nc", "r")
myfile83 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/frst-ctl/Eyr/cLitter/cLitter_frst-ctl_climate_beforeinter.nc", "r")
myfile73 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/frst-ctl/Eyr/cLand/cLand_frst-ctl_climate_beforeinter.nc", "r")

myfile133 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/irr-crop/Eyr/cVeg/cVeg_irr-crop_climate_beforeinter.nc", "r")
myfile143 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/irr-crop/Eyr/cSoil/cSoil_irr-crop_climate_beforeinter.nc", "r")
myfile153 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/irr-crop/Eyr/cLitter/cLitter_irr-crop_climate_beforeinter.nc", "r")
myfile163 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/irr-crop/Eyr/cLand/cLand_irr-crop_climate_beforeinter.nc", "r")

indexfile= addfile("/work/bm1147/b380949/EC-earth/dataset/anlysisneed/region_index_globalgrids.nc", "r")
area_ecearth= addfile("/work/bm1147/b380949/EC-earth/dataset/areas4_renamed.nc", "r")
area_ec= area_ecearth->area(0,:)
range_index_NAo= indexfile->range_index_NA
range_index_AMo=indexfile->range_index_amazon
range_index_CGo=indexfile->range_index_congo
range_index_NEo=indexfile->range_index_NE
range_index_EAo=indexfile->range_index_EA
range_index_SEAo=indexfile->range_index_SEA
range_index_AUo=indexfile->range_index_Australia
;print(range_index_AU)
;printVarSummary(range_index_SEA)
range_index_NAo=mask(range_index_NAo,range_index_NAo.gt.90000,False)
range_index_AMo=mask(range_index_AMo,range_index_AMo.gt.90000,False)
range_index_CGo=mask(range_index_CGo,range_index_CGo.gt.90000,False)
range_index_NEo=mask(range_index_NEo,range_index_NEo.gt.90000,False)
range_index_EAo=mask(range_index_EAo,range_index_EAo.gt.90000,False)
range_index_SEAo=mask(range_index_SEAo,range_index_SEAo.gt.90000,False)
range_index_AUo=mask(range_index_AUo,range_index_AUo.gt.90000,False)

igood_NA = ind(.not.ismissing(range_index_NAo))
range_index_NA = doubletoint(range_index_NAo(igood_NA))
igood_AM = ind(.not.ismissing(range_index_AMo))
range_index_AM = doubletoint(range_index_AMo(igood_AM))
igood_CG = ind(.not.ismissing(range_index_CGo))
range_index_CG = doubletoint(range_index_CGo(igood_CG))
igood_NE = ind(.not.ismissing(range_index_NEo))
range_index_NE = doubletoint(range_index_NEo(igood_NE))
igood_EA = ind(.not.ismissing(range_index_EAo))
range_index_EA = doubletoint(range_index_EAo(igood_EA))
igood_SEA = ind(.not.ismissing(range_index_SEAo))
range_index_SEA = doubletoint(range_index_SEAo(igood_SEA))
igood_AU = ind(.not.ismissing(range_index_AUo))
range_index_AU = doubletoint(range_index_AUo(igood_AU))

pool4eo=new((/4,88838/),double) 
pool4eo(0,:) = myfile13n->cVeg_nonlocal_withoutinter(0,:)
pool4eo(1,:) =myfile23->cSoil_nonlocal_withoutinter(0,:)
pool4eo(2,:) =myfile43->cLitter_nonlocal_withoutinter(0,:)
pool4eo(3,:) =myfile33->cLand_nonlocal_withoutinter(0,:)
;.............afforestation..........
flux4eo=new((/4,88838/),double) 
flux4eo(0,:) = myfile53->cVeg_nonlocal_withoutinter(0,:)
flux4eo(1,:) =myfile63->cSoil_nonlocal_withoutinter(0,:)
flux4eo(2,:) =myfile83->cLitter_nonlocal_withoutinter(0,:)
flux4eo(3,:) =myfile73->cLand_nonlocal_withoutinter(0,:)
;................irri...................
pool43eo=new((/4,88838/),double) 
pool43eo(0,:)= myfile133->cVeg_nonlocal_withoutinter(0,:)
pool43eo(1,:)=myfile143->cSoil_nonlocal_withoutinter(0,:)
pool43eo(2,:)= myfile153->cLitter_nonlocal_withoutinter(0,:)
pool43eo(3,:) =myfile163->cLand_nonlocal_withoutinter(0,:)
;.......deal with missing value...........
pool4eo=mask(pool4eo,isnan_ieee(pool4eo),False)
flux4eo=mask(flux4eo,isnan_ieee(flux4eo),False)
pool43eo=mask(pool43eo,isnan_ieee(pool43eo),False)
print(any(isnan_ieee(pool4eo)))
print(any(isnan_ieee(flux4eo)))
print(any(isnan_ieee(pool43eo)))
;.........total...............
pool5eo=new((/4,88838/),double) 
pool5eo(0,:) = myfile13n->cVeg_total_withoutinter(0,:)
pool5eo(1,:) =myfile23->cSoil_total_withoutinter(0,:)
pool5eo(2,:) =myfile43->cLitter_total_withoutinter(0,:)
pool5eo(3,:) =myfile33->cLand_total_withoutinter(0,:)
;.............afforestation..........
flux5eo=new((/4,88838/),double) 
flux5eo(0,:) = myfile53->cVeg_total_withoutinter(0,:)
flux5eo(1,:) =myfile63->cSoil_total_withoutinter(0,:)
flux5eo(2,:) =myfile83->cLitter_total_withoutinter(0,:)
flux5eo(3,:) =myfile73->cLand_total_withoutinter(0,:)
;................irri...................
pool53eo=new((/4,88838/),double) 
pool53eo(0,:)= myfile133->cVeg_total_withoutinter(0,:)
pool53eo(1,:)=myfile143->cSoil_total_withoutinter(0,:)
pool53eo(2,:)= myfile153->cLitter_total_withoutinter(0,:)
pool53eo(3,:) =myfile163->cLand_total_withoutinter(0,:)
;.......deal with missing value...........
pool5eo=mask(pool5eo,isnan_ieee(pool5eo),False)
flux5eo=mask(flux5eo,isnan_ieee(flux5eo),False)
pool53eo=mask(pool53eo,isnan_ieee(pool53eo),False)
print(any(isnan_ieee(pool5eo)))
print(any(isnan_ieee(flux5eo)))
print(any(isnan_ieee(pool53eo)))
;.......total:crop:pool5eo(4,88838),frst:flux5eo,irri:pool53eo .................................
;.......nonlocal:pool4eo,flux4eo,pool43eo
total_ecearth_NA=new((/3,4/),double)
total_ecearth_AM=new((/3,4/),double)
total_ecearth_CG=new((/3,4/),double)
total_ecearth_NE=new((/3,4/),double)
total_ecearth_EA=new((/3,4/),double)
total_ecearth_SEA=new((/3,4/),double)
total_ecearth_AU=new((/3,4/),double)
total_ecearth_g=new((/3,4/),double)

nonlocal_ecearth_NA=new((/3,4/),double)
nonlocal_ecearth_AM=new((/3,4/),double)
nonlocal_ecearth_CG=new((/3,4/),double)
nonlocal_ecearth_NE=new((/3,4/),double)
nonlocal_ecearth_EA=new((/3,4/),double)
nonlocal_ecearth_SEA=new((/3,4/),double)
nonlocal_ecearth_AU=new((/3,4/),double)
nonlocal_ecearth_g=new((/3,4/),double)

total_ecearth_g(2,:)=dim_sum_wgt_n_Wrap(pool53eo(:,:),area_ec(:), 1, 1)
nonlocal_ecearth_g(2,:)=dim_sum_wgt_n_Wrap(pool43eo(:,:),area_ec(:), 1, 1)
total_ecearth_NA(2,:)=dim_sum_wgt_n_Wrap(pool53eo(:,range_index_NA),area_ec(range_index_NA), 1, 1)
nonlocal_ecearth_NA(2,:)=dim_sum_wgt_n_Wrap(pool43eo(:,range_index_NA),area_ec(range_index_NA), 1, 1)
total_ecearth_AM(2,:)=dim_sum_wgt_n_Wrap(pool53eo(:,range_index_AM),area_ec(range_index_AM), 1, 1)
nonlocal_ecearth_AM(2,:)=dim_sum_wgt_n_Wrap(pool43eo(:,range_index_AM),area_ec(range_index_AM), 1, 1)
total_ecearth_CG(2,:)=dim_sum_wgt_n_Wrap(pool53eo(:,range_index_CG),area_ec(range_index_CG), 1, 1)
nonlocal_ecearth_CG(2,:)=dim_sum_wgt_n_Wrap(pool43eo(:,range_index_CG),area_ec(range_index_CG), 1, 1)
total_ecearth_NE(2,:)=dim_sum_wgt_n_Wrap(pool53eo(:,range_index_NE),area_ec(range_index_NE), 1, 1)
nonlocal_ecearth_NE(2,:)=dim_sum_wgt_n_Wrap(pool43eo(:,range_index_NE),area_ec(range_index_NE), 1, 1)
total_ecearth_EA(2,:)=dim_sum_wgt_n_Wrap(pool53eo(:,range_index_EA),area_ec(range_index_EA), 1, 1)
nonlocal_ecearth_EA(2,:)=dim_sum_wgt_n_Wrap(pool43eo(:,range_index_EA),area_ec(range_index_EA), 1, 1)
total_ecearth_SEA(2,:)=dim_sum_wgt_n_Wrap(pool53eo(:,range_index_SEA),area_ec(range_index_SEA), 1, 1)
nonlocal_ecearth_SEA(2,:)=dim_sum_wgt_n_Wrap(pool43eo(:,range_index_SEA),area_ec(range_index_SEA), 1, 1)
total_ecearth_AU(2,:)=dim_sum_wgt_n_Wrap(pool53eo(:,range_index_AU),area_ec(range_index_AU), 1, 1)
nonlocal_ecearth_AU(2,:)=dim_sum_wgt_n_Wrap(pool43eo(:,range_index_AU),area_ec(range_index_AU), 1, 1)

total_ecearth_g(0,:)=dim_sum_wgt_n_Wrap(pool5eo(:,:),area_ec(:), 1, 1)
nonlocal_ecearth_g(0,:)=dim_sum_wgt_n_Wrap(pool4eo(:,:),area_ec(:), 1, 1)
total_ecearth_NA(0,:)=dim_sum_wgt_n_Wrap(pool5eo(:,range_index_NA),area_ec(range_index_NA), 1, 1)
nonlocal_ecearth_NA(0,:)=dim_sum_wgt_n_Wrap(pool4eo(:,range_index_NA),area_ec(range_index_NA), 1, 1)
total_ecearth_AM(0,:)=dim_sum_wgt_n_Wrap(pool5eo(:,range_index_AM),area_ec(range_index_AM), 1, 1)
nonlocal_ecearth_AM(0,:)=dim_sum_wgt_n_Wrap(pool4eo(:,range_index_AM),area_ec(range_index_AM), 1, 1)
total_ecearth_CG(0,:)=dim_sum_wgt_n_Wrap(pool5eo(:,range_index_CG),area_ec(range_index_CG), 1, 1)
nonlocal_ecearth_CG(0,:)=dim_sum_wgt_n_Wrap(pool4eo(:,range_index_CG),area_ec(range_index_CG), 1, 1)
total_ecearth_NE(0,:)=dim_sum_wgt_n_Wrap(pool5eo(:,range_index_NE),area_ec(range_index_NE), 1, 1)
nonlocal_ecearth_NE(0,:)=dim_sum_wgt_n_Wrap(pool4eo(:,range_index_NE),area_ec(range_index_NE), 1, 1)
total_ecearth_EA(0,:)=dim_sum_wgt_n_Wrap(pool5eo(:,range_index_EA),area_ec(range_index_EA), 1, 1)
nonlocal_ecearth_EA(0,:)=dim_sum_wgt_n_Wrap(pool4eo(:,range_index_EA),area_ec(range_index_EA), 1, 1)
total_ecearth_SEA(0,:)=dim_sum_wgt_n_Wrap(pool5eo(:,range_index_SEA),area_ec(range_index_SEA), 1, 1)
nonlocal_ecearth_SEA(0,:)=dim_sum_wgt_n_Wrap(pool4eo(:,range_index_SEA),area_ec(range_index_SEA), 1, 1)
total_ecearth_AU(0,:)=dim_sum_wgt_n_Wrap(pool5eo(:,range_index_AU),area_ec(range_index_AU), 1, 1)
nonlocal_ecearth_AU(0,:)=dim_sum_wgt_n_Wrap(pool4eo(:,range_index_AU),area_ec(range_index_AU), 1, 1)

total_ecearth_g(1,:)=dim_sum_wgt_n_Wrap(flux5eo(:,:),area_ec(:), 1, 1)
nonlocal_ecearth_g(1,:)=dim_sum_wgt_n_Wrap(flux4eo(:,:),area_ec(:), 1, 1)
total_ecearth_NA(1,:)=dim_sum_wgt_n_Wrap(flux5eo(:,range_index_NA),area_ec(range_index_NA), 1, 1)
nonlocal_ecearth_NA(1,:)=dim_sum_wgt_n_Wrap(flux4eo(:,range_index_NA),area_ec(range_index_NA), 1, 1)
total_ecearth_AM(1,:)=dim_sum_wgt_n_Wrap(flux5eo(:,range_index_AM),area_ec(range_index_AM), 1, 1)
nonlocal_ecearth_AM(1,:)=dim_sum_wgt_n_Wrap(flux4eo(:,range_index_AM),area_ec(range_index_AM), 1, 1)
total_ecearth_CG(1,:)=dim_sum_wgt_n_Wrap(flux5eo(:,range_index_CG),area_ec(range_index_CG), 1, 1)
nonlocal_ecearth_CG(1,:)=dim_sum_wgt_n_Wrap(flux4eo(:,range_index_CG),area_ec(range_index_CG), 1, 1)
total_ecearth_NE(1,:)=dim_sum_wgt_n_Wrap(flux5eo(:,range_index_NE),area_ec(range_index_NE), 1, 1)
nonlocal_ecearth_NE(1,:)=dim_sum_wgt_n_Wrap(flux4eo(:,range_index_NE),area_ec(range_index_NE), 1, 1)
total_ecearth_EA(1,:)=dim_sum_wgt_n_Wrap(flux5eo(:,range_index_EA),area_ec(range_index_EA), 1, 1)
nonlocal_ecearth_EA(1,:)=dim_sum_wgt_n_Wrap(flux4eo(:,range_index_EA),area_ec(range_index_EA), 1, 1)
total_ecearth_SEA(1,:)=dim_sum_wgt_n_Wrap(flux5eo(:,range_index_SEA),area_ec(range_index_SEA), 1, 1)
nonlocal_ecearth_SEA(1,:)=dim_sum_wgt_n_Wrap(flux4eo(:,range_index_SEA),area_ec(range_index_SEA), 1, 1)
total_ecearth_AU(1,:)=dim_sum_wgt_n_Wrap(flux5eo(:,range_index_AU),area_ec(range_index_AU), 1, 1)
nonlocal_ecearth_AU(1,:)=dim_sum_wgt_n_Wrap(flux4eo(:,range_index_AU),area_ec(range_index_AU), 1, 1)
;..............................................................
total_ecearth_g!0="scenario"
nonlocal_ecearth_g!0="scenario"
;.region:North_America{28:54},{260:300} ;Amazon{-16:3},{283:308} ;Congo{-9:9},{13:29} 
;North Eurasia{56:67},{32:92};East Asia{46:62},{106:140};Southeast Asia {22:36},{101:119}
;Australia {-23:-12},{115:150}

ratio_ecearth=new((/8,3,4/),double);8 region;3:crop-ctl;frst-ctl;irri-crop;4:cVeg,cSoil,cLitter,cLand
do i=0,3
ratio_ecearth(0,:,i)=nonlocal_ecearth_g(:,i)/total_ecearth_g(:,i)
ratio_ecearth(1,:,i)=nonlocal_ecearth_NA(:,i)/total_ecearth_NA(:,i)
ratio_ecearth(2,:,i)=nonlocal_ecearth_AM(:,i)/total_ecearth_AM(:,i)
ratio_ecearth(3,:,i)=nonlocal_ecearth_CG(:,i)/total_ecearth_CG(:,i)
ratio_ecearth(4,:,i)=nonlocal_ecearth_NE(:,i)/total_ecearth_NE(:,i)
ratio_ecearth(5,:,i)=nonlocal_ecearth_EA(:,i)/total_ecearth_EA(:,i)
ratio_ecearth(6,:,i)=nonlocal_ecearth_SEA(:,i)/total_ecearth_SEA(:,i)
ratio_ecearth(7,:,i)=nonlocal_ecearth_AU(:,i)/total_ecearth_AU(:,i)
end do
;.........................................................
wks_type = "pdf"
wks_type@wkPaperWidthF =  8
wks_type@wkPaperHeightF = 8
wks = gsn_open_wks(wks_type,"nonlocal_ratio_soil_21oct")          ; send graphics to PNG file
;.....sres......................................................
     sres = True
     sres@vpWidthF = 0.6
     sres@vpHeightF = 0.3
     sres@gsnPaperMargin = 0
     ;sres@vpXF = 0.2
     ;sres@vpYF=0.9
     sres@trXMinF = 0.
     sres@trXMaxF = 10.75
     sres@gsnYRefLine=0.
     sres@gsnYRefLineColor="grey40"
     sres@gsnDraw = False
     sres@gsnFrame =False
     sres@gsnXYBarChart = True
     sres@gsnXYBarChartBarWidth = 0.25           ; change bar widths
     sres@gsnXYBarChartOutlineThicknessF = 0.5
     sres@tmXBMode          = "Explicit"         ; explicit labels
     sres@tmXBValues        = (/0.75,2,3.25,4.5,5.75,7,8.25,9.5/);fspan(0.75, 7.75, 1) ;(/0.75,1.75,2.75,3.75,4.75,5.75,6.75,7.75/)
     sres@tmXBLabels = (/" "," "," "," "," "," "," "," "/)
     sres@tiXAxisFontHeightF=0.02
     sres@tmYLMajorThicknessF = 0.5
    sres@tmXBMajorThicknessF = 0.5
sres@tmYLMinorOn =False
sres@tmXBMajorOutwardLengthF = 0.0
sres@tmYLMajorOutwardLengthF=0.0
sres@tmXBMajorLengthF=0.0
sres@tmYLMajorLengthF=0.01
sres@tmXTOn        = False ; turn off top tickmarks
sres@tmYROn        = False ; turn off right tickmarks
sres@tmXTBorderOn  = True ; turn off top axis
sres@tmYRBorderOn  = True ; turn off right axis
sres@tmBorderThicknessF= 1.0
     sres@tmXBLabelFontHeightF = 0.02
     sres@tmYLLabelFontHeightF = 0.02     
     sres@tiMainFontHeightF = 0.024
     sres@tiMainOffsetYF = -0.024
     sres@tiMainOffsetXF = 0.
     sres@tiMainFont = "times-roman"
     sres@tmXBLabelFont="times-roman"
     sres@tmYLLabelFont="times-roman"
     sres@tiXAxisFont ="times-roman"
     sres@tiYAxisFont="times-roman"
     sres@tiYAxisFontHeightF=0.024
sres@gsnLeftStringFontHeightF=0.024
sres@gsnLeftStringOrthogonalPosF=0.009
sres@gsnLeftStringParallelPosF= 0.0
sres@gsnStringFont="times-Bold"
sres@tmYLLabelDeltaF = -0.8   
     ;...............................................................
     plot=new((/3,7/),graphic) ;first epoch  second pre or wind
     plot1=new(3,graphic)     
     sres@trYMinF = -10.   ;min(ratio_ecearth(:,0,1)*100)-1
     sres@trYMaxF = max(ratio_mpi(:,0)*100)+1
    sres@gsnLeftString="(a)"
    sres@tiMainString = "Cropland expansion (CROP-CTL)"
    sres@tiYAxisString = "Ratio of non-local/total effect (%)"
    sres@tiYAxisOffsetYF= -0.05
    sres@tiYAxisOffsetXF= -0.004
     sres@gsnXYBarChartColors = (/0.51,0.68,0.50,0.5/)
     plot(0,0) = gsn_csm_xy(wks,(/0.625,1.875,3.125,4.375,5.625,6.875,8.125,9.375/),ratio_mpi(:,0)*100,sres)   ; draw each time series
     sres@gsnXYBarChartColors = (/0.49,0.64,0.82,1.0/)
     plot(1,0) = gsn_csm_xy(wks,(/0.875,2.125,3.375,4.625,5.875,7.125,8.375,9.625/),ratio_cesm(:,0)*100,sres)   ; advancing the frame
     sres@gsnXYBarChartColors = (/0.86,0.66,0.41,0.5/)
     plot(2,0) = gsn_csm_xy(wks,(/1.125,2.375,3.625,4.875,6.125,7.375,8.625,9.875/),ratio_ecearth(:,0,1)*100,sres)   ; advancing the frame
     overlay(plot(0,0),plot(1,0))
     overlay(plot(0,0),plot(2,0))

     sres@trYMinF = -60.         ;min(ratio_ecearth(:,1,1)*100)-1
     sres@trYMaxF = max(ratio_ecearth(:,1,1)*100)+1
      sres@tiYAxisString = " "
      sres@gsnLeftString="(b)"
      sres@tiMainOffsetXF = -0.022
      sres@tiMainString = "Afforestation (FRST-CTL)"
       sres@gsnXYBarChartColors = (/0.51,0.68,0.50,0.5/)
     plot(0,1) = gsn_csm_xy(wks,(/0.625,1.875,3.125,4.375,5.625,6.875,8.125,9.375/),ratio_mpi(:,1)*100,sres)   ; draw each time series
     sres@gsnXYBarChartColors = (/0.49,0.64,0.82,1.0/)
     plot(1,1) = gsn_csm_xy(wks,(/0.875,2.125,3.375,4.625,5.875,7.125,8.375,9.625/),ratio_cesm(:,1)*100,sres)   ; advancing the frame
     sres@gsnXYBarChartColors = (/0.86,0.66,0.41,0.5/)
     plot(2,1) = gsn_csm_xy(wks,(/1.125,2.375,3.625,4.875,6.125,7.375,8.625,9.875/),ratio_ecearth(:,1,1)*100,sres)   ; advancing the frame
     overlay(plot(0,1),plot(1,1))
     overlay(plot(0,1),plot(2,1))

    sres@trYMinF = 3240   ;   min(ratio_mpi(:,2)*100)-10
    sres@trYMaxF = 3260    ; max(ratio_mpi(:,2)*100)+10
    sres@tmYLMode = "Explicit"
    sres@tmYLValues  = (/3240,3240,3250,3260/)
    sres@tmYLLabels  = (/" ","3240","3250","3260"/)
     sres@vpHeightF = 0.1
      sres@gsnLeftString="(c)"
      sres@tiMainOffsetXF = 0
      sres@tiMainOffsetYF = -0.025
      sres@gsnLeftStringOrthogonalPosF=-0.05
      sres@tiMainString = "Irrigation expansion (IRR-CROP)"
    sres@tiYAxisString = " "
     sres@gsnXYBarChartColors =(/0.51,0.68,0.50,0.5/) ;(/"steelblue"/);MPIESM
     plot(0,2) = gsn_csm_xy(wks,(/0.625,1.875,3.125,4.375,5.625,6.875,8.125,9.375/),ratio_mpi(:,2)*100,sres)   ; draw each time series
     sres@gsnXYBarChartColors = (/0.49,0.64,0.82,1.0/);(/"orangered4"/) ;CESM
     plot(1,2) = gsn_csm_xy(wks,(/0.875,2.125,3.375,4.625,5.875,7.125,8.375,9.625/),ratio_cesm(:,2)*100,sres)   ; advancing the frame
     sres@gsnXYBarChartColors = (/0.86,0.66,0.41,0.5/);(/"sandybrown"/) ;EC-Earth
     plot(2,2) = gsn_csm_xy(wks,(/1.125,2.375,3.625,4.875,6.125,7.375,8.625,9.875/),ratio_ecearth(:,2,1)*100,sres)   ; advancing the frame
     overlay(plot(0,2),plot(1,2))
     overlay(plot(0,2),plot(2,2))

    sres@vpWidthF = 0.6
    sres@vpHeightF = 0.2
    sres@tmYLMode = "Automatic"
     sres@tmXBLabelAngleF = 23
     sres@tmXBLabels = (/"Global","ENA","WA","CCB","NE","NEA","SSEA","NAU"/)
    sres@tmXBMajorLengthF=0.01
    sres@gsnLeftString=" "
    sres@tiMainString = " "
    sres@tiYAxisString = " "
     sres@trYMinF = -110.   ;min(ratio_ecearth(:,0,1)*100)-1
     sres@trYMaxF = -30.
     sres@gsnXYBarChartColors = (/0.51,0.68,0.50,0.5/)
     plot(0,3) = gsn_csm_xy(wks,(/0.625,1.875,3.125,4.375,5.625,6.875,8.125,9.375/),ratio_mpi(:,0)*100,sres)   ; draw each time series
     sres@gsnXYBarChartColors = (/0.49,0.64,0.82,1.0/)
     plot(1,3) = gsn_csm_xy(wks,(/0.875,2.125,3.375,4.625,5.875,7.125,8.375,9.625/),ratio_cesm(:,0)*100,sres)   ; advancing the frame
     sres@gsnXYBarChartColors = (/0.86,0.66,0.41,0.5/)
     plot(2,3) = gsn_csm_xy(wks,(/1.125,2.375,3.625,4.875,6.125,7.375,8.625,9.875/),ratio_ecearth(:,0,1)*100,sres)   ; advancing the frame
     overlay(plot(0,3),plot(1,3))
     overlay(plot(0,3),plot(2,3))
     sres@tmXBLabelAngleF = 23
     sres@tmXBLabels = (/"Global","ENA","WA","CCB","NE","NEA","SSEA","NAU"/)
    sres@tmXBMajorLengthF=0.01
    sres@tmYLMajorLengthF=0.01
    sres@gsnLeftString=" "
    sres@tiMainString = " "
    sres@tiYAxisString = " "
     sres@trYMinF = -255.   ;min(ratio_ecearth(:,0,1)*100)-1
     sres@trYMaxF = -238.
    sres@tmYLMode = "Explicit"
    sres@tmYLValues  = (/-255,-250,-245,-240/)
    sres@tmYLLabels  = (/"-255","-250","-245","-240"/)
     sres@gsnXYBarChartColors = (/0.51,0.68,0.50,0.5/)
     plot(0,4) = gsn_csm_xy(wks,(/0.625,1.875,3.125,4.375,5.625,6.875,8.125,9.375/),ratio_mpi(:,1)*100,sres)   ; draw each time series
     sres@gsnXYBarChartColors = (/0.49,0.64,0.82,1.0/)
     plot(1,4) = gsn_csm_xy(wks,(/0.875,2.125,3.375,4.625,5.875,7.125,8.375,9.625/),ratio_cesm(:,1)*100,sres)   ; advancing the frame
     sres@gsnXYBarChartColors = (/0.86,0.66,0.41,0.5/)
     plot(2,4) = gsn_csm_xy(wks,(/1.125,2.375,3.625,4.875,6.125,7.375,8.625,9.875/),ratio_ecearth(:,1,1)*100,sres)   ; advancing the frame
     overlay(plot(0,4),plot(1,4))
     overlay(plot(0,4),plot(2,4))

    sres@trYMinF = -35   ;   min(ratio_mpi(:,2)*100)-10
    sres@trYMaxF = 135    ; max(ratio_mpi(:,2)*100)+10
     sres@tmYLMode = "Automatic"
     sres@vpHeightF = 0.2
      sres@gsnLeftString=" "
      sres@tiMainOffsetXF = 0
      sres@tiMainOffsetYF = -0.025
      sres@gsnLeftStringOrthogonalPosF=-0.008

      sres@tmXBLabels = (/" "," "," "," "," "," "," "," "/)
      sres@tmXBMajorLengthF=0.00
      sres@tiMainString = " "
    sres@tiYAxisString = " "
     sres@gsnXYBarChartColors =(/0.51,0.68,0.50,0.5/) ;(/"steelblue"/);MPIESM
     plot(0,5) = gsn_csm_xy(wks,(/0.625,1.875,3.125,4.375,5.625,6.875,8.125,9.375/),ratio_mpi(:,2)*100,sres)   ; draw each time series
     sres@gsnXYBarChartColors = (/0.49,0.64,0.82,1.0/);(/"orangered4"/) ;CESM
     plot(1,5) = gsn_csm_xy(wks,(/0.875,2.125,3.375,4.625,5.875,7.125,8.375,9.625/),ratio_cesm(:,2)*100,sres)   ; advancing the frame
     sres@gsnXYBarChartColors = (/0.86,0.66,0.41,0.5/);(/"sandybrown"/) ;EC-Earth
     plot(2,5) = gsn_csm_xy(wks,(/1.125,2.375,3.625,4.875,6.125,7.375,8.625,9.875/),ratio_ecearth(:,2,1)*100,sres)   ; advancing the frame
     overlay(plot(0,5),plot(1,5))
     overlay(plot(0,5),plot(2,5))

    sres@trYMinF = -700   ;   min(ratio_mpi(:,2)*100)-10
    sres@trYMaxF = -60   ; max(ratio_mpi(:,2)*100)+10
    sres@tmYLMode = "Explicit"
    sres@tmYLValues  = (/-100,-300,-600,-600/)
    sres@tmYLLabels  = (/"-100","-300","-600"," "/)
     sres@vpHeightF = 0.189
      sres@gsnLeftString=" "
      sres@tiMainOffsetXF = 0
      sres@tiMainOffsetYF = -0.025
      sres@gsnLeftStringOrthogonalPosF=-0.008
      sres@tmXBMajorLengthF=0.01
      sres@tmXBLabelAngleF = 23
     sres@tmXBLabels = (/"Global","ENA","WA","CCB","NE","NEA","SSEA","NAU"/)
      sres@tiMainString = " "
    sres@tiYAxisString = " "
     sres@gsnXYBarChartColors =(/0.51,0.68,0.50,0.5/) ;(/"steelblue"/);MPIESM
     plot(0,6) = gsn_csm_xy(wks,(/0.625,1.875,3.125,4.375,5.625,6.875,8.125,9.375/),ratio_mpi(:,2)*100,sres)   ; draw each time series
     sres@gsnXYBarChartColors = (/0.49,0.64,0.82,1.0/);(/"orangered4"/) ;CESM
     plot(1,6) = gsn_csm_xy(wks,(/0.875,2.125,3.375,4.625,5.875,7.125,8.375,9.625/),ratio_cesm(:,2)*100,sres)   ; advancing the frame
     sres@gsnXYBarChartColors = (/0.86,0.66,0.41,0.5/);(/"sandybrown"/) ;EC-Earth
     plot(2,6) = gsn_csm_xy(wks,(/1.125,2.375,3.625,4.875,6.125,7.375,8.625,9.875/),ratio_ecearth(:,2,1)*100,sres)   ; advancing the frame
     overlay(plot(0,6),plot(1,6))
     overlay(plot(0,6),plot(2,6))

  getvalues plot(0,0)                      ; Get plot size for use in
    "vpHeightF" : vph                 ; determining size of
    "vpWidthF"  : vpw                 ; labelbar.
  end getvalues   

  lbres                           = True
  lbres@lbAutoManage              = False         ; we control label bar
  lbres@lbFillColors              = (/(/0.86,0.66,0.41,1.0/),(/0.49,0.64,0.82,1.0/),(/0.51,0.68,0.50,1.0/)/)   ; use nice strided colors  
  lbres@lbPerimOn                 = False         ; Turn off labelbar perimeter.
  lbres@lbMonoFillPattern         = True          ; one pattern, all solid
  lbres@vpWidthF           = 0.25* vpw   ; labelbar width
  lbres@vpHeightF          = 0.2*vph    ; labelbar height
  lbres@lbBoxMajorExtentF  = 0.3          ; puts space between color boxes
 lbres@lbLabelFont = "times-roman"
  lbres@lbLabelFontHeightF = 0.01       ; font height. default is small
  lbres@lbLabelJust        = "CenterLeft"  ; left justify labels
  lbres@lgPerimColor       = "white"
lables=(/" MPI_ESM","   CESM","EC-Earth"/)
colors       = (/(/0.51,0.68,0.50,0.5/),(/0.49,0.64,0.82,1.0/),(/0.86,0.66,0.41,0.5/)/)
gsn_labelbar_ndc(wks,3,(/"EC-Earth","CESM","MPI-ESM"/),0.18,0.4,lbres)
     pnlres                = True
     pnlres@gsnMaximize    = True   ; Maximize in frame
     pnlres@gsnPaperMargin =0 
  pnlres@gsnBoxMargin  = 0
pnlres@gsnPanelRight = 1.0
pnlres@gsnPanelLeft =0.01
     pnlres@gsnPanelYWhiteSpacePercent = 0.00
     pnlres@gsnPanelXWhiteSpacePercent = 0.00
     pnlres@gsnPanelFigureStringsFontHeightF=0.02
     pnlres@amJust="TopLeft"
     pnlres@gsnPanelFigureStringsPerimOn      = False  
     ;drawNDCGrid(wks)
     pnlres@gsnPanelDebug = True
    pnlres@gsnPanelXF = (/0.045,0.372,0.7,0.045,0.372,0.7,0.7/) ;0.4
    pnlres@gsnPanelYF = (/0.57,0.57,0.57,0.422,0.422,0.5178,0.418/) 
      pnlres@gsnPanelRowSpec = True                   ; tell panel what order to plot
     gsn_panel(wks,plot(0,:),(/3,3,1/),pnlres)
end

