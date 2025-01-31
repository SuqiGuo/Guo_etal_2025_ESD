load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin
;..........input data......................................................
mpiesm = addfile("/home/b/b380949/gridarea/gridarea_sea0.nc","r")
cesm=addfile("/home/b/b380949/gridarea/AREA_CESM.nc","r")
ecearth=addfile("/work/bm1147/b380949/EC-earth/dataset/anlysisneed/gridarea_ecearth.nc","r")
;...................mpiesm..................................................
myfile1 = addfile("/work/bm1147/b380949/web-monitoring/crop-ctl/Eyr/cVeg/cVeg_crop-ctl_mpiesm_signal-separated_revised.nc", "r")
myfile2 = addfile("/work/bm1147/b380949/web-monitoring/crop-ctl/Eyr/cSoil/cSoil_crop-ctl_mpiesm_signal-separated_revised.nc", "r")
myfile4 = addfile("/work/bm1147/b380949/web-monitoring/crop-ctl/Eyr/cLitter/cLitter_crop-ctl_mpiesm_signal-separated_revised.nc", "r")
myfile3 = addfile("/work/bm1147/b380949/web-monitoring/crop-ctl/Eyr/cProduct/cProduct_crop-ctl_mpiesm_signal-separated_revised.nc", "r")

myfile5 = addfile("/work/bm1147/b380949/web-monitoring/frst-ctl/Eyr/cVeg/cVeg_frst-ctl_mpiesm_signal-separated_revised.nc", "r")
myfile6 = addfile("/work/bm1147/b380949/web-monitoring/frst-ctl/Eyr/cSoil/cSoil_frst-ctl_mpiesm_signal-separated_revised.nc", "r")
myfile8 = addfile("/work/bm1147/b380949/web-monitoring/frst-ctl/Eyr/cLitter/cLitter_frst-ctl_mpiesm_signal-separated_revised.nc", "r")
myfile7 = addfile("/work/bm1147/b380949/web-monitoring/frst-ctl/Eyr/cProduct/cProduct_frst-ctl_mpiesm_signal-separated_revised.nc", "r")

myfile13 = addfile("/work/bm1147/b380949/web-monitoring/irri-crop/Eyr/cVeg/cVeg_irri-crop_mpiesm_signal-separated_revised.nc", "r")
myfile14= addfile("/work/bm1147/b380949/web-monitoring/irri-crop/Eyr/cSoil/cSoil_irri-crop_mpiesm_signal-separated_revised.nc", "r")
myfile15= addfile("/work/bm1147/b380949/web-monitoring/irri-crop/Eyr/cLitter/cLitter_irri-crop_mpiesm_signal-separated_revised.nc", "r")
myfile16= addfile("/work/bm1147/b380949/web-monitoring/irri-crop/Eyr/cProduct/cProduct_irri-crop_mpiesm_signal-separated_revised.nc", "r")
;..........CESM................................................................
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
;.........EC-Earth...........................................
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
;ecearth--sea is missing value
cesmlsm=addfile("/home/b/b380949/gridarea/landmask_glacier_cesm.nc", "r")
slm = cesmlsm->landmask

ti=cd_calendar(myfile1->time,-1)
ti2=cd_calendar(myfile5->time,-1)
ti4=cd_calendar(myfile13->time,-1)
;........cell area............................
gridarea = mpiesm->slm ;mpiesm
gridarea1= cesm->area ;cesm
gridarea2=todouble(gridarea1)
copy_VarMeta(gridarea1, gridarea2)
gridarea0=where(ismissing(gridarea2), 0, gridarea2)
copy_VarMeta(gridarea2, gridarea0)
;print(gridarea0)
gridarea2n=gridarea0*1000000 ;km^2->m^2
copy_VarMeta(gridarea2, gridarea2n)

gridarea_ecearth=ecearth->cell_area ;ecearth
printVarSummary(gridarea_ecearth)
;................MPI
;.........deforestation..... ..................
pool=new((/4,160,96,192/),double) 
pool(0,:,:,:) = myfile1->cVeg_nonlocal(0:159,:,:)
pool&time=ti4
pool(1,:,:,:) =myfile2->cSoil_nonlocal(0:159,:,:)
pool(2,:,:,:) =myfile4->cLitter_nonlocal(0:159,:,:)
pool(3,:,:,:) =myfile3->cProduct_nonlocal(0:159,:,:)
;.............afforestation..........
pool_mpiesm_frst=new((/4,160,96,192/),double) 
pool_mpiesm_frst(0,:,:,:) = myfile5->cVeg_nonlocal(0:159,:,:)
pool_mpiesm_frst&time=ti4
pool_mpiesm_frst(1,:,:,:) =myfile6->cSoil_nonlocal(0:159,:,:)
pool_mpiesm_frst(2,:,:,:) =myfile8->cLitter_nonlocal(0:159,:,:)
pool_mpiesm_frst(3,:,:,:) =myfile7->cProduct_nonlocal(0:159,:,:)
;print(pool_mpiesm_frst(0,0,159,:,:))
;................irri...................
pool_mpiesm_irr=new((/4,160,96,192/),double) 
pool_mpiesm_irr(0,:,:,:)= myfile13->cVeg_nonlocal
pool_mpiesm_irr&time=ti4
pool_mpiesm_irr(1,:,:,:)=myfile14->cSoil_nonlocal
pool_mpiesm_irr(2,:,:,:)= myfile15->cLitter_nonlocal
pool_mpiesm_irr(3,:,:,:) =myfile16->cProduct_nonlocal
;................CESM
pool_cesm_crop=new((/4,160,192,288/),double) ;crop
pool_cesm_crop(0,:,:,:) = month_to_annual(myfile12->TOTVEGC_nonlocal(0:1919,:,:),1)
pool_cesm_crop(1,:,:,:) =month_to_annual(myfile22->TOTSOMC_nonlocal(0:1919,:,:),1)
pool_cesm_crop(2,:,:,:) =month_to_annual(myfile42->TOTLITC_nonlocal(0:1919,:,:),1)
pool_cesm_crop(3,:,:,:) =month_to_annual(myfile32->TOTECOSYSC_nonlocal(0:1919,:,:),1)  ;cLand 

pool_cesm_crop&year=ti4
;.............afforestation..........
pool_cesm_frst=new((/4,160,192,288/),double) ;frst
pool_cesm_frst(0,:,:,:) = month_to_annual(myfile52->TOTVEGC_nonlocal,1)
pool_cesm_frst&year=ti4
pool_cesm_frst(1,:,:,:) =month_to_annual(myfile62->TOTSOMC_nonlocal,1)
pool_cesm_frst(2,:,:,:) =month_to_annual(myfile82->TOTLITC_nonlocal,1)
pool_cesm_frst(3,:,:,:) =month_to_annual(myfile72->TOTECOSYSC_nonlocal,1)
;print(pool_cesm_frst(0,0,159,:,:))
;................irri...................
pool_cesm_irr=new((/4,160,192,288/),double) ;irri
pool_cesm_irr(0,:,:,:)= month_to_annual(myfile132->TOTVEGC_nonlocal(0:1919,:,:),1)
pool_cesm_irr&year=ti4
pool_cesm_irr(1,:,:,:)=month_to_annual(myfile142->TOTSOMC_nonlocal(0:1919,:,:),1)
pool_cesm_irr(2,:,:,:)= month_to_annual(myfile152->TOTLITC_nonlocal(0:1919,:,:),1)
pool_cesm_irr(3,:,:,:)= month_to_annual(myfile162->TOTECOSYSC_nonlocal(0:1919,:,:),1)

pool_cesm_crop_n=mask(pool_cesm_crop,  slm.eq.0, False)
copy_VarMeta(pool_cesm_crop, pool_cesm_crop_n)

pool_cesm_irr_n=mask(pool_cesm_irr, slm.eq.0, False)
copy_VarMeta(pool_cesm_irr, pool_cesm_irr_n)

pool_cesm_frst_n=mask(pool_cesm_frst, slm.eq.0, False)
copy_VarMeta(pool_cesm_frst, pool_cesm_frst_n)

;...ecearth............
;..crop........
pool_ecearth_crop=new((/4,160,192,288/),double) 
pool_ecearth_crop(0,:,:,:) = myfile13n->cVeg_nonlocal
pool_ecearth_crop&time=ti4
pool_ecearth_crop(1,:,:,:) =myfile23->cSoil_nonlocal
pool_ecearth_crop(2,:,:,:) =myfile43->cLitter_nonlocal
pool_ecearth_crop(3,:,:,:) =myfile33->cLand_nonlocal
printVarSummary(pool_ecearth_crop)
;.............afforestation..........
pool_ecearth_frst=new((/4,160,192,288/),double) 
pool_ecearth_frst(0,:,:,:) = myfile53->cVeg_nonlocal
pool_ecearth_frst&time=ti4
pool_ecearth_frst(1,:,:,:) =myfile63->cSoil_nonlocal
pool_ecearth_frst(2,:,:,:) =myfile83->cLitter_nonlocal
pool_ecearth_frst(3,:,:,:) =myfile73->cLand_nonlocal
;print(pool_ecearth_frst(0,0,159,:,:))
;................irri...................
pool_ecearth_irr=new((/4,160,192,288/),double) 
pool_ecearth_irr(0,:,:,:)= myfile133->cVeg_nonlocal
pool_ecearth_irr&time=ti4
pool_ecearth_irr(1,:,:,:)=myfile143->cSoil_nonlocal
pool_ecearth_irr(2,:,:,:)= myfile153->cLitter_nonlocal
pool_ecearth_irr(3,:,:,:) =myfile163->cLand_nonlocal
;.................................
pool_ecearth_crop=mask(pool_ecearth_crop,isnan_ieee(pool_ecearth_crop),False)
pool_ecearth_frst=mask(pool_ecearth_frst,isnan_ieee(pool_ecearth_frst),False)
pool_ecearth_irr=mask(pool_ecearth_irr,isnan_ieee(pool_ecearth_irr),False)
;.............global integral.....mpi...........................
Gpool=new((/5,160/),double) 
Gpool_mpiesm_frst=new((/5,160/),double) 
Gpool_mpiesm_irr=new((/5,160/),double) 
do i=0,3
Gpool(i,:)=wgt_areasum2(pool(i,:,:,:),gridarea,0)*0.000000000001  ;kg/m^2->GtC
Gpool_mpiesm_frst(i,:)=wgt_areasum2(pool_mpiesm_frst(i,:,:,:),gridarea,0)*0.000000000001  ;kg/m^2->GtC since gridarerevised unit2 km
Gpool_mpiesm_irr(i,:)=wgt_areasum2(pool_mpiesm_irr(i,:,:,:),gridarea,0)*0.000000000001  ;kg/m^2->GtC
end do
Gpool!0="carbonpool"
Gpool!1="time"
Gpool&time=ti4
Gpool(3,:)=dim_sum_n_Wrap(Gpool(0:3,:),0)

Gpool_mpiesm_frst!0="carbonpool"
Gpool_mpiesm_frst!1="time"
Gpool_mpiesm_frst&time=ti4
Gpool_mpiesm_frst(3,:)=dim_sum_n_Wrap(Gpool_mpiesm_frst(0:3,:),0)

Gpool_mpiesm_irr!0="carbonpool"
Gpool_mpiesm_irr!1="time"
Gpool_mpiesm_irr&time=ti4
Gpool_mpiesm_irr(3,:)=dim_sum_n_Wrap(Gpool_mpiesm_irr(0:3,:),0)
;...........CESM.........
Gpool_cesm_crop_n=new((/5,160/),double) 
Gpool_cesm_frst_n=new((/5,160/),double) 
Gpool_cesm_irr_n=new((/5,160/),double) 
do i=0,3
Gpool_cesm_crop_n(i,:)=wgt_areasum2(pool_cesm_crop_n(i,:,:,:),gridarea2n,0)*0.000000000000001  ;kg/m^2->GtC
Gpool_cesm_frst_n(i,:)=wgt_areasum2(pool_cesm_frst_n(i,:,:,:),gridarea2n,0)*0.000000000000001  ;kg/m^2->GtC since gridarerevised unit2 km
Gpool_cesm_irr_n(i,:)=wgt_areasum2(pool_cesm_irr_n(i,:,:,:),gridarea2n,0)*0.000000000000001  ;kg/m^2->GtC
end do
Gpool_cesm_crop_n!0="carbonpool"
Gpool_cesm_crop_n!1="time"
Gpool_cesm_crop_n&time=ti4

Gpool_cesm_frst_n!0="carbonpool"
Gpool_cesm_frst_n!1="time"
Gpool_cesm_frst_n&time=ti4

Gpool_cesm_irr_n!0="carbonpool"
Gpool_cesm_irr_n!1="time"
Gpool_cesm_irr_n&time=ti4
printVarSummary(Gpool_cesm_crop_n)
printMinMax(Gpool_cesm_crop_n, 1)
;..ecearth.....
Gpool_ecearth_crop=new((/5,160/),double) 
Gpool_ecearth_frst=new((/5,160/),double) 
Gpool_ecearth_irr=new((/5,160/),double) 
do i=0,3
Gpool_ecearth_crop(i,:)=wgt_areasum2(pool_ecearth_crop(i,:,:,:),gridarea2n,0)*0.000000000001  ;kg/m^2->GtC
Gpool_ecearth_frst(i,:)=wgt_areasum2(pool_ecearth_frst(i,:,:,:),gridarea2n,0)*0.000000000001  ;kg/m^2->GtC since gridarerevised unit2 km
Gpool_ecearth_irr(i,:)=wgt_areasum2(pool_ecearth_irr(i,:,:,:),gridarea2n,0)*0.000000000001  ;kg/m^2->GtC
end do
Gpool_ecearth_crop!0="carbonpool"
Gpool_ecearth_crop!1="time"
Gpool_ecearth_crop&time=ti4

Gpool_ecearth_frst!0="carbonpool"
Gpool_ecearth_frst!1="time"
Gpool_ecearth_frst&time=ti4

Gpool_ecearth_irr!0="carbonpool"
Gpool_ecearth_irr!1="time"
Gpool_ecearth_irr&time=ti4
printVarSummary(Gpool_ecearth_crop)
printMinMax(Gpool_ecearth_crop, 1)
printVarSummary(Gpool)
   fname = "cropland_expansion_timeseries.nc"
   system("/bin/rm -f "+"/work/bm1147/b380949/paper1-figure-scripts/"+fname)
   a = addfile("/work/bm1147/b380949/paper1-figure-scripts/"+fname,"c")
   a->Gpool = Gpool
   a->Gpool_cesm_crop_n = Gpool_cesm_crop_n
   a->Gpool_ecearth_crop = Gpool_ecearth_crop
   a->carbonpool = (/',xmax,xmax,xmin,xmin/)

  fname2 = "Afforestation_timeseries.nc"
   system("/bin/rm -f "+"/work/bm1147/b380949/paper1-figure-scripts/"+fname2)
   b = addfile("/work/bm1147/b380949/paper1-figure-scripts/"+fname2,"c")
   b->Gpool_mpiesm_frst = Gpool_mpiesm_frst
   b->Gpool_cesm_frst_n = Gpool_cesm_frst_n
   b->Gpool_ecearth_frst = Gpool_ecearth_frst

   fname3 = "Irrigation_timeseries.nc"
   system("/bin/rm -f "+"/work/bm1147/b380949/paper1-figure-scripts/"+fname3)
   c = addfile("/work/bm1147/b380949/paper1-figure-scripts/"+fname3,"c")
   c->Gpool_mpiesm_irr = Gpool_mpiesm_irr
   c->Gpool_cesm_irr_n = Gpool_cesm_irr_n
   c->Gpool_ecearth_irr = Gpool_ecearth_irr
end