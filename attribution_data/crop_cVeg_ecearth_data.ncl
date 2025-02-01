
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin
;......................................................................
myfile = addfile("/home/b/b380949/gridarea/gridarea_sea0.nc","r")
my=addfile("/home/b/b380949/gridarea/AREA_CESM.nc","r")
myfilen = addfile("/work/bm1147/b380949/web-monitoring/irri-crop/Eyr/cVeg/cVeg_irri-crop_mpiesm_signal-separated_revised.nc", "r")

myfile1 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/crop-ctl/Eyr/cVeg/interped_cVeg_crop-ctl_ecearth_signal-separated.nc", "r")
myfile2 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/crop-ctl/Eyr/cSoil/interped_cSoil_crop-ctl_ecearth_signal-separated.nc", "r")

myfile5 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/frst-ctl/Eyr/cVeg/interped_cVeg_frst-ctl_ecearth_signal-separated.nc", "r")
myfile6 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/frst-ctl/Eyr/cSoil/interped_cSoil_frst-ctl_ecearth_signal-separated.nc", "r")

myfile13 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/irr-crop/Eyr/cVeg/interped_cVeg_irr-crop_ecearth_signal-separated.nc", "r")
myfile14 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/irr-crop/Eyr/cSoil/interped_cSoil_irr-crop_ecearth_signal-separated.nc", "r")

myfile1_2 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/crop-ctl/Eyr/mrsos/interped_mrsos_crop-ctl_ecearth_signal-separated.nc", "r")
;ecearth mrsos unit:m^3/m^3-->*10kg/m^2i.e.1m^3/m^3=10kg/m^2
;kg m-2 mrsos
;(time,depth,lat,lon) kg m-2 mrsol
myfile1_4 = addfile("/work/bm1147/b380949/web-monitoring/EC-earth/crop-ctl/Eyr/tas/interped_tas_crop-ctl_ecearth_signal-separated.nc", "r")
;ecearth tas unit: K
;myfile1_2 = addfile("/work/bm1147/b380949/web-monitoring/CESM/frst-ctl/Lmon/SOILWATER_10CM/SOILWATER_10CM_frst-ctl_cesm_signal-separated.nc", "r")
;myfile1_4 = addfile("/work/bm1147/b380949/web-monitoring/CESM/frst-ctl/Lmon/TREFHT/TREFHT_frst-ctl_cesm_signal-separated.nc", "r")
;myfile1_2 = addfile("/work/bm1147/b380949/web-monitoring/CESM/irri-crop/Lmon/SOILWATER_10CM/SOILWATER_10CM_irri-crop_cesm_signal-separated.nc", "r")

ti4=cd_calendar(myfilen->time(10::),-1)
;........read in data............................
gridarea = myfile->slm
landsea_ecearth=addfile("/work/bm1147/b380949/EC-earth/dataset/anlysisneed/landmask_ecearth_regridded.nc","r")
landmask_ecearth=landsea_ecearth->landmask(0,:,:)

;................CESM
poolo=new((/4,150,192,288/),double) ;crop
poolo(0,:,:,:) = myfile1->cVeg_nonlocal(10:159,:,:)
poolo(1,:,:,:) =myfile2->cSoil_nonlocal(10:159,:,:)
printVarSummary(poolo)
poolo&time=ti4
;.............afforestation..........
flux=new((/4,150,192,288/),double) ;frst
flux(0,:,:,:) = myfile5->cVeg_nonlocal(10:159,:,:)
flux&time=ti4
flux(1,:,:,:) =myfile6->cSoil_nonlocal(10:159,:,:)
;print(flux(0,0,159,:,:))
;................irri...................
pool3=new((/4,150,192,288/),double) ;irri
pool3(0,:,:,:)= myfile13->cVeg_nonlocal(10:159,:,:)
pool3&time=ti4
pool3(1,:,:,:)=myfile14->cSoil_nonlocal(10:159,:,:)
printVarSummary(poolo(0,1:149,{-60:80},:))
printVarSummary(pool3)
;.............cVeg_nonlocal_crop_expansion...............
poolo=mask(poolo,isnan_ieee(poolo),False)
print(any(isnan_ieee(poolo)))
pool=mask(poolo, landmask_ecearth.eq.1, True)
copy_VarMeta(poolo, pool);......................................
Gpool_cveg_diff=new((/149,149,288/),double) 
Gpool_cveg_diff=pool(0,1:149,{-60:80},:)-pool(0,0:148,{-60:80},:)
copy_VarMeta(pool(0,1:149,{-60:80},:), Gpool_cveg_diff)
Gpool_cveg_diff&time=ti4(1:149)
Gpool_cveg_diff_avg=dim_avg_n_Wrap(Gpool_cveg_diff, 0)
;printVarSummary(Gpool_cveg_diff)
;printVarSummary(Gpool_cveg_diff)
;..............climate........................................
tas=new((/149,149,288/),double) ;3varia, 3 simulation
mrsosn=new((/149,149,288/),double) ;
tas= myfile1_4->tas_nonlocal(1::,{-60:80},:)
mrsosn = myfile1_2->mrsos_nonlocal(1::,{-60:80},:)
mrsos=mrsosn*10
copy_VarMeta(mrsosn, mrsos)

;....when use pr as a independent variable 
;mrsos_new = myfile1_2->pr_nonlocal(12::,{-60:80},:)
;mrsos=mrsos_new*86400 ;/day
;copy_VarMeta(mrsos_new, mrsos)
;..........................................
variables=new((/3,149,149,288/),double) 
variables(0,:,:,:)=1
variables(1,:,:,:)=tas
variables(2,:,:,:)=mrsos
Var_avg=dim_avg_n_Wrap(variables,1)
bstd=new((/3,149,288/),double)  ;standardized b
b=new((/3,149,288/),double)  ;standardized b
pval=new((/3,149,288/),string)   ;the reliability of t test result of b
Fr=new((/149,288/),double)  ;f test result
F_pval=new((/5,192/),string) ; the reliability of the whole regression
index=new(149,long)
a=new((/3,149/),double)
  do i=0,148
  do j=0,287
;.......................
    index=0
;........................................................
index=ind((.not.ismissing(Gpool_cveg_diff(:,i,j))).and.(.not.ismissing(variables(1,:,i,j)))\
          .and.(.not.ismissing(variables(2,:,i,j))))
if (.not.all(ismissing(index)))   ;???
;..........................................
b(:,i,j)=reg_multlin(Gpool_cveg_diff(index,i,j),variables(:,index,i,j),False) ; partial regression coef
end if
end do 
end do
copy_VarMeta(Var_avg, b)
deltacveg_time=new((/3,149,149,288/),double)  ;standardized b
do i= 0,148
deltacveg_time(:,i,:,:)=b(:,:,:)*variables(:,i,:,:) ;1 t2m,2 mrso
end do
copy_VarMeta(variables,deltacveg_time)
deltacveg_time_new=mask(deltacveg_time, deltacveg_time.eq.0, False)
;;;;;;;;;;;;;;;Ftest..........
prediction=dim_sum_n_Wrap(deltacveg_time_new, 0);3,149,149,288
syy=dim_variance_n_Wrap(Gpool_cveg_diff, 0); Gpool_cveg_diff=new((/149,149,288/),double) 
u_inter=new((/149,149,288/),double)  ;standardized b
do i=0,148
  u_inter(i,:,:)=(prediction(i,:,:)-Gpool_cveg_diff_avg)^2
end do 
copy_VarMeta(prediction, u_inter)
u=dim_sum_n_Wrap(u_inter, 0); 4,149,288
r2=new((/149,288/),double)  ;standardized b
F=new((/149,288/),double)  ;standardized b
syy_new=mask(syy, syy.eq.0, False)
copy_VarMeta(syy, syy_new)

r2(:,:)=u(:,:)/syy_new
copy_VarMeta(u, r2)

inter=1-r2
copy_VarMeta(r2, inter)
inter_new=mask(inter, inter.eq.0, False)
copy_VarMeta(inter, inter_new)

F=(79*r2)/inter_new
copy_VarMeta(syy, F)
Q=new((/149,288/),double)  ;standardized b
Q=syy-u

copy_VarMeta(syy, Q)
;...........correlation between predication and Gpool_cveg_diff=r2
correlation=new((/149,288/),double)  ;standardized b
correlation=escorc_n(prediction, Gpool_cveg_diff, 0, 0)
copy_VarMeta(syy, correlation)
;.............data output....................
b@long_name="Regression coefficients"
correlation@long_name="correlation of prediction and original dataset(R-squared)"

diro = "/home/b/b380949/mulreg/"                     ; Output directory
        filo ="attri_crop_cVeg_ecearth.nc"             ; Output file
  system("/bin/rm -f " + diro + filo)    ; remove if exists
  fout  = addfile (diro + filo, "c")  ; open output file

   fileAtt                           = True 
   fileAtt@creation_date = systemfunc("date")
   fileAtt@contact = "Suqi.Guo@geographie.uni-muenchen.de ;Ludwig-Maximilians-University Munich"
   fileAtt@comment = "Multiple linear regression analysis between carbon pools and climate(2m-temperature and surface soil moisture);~C~The first dimention of b indicate predictors; b=0-->constant term; b=1-->2m-temperature; b=2-->surface soil moisture;~C~Correlation is R-squared which is the correlation of prediction and original dataset.~C~It also the percentage of the dependent variable variation the regression model explains."
   fileattdef (fout, fileAtt)
   
   setfileoption(fout,"DefineMode",False) 
      fout->b  =  b
      fout->correlation  = correlation

end

