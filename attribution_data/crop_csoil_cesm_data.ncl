
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin
;......................................................................
myfile = addfile("/pf/b/b380949/gridarea/gridarea_sea0.nc","r")
my=addfile("/pf/b/b380949/gridarea/AREA_CESM.nc","r")
myfilen = addfile("/work/bm1147/b380949/web-monitoring/irri-crop/Eyr/cVeg/cVeg_irri-crop_mpiesm_signal-separated_revised.nc", "r")

myfile1 = addfile("/work/bm1147/b380949/web-monitoring/CESM/crop-ctl/Lmon/TOTVEGC/TOTVEGC_crop-ctl_cesm_signal-separated.nc", "r")
myfile2 = addfile("/work/bm1147/b380949/web-monitoring/CESM/crop-ctl/Lmon/TOTSOMC/TOTSOMC_crop-ctl_cesm_signal-separated.nc", "r")

myfile5 = addfile("/work/bm1147/b380949/web-monitoring/CESM/frst-ctl/Lmon/TOTVEGC/TOTVEGC_frst-ctl_cesm_signal-separated.nc", "r")
myfile6 = addfile("/work/bm1147/b380949/web-monitoring/CESM/frst-ctl/Lmon/TOTSOMC/TOTSOMC_frst-ctl_cesm_signal-separated.nc", "r")

myfile13 = addfile("/work/bm1147/b380949/web-monitoring/CESM/irri-crop/Lmon/TOTVEGC/TOTVEGC_irri-crop_cesm_signal-separated.nc", "r")
myfile14 = addfile("/work/bm1147/b380949/web-monitoring/CESM/irri-crop/Lmon/TOTSOMC/TOTSOMC_irri-crop_cesm_signal-separated.nc", "r")
myfile1_2 = addfile("/work/bm1147/b380949/web-monitoring/CESM/crop-ctl/Lmon/SOILWATER_10CM/SOILWATER_10CM_crop-ctl_cesm_signal-separated.nc", "r")
;;kg m-2 mrsos
;(time,depth,lat,lon) kg m-2 mrsol
myfile1_4 = addfile("/work/bm1147/b380949/web-monitoring/CESM/crop-ctl/Lmon/TREFHT/TREFHT_crop-ctl_cesm_signal-separated.nc", "r")
;myfile1_2 = addfile("/work/bm1147/b380949/web-monitoring/CESM/frst-ctl/Lmon/SOILWATER_10CM/SOILWATER_10CM_frst-ctl_cesm_signal-separated.nc", "r")
;myfile1_4 = addfile("/work/bm1147/b380949/web-monitoring/CESM/frst-ctl/Lmon/TREFHT/TREFHT_frst-ctl_cesm_signal-separated.nc", "r")
;myfile1_2 = addfile("/work/bm1147/b380949/web-monitoring/CESM/irri-crop/Lmon/SOILWATER_10CM/SOILWATER_10CM_irri-crop_cesm_signal-separated.nc", "r")

ti4=cd_calendar(myfilen->time(10::),-1)
;........read in data............................
gridarea = myfile->slm
;................CESM
pooln=new((/4,150,192,288/),double) ;crop
pooln(0,:,:,:) = month_to_annual(myfile1->TOTVEGC_nonlocal(120:1919,:,:),1)
pooln(1,:,:,:) =month_to_annual(myfile2->TOTSOMC_nonlocal(120:1919,:,:),1)
printVarSummary(pooln)
pooln&year=ti4
;.............afforestation..........
fluxn=new((/4,150,192,288/),double) ;frst
fluxn(0,:,:,:) = month_to_annual(myfile5->TOTVEGC_nonlocal(120:1919,:,:),1)
fluxn&year=ti4
fluxn(1,:,:,:) =month_to_annual(myfile6->TOTSOMC_nonlocal(120:1919,:,:),1)
;print(flux(0,0,159,:,:))
;................irri...................
pool3n=new((/4,150,192,288/),double) ;irri
pool3n(0,:,:,:)= month_to_annual(myfile13->TOTVEGC_nonlocal(120:1919,:,:),1)
pool3n&year=ti4
pool3n(1,:,:,:)=month_to_annual(myfile14->TOTSOMC_nonlocal(120:1919,:,:),1)
printVarSummary(pooln(0,1:149,{-60:80},:))

pool=pooln*0.001
flux=fluxn*0.001
pool3=pool3n*0.001
copy_VarMeta(pooln, pool)
copy_VarMeta(fluxn, flux)
copy_VarMeta(pool3n, pool3)
printVarSummary(pool3)
;.............cVeg_nonlocal_crop_expansion...............
Gpool_cveg_diff=new((/149,149,288/),double) 
Gpool_cveg_diff=pool(1,1:149,{-60:80},:)-pool(1,0:148,{-60:80},:)
copy_VarMeta(pool(1,1:149,{-60:80},:), Gpool_cveg_diff)
Gpool_cveg_diff&year=ti4(1:149)
Gpool_cveg_diff_avg=dim_avg_n_Wrap(Gpool_cveg_diff, 0)
;..............climate........................................
tas=new((/149,149,288/),double) ;3varia, 3 simulation
mrsos=new((/149,149,288/),double) ;
tas= month_to_annual(myfile1_4->TREFHT_nonlocal(12::,{-60:80},:),1)
mrsos = month_to_annual(myfile1_2->SOILWATER_10CM_nonlocal(132:1919,{-60:80},:),1)
;..........................................
variables=new((/3,149,149,288/),double) 
variables(0,:,:,:)=1
variables(1,:,:,:)=tas
variables(2,:,:,:)=mrsos
Var_avg=dim_avg_n_Wrap(variables,1)
bstd=new((/3,149,288/),double)  ;standardized b
b=new((/3,149,288/),double)  ;standardized b
pval=new((/3,149,288/),string)   ;The reliability of the t-test result of b.
Fr=new((/149,288/),double)  ;f test result
F_pval=new((/5,192/),string) ; the reliability of the whole regression model
index=new(149,long)
a=new((/3,149/),double)
  do i=0,148
  do j=0,287
;.......................
    index=0
;........................................................
index=ind((.not.ismissing(Gpool_cveg_diff(:,i,j))).and.(.not.ismissing(variables(0,:,i,j)))\
          .and.(.not.ismissing(variables(1,:,i,j))))
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

diro = "/pf/b/b380949/"                     ; Output directory
        filo ="attri_crop_cSoil_CESM.nc"             ; Output file
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

