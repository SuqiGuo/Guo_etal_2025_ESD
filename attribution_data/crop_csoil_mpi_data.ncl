
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin
;......................................................................
myfile = addfile("/pf/b/b380949/gridarea/gridarea_sea0.nc","r")
my=addfile("/pf/b/b380949/gridarea/AREA_CESM.nc","r")
myfile1 = addfile("/work/bm1147/b380949/web-monitoring/crop-ctl/Eyr/cVeg/cVeg_crop-ctl_mpiesm_signal-separated_revised.nc", "r")
myfile2 = addfile("/work/bm1147/b380949/web-monitoring/crop-ctl/Eyr/cSoil/cSoil_crop-ctl_mpiesm_signal-separated_revised.nc", "r")
myfile4 = addfile("/work/bm1147/b380949/web-monitoring/crop-ctl/Eyr/cLitter/cLitter_crop-ctl_mpiesm_signal-separated_revised.nc", "r")
myfile3 = addfile("/work/bm1147/b380949/web-monitoring/crop-ctl/Eyr/cProduct/cProduct_crop-ctl_mpiesm_signal-separated_revised.nc", "r")
myfile1_2 = addfile("/work/bm1147/b380949/web-monitoring/crop-ctl/Lmon/mrsos/mrsos_crop-ctl_mpiesm_signal-separated.nc", "r")
;;kg m-2 mrsos;
myfile1_4 = addfile("/work/bm1147/b380949/web-monitoring/crop-ctl/Amon/tas/tas_crop-ctl_mpiesm_signal-separated.nc", "r")
;tas:2-m air tem  K
;myfile1_2 = addfile("/work/bm1147/b380949/web-monitoring/crop-ctl/Amon/pr/pr_crop-ctl_mpiesm_signal-separated.nc", "r")

myfile5 = addfile("/work/bm1147/b380949/web-monitoring/frst-ctl/Eyr/cVeg/cVeg_frst-ctl_mpiesm_signal-separated_revised.nc", "r")
myfile6 = addfile("/work/bm1147/b380949/web-monitoring/frst-ctl/Eyr/cSoil/cSoil_frst-ctl_mpiesm_signal-separated_revised.nc", "r")
myfile8 = addfile("/work/bm1147/b380949/web-monitoring/frst-ctl/Eyr/cLitter/cLitter_frst-ctl_mpiesm_signal-separated_revised.nc", "r")
myfile7 = addfile("/work/bm1147/b380949/web-monitoring/frst-ctl/Eyr/cProduct/cProduct_frst-ctl_mpiesm_signal-separated_revised.nc", "r")

myfile13 = addfile("/work/bm1147/b380949/web-monitoring/irri-crop/Eyr/cVeg/cVeg_irri-crop_mpiesm_signal-separated_revised.nc", "r")
myfile14= addfile("/work/bm1147/b380949/web-monitoring/irri-crop/Eyr/cSoil/cSoil_irri-crop_mpiesm_signal-separated_revised.nc", "r")
myfile15= addfile("/work/bm1147/b380949/web-monitoring/irri-crop/Eyr/cLitter/cLitter_irri-crop_mpiesm_signal-separated_revised.nc", "r")
myfile16= addfile("/work/bm1147/b380949/web-monitoring/irri-crop/Eyr/cProduct/cProduct_irri-crop_mpiesm_signal-separated_revised.nc", "r")

ti=cd_calendar(myfile1->time,-1)
ti2=cd_calendar(myfile5->time,-1)
ti4=cd_calendar(myfile13->time,-1)
;........read in data............................
gridarea = myfile->slm
;.............afforestation..........
flux=new((/4,161,96,192/),double) 
flux(0,:,:,:) = myfile5->cVeg_nonlocal(:,::-1,:)
flux&time=ti2
flux(1,:,:,:) =myfile6->cSoil_nonlocal(:,::-1,:)
flux(2,:,:,:) =myfile8->cLitter_nonlocal(:,::-1,:)
flux(3,:,:,:) =myfile7->cProduct_nonlocal(:,::-1,:)
;print(flux(0,0,159,:,:))
;................irri...................
pool3=new((/4,160,96,192/),double) 
pool3(0,:,:,:)= myfile13->cVeg_nonlocal(:,::-1,:)
pool3&time=ti4
pool3(1,:,:,:)=myfile14->cSoil_nonlocal(:,::-1,:)
pool3(2,:,:,:)= myfile15->cLitter_nonlocal(:,::-1,:)
pool3(3,:,:,:) =myfile16->cProduct_nonlocal(:,::-1,:)
;.........deforestation..... ..................
pool=new((/4,162,96,192/),double) 
pool(0,:,:,:) = myfile1->cVeg_nonlocal(:,::-1,:)
pool&time=ti
pool(1,:,:,:) =myfile2->cSoil_nonlocal(:,::-1,:)
pool(2,:,:,:) =myfile4->cLitter_nonlocal(:,::-1,:)
pool(3,:,:,:) =myfile3->cProduct_nonlocal(:,::-1,:)
;.............cVeg_nonlocal_crop_expansion...............
Gpool_cveg_diff=new((/161,75,192/),double) 
Gpool_cveg_diff=pool(1,1:161,{-60:80},:)-pool(1,0:160,{-60:80},:)
copy_VarMeta(pool(1,1:161,{-60:80},:), Gpool_cveg_diff)
Gpool_cveg_diff&time=ti(1:161)
Gpool_cveg_diff_avg=dim_avg_n_Wrap(Gpool_cveg_diff, 0)
;printVarSummary(Gpool_cveg_diff)
;printVarSummary(Gpool_cveg_diff)
;..............climate........................................
tas=new((/161,75,192/),double) ;3varia, 3 simulation
mrsos=new((/161,75,192/),double) ;
tas= month_to_annual(myfile1_4->tas_nonlocal(12::,{-60:80},:),1)
mrsos = month_to_annual(myfile1_2->mrsos_nonlocal(12::,{-60:80},:),1)
;....when use pr as a independent variable 
;mrsos_new = month_to_annual(myfile1_2->pr_nonlocal(12::,{-60:80},:),1)
;mrsos=mrsos_new*86400 ;/day
;copy_VarMeta(mrsos_new, mrsos)
;........................
;..........................................
variables=new((/3,161,75,192/),double) 
variables(0,:,:,:)=1
variables(1,:,:,:)=tas
variables(2,:,:,:)=mrsos
Var_avg=dim_avg_n_Wrap(variables,1)
bstd=new((/3,75,192/),double)  ;standardized b
b=new((/3,75,192/),double)  ;standardized b
pval=new((/3,75,192/),string)   ;the reliability of t test result of b
Fr=new((/75,192/),double)  ;f test result
F_pval=new((/5,192/),string) ; the reliability of the whole regression
index=new(161,long)
a=new((/3,161/),double)
  do i=0,74
  do j=0,191
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
deltacveg_time=new((/3,161,75,192/),double)  ;standardized b
do i= 0,160
deltacveg_time(:,i,:,:)=b(:,:,:)*variables(:,i,:,:) ;1 t2m,2 mrso
end do
copy_VarMeta(variables,deltacveg_time)
deltacveg_time_new=mask(deltacveg_time, deltacveg_time.eq.0, False)
;;;;;;;;;;;;;;;Ftest..........
prediction=dim_sum_n_Wrap(deltacveg_time_new, 0);3,161,75,192
syy=dim_variance_n_Wrap(Gpool_cveg_diff, 0); Gpool_cveg_diff=new((/161,75,192/),double) 
u_inter=new((/161,75,192/),double)  ;standardized b
do i=0,160
  u_inter(i,:,:)=(prediction(i,:,:)-Gpool_cveg_diff_avg)^2
end do 
copy_VarMeta(prediction, u_inter)
u=dim_sum_n_Wrap(u_inter, 0); 4,75,192
r2=new((/75,192/),double)  ;standardized b
F=new((/75,192/),double)  ;standardized b
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
Q=new((/75,192/),double)  ;standardized b
Q=syy-u

copy_VarMeta(syy, Q)
;...........correlation between predication and Gpool_cveg_diff=r2
correlation=new((/75,192/),double)  ;standardized b
correlation=escorc_n(prediction, Gpool_cveg_diff, 0, 0)
copy_VarMeta(syy, correlation)
;.............data output....................
b@long_name="Regression coefficients"
correlation@long_name="correlation of prediction and original dataset(R-squared)"

diro = "/pf/b/b380949/"                     ; Output directory
        filo ="attri_crop_cSoil.nc"             ; Output file
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

