#!/bin/bash

# THIS SCRIPT DOES THE PREPROCESSING AND EXECUTES THE PYTHON SCRIPT TO
# SEPARATE LOCAL, NON-LOCAL AND TOTAL IMPACTS OF LCLM CHANGE FOR
# LAMACLIMA WP1 SIMULATIONS: CTL, CROP, IRR, FRST, HARV

# BEFORE EXECUTING THE SCRIPT DEFINE THE FOLLOWING VARIABLES BELOW :
# 1. CMOR TABLE ID (--> TEMPORAL RESOLUTION)
# 2. CMORvar
# 3. SCENARIO COMBINATION

# EXECUTE THIS SCRIPT VIA: sbatch signal_seperation.sh

#SBATCH --job-name=signal_separation      # Specify job name
#SBATCH --partition=shared                # Specify partition name
#SBATCH --ntasks=2                        # Specify max. number of tasks to be invoked
#SBATCH --cpus-per-task=16                # Specify number of CPUs per task
#SBATCH --time=02:00:00                   # Set a limit on the total run time
#SBATCH --account=bm1147                  # Charge resources on this project account
#SBATCH --output=signal_separation.out    # File name for standard output
#SBATCH --error=signal_separation.out     # File name for standard error output

# Bind your OpenMP threads --> I have to admit I have no idea if that's needed or what is happening here....
export OMP_NUM_THREADS=8
export KMP_AFFINITY=verbose,granularity=core,compact,1
export KMP_STACKSIZE=64m

module unload python
module load python/3.5.2

MODEL="mpiesm" #, "cesm" "ecearth" "mpiesem"
CMOR_TABLE="Amon" # "Lmon" "LImon" "Emon" "Lmon" "Eyr" "Lmon"; choose "cesm" for all cesm variables as they are currently stored in /scratch/b/b380948/signal_separation/cesm_output/
# CMOR_VAR_LIST="cSoil cLitter cVeg cProduct"
# CMOR_VAR_LIST="gpp npp nbp ra rh"
# CMOR_VAR_LIST="mrsol nep"
# CMOR_VAR_LIST="TSA TG EFLX_LH_TOT FSH PRECC PRECL WIND Q2M" # cesm model output data used with cesm table
CMOR_VAR_LIST="pr" # cesm model output data used with cesm table
# CMOR_VAR_LIST="PRECC PRECL" # still missing for frst-ctl and crop-ctl 
# CMOR_VAR_LIST="tas"
SIM1="frst" # "frst" "crop" "irr" "harv"
SIM2="ctl" # "ctl" "crop" "frst"
SCENARIO_COMBNATION="${SIM1}-${SIM2}" #, "irr-crop" "irr-ctl" "frst-ctl" "harv-frst" "harv-ctl"


for CMOR_VAR in ${CMOR_VAR_LIST}; do

if  [ "${MODEL}" == "mpiesm" ]; then
  # DIR1=/work/bm1147/b380949/mpiesm-1.2.01-release/experiments/ssp245_r1i1p1f1-LR_suq09/archive/CMIP6/ScenarioMIP/MPI-M/MPI-ESM1-2-LR/ssp245/r1i1p1f1
    DIR1=/work/bm1147/b380949/mpiesm-1.2.01-release/experiments/FRST/archive/CMIP6/ScenarioMIP/MPI-M/MPI-ESM1-2-LR/ssp245/r1i1p1f1

    DIR2=/work/bm1147/b380949/mpiesm-1.2.01-release/experiments/CTL/archive/CMIP6/ScenarioMIP/MPI-M/MPI-ESM1-2-LR/ssp245/r1i1p1f1


elif [ "$SCENARIO_COMBNATION" == "crop-ctl" ] && [ "${MODEL}" == "cesm" ]; then
    DIR1=/scratch/b/b380948/signal_separation/cesm_output/crop
    DIR2=/scratch/b/b380948/signal_separation/cesm_output/ctl
    

elif [ "$SCENARIO_COMBNATION" == "frst-ctl" ] && [ "${MODEL}" == "cesm" ]; then
    DIR1=/scratch/b/b380948/signal_separation/cesm_output/frst
    DIR2=/scratch/b/b380948/signal_separation/cesm_output/ctl

# 
# elif [ "$SCENARIO" == "irr" ]; then
#     DIR=
# 
# elif [ "$SCENARIO" == "frst" ]; then
#     DIR=
#     
# elif [ "$SCENARIO" == "harv" ]; then
#     DIR=    
fi

# SCRATCH=/scratch/b/b380948/prepared_for_arch/${SCENARIO}/
SCRATCH="/work/bm1147/b380949/web-monitoring"

# CREATE TEMPORARY FOLDERS ON SCRATCH TO SAFE THE DATA
if [ ! -d "${SCRATCH}" ]; then
  echo "Folder ${SCRATCH} does not exist and will be made."
  mkdir $SCRATCH
fi
cd $SCRATCH

# if [ ${MODEL} == "cesm" ] && [ ! -d "cesm_output" ]; then  
#   echo "Folder ${SCRATCH}/cesm_output does not exist and will be made."
#   mkdir cesm_output  
# else
#   echo "Folder ${SCRATCH}/cesm_output already exists."
# fi  
# cd cesm_output

if [ ! -d "${SCENARIO_COMBNATION}" ]; then
  echo "Folder ${SCRATCH}/${SCENARIO_COMBNATION} does not exist and will be made."
  mkdir ${SCENARIO_COMBNATION}
else
  echo "Folder ${SCRATCH}/${SCENARIO_COMBNATION} already exists."
fi
cd ${SCENARIO_COMBNATION}

if [ ! -d "${CMOR_TABLE}" ]; then
  echo "Folder ${CMOR_TABLE} does not exist and will be made."
  mkdir ${CMOR_TABLE}
else 
  echo "Folder ${CMOR_TABLE} already exists."  
fi
cd ${CMOR_TABLE}

if [ ! -d "${CMOR_VAR}" ]; then
  echo "Folder ${CMOR_VAR} does not exist and will be made."
  mkdir ${CMOR_VAR}
else
  echo "Folder ${CMOR_VAR} already exists."
fi
cd ${CMOR_VAR}


# COPY ALL MODEL OUTPUT FILES, MERGE AND CREATE DIFFERENCE FILES
if [ "${MODEL}" == "mpiesm" ]; then
  FILES1=`find ${DIR1}/${CMOR_TABLE}/${CMOR_VAR}/* -name "${CMOR_VAR}*"`
  FILES2=`find ${DIR2}/${CMOR_TABLE}/${CMOR_VAR}/* -name "${CMOR_VAR}*"`  

  cdo mergetime ${FILES1} ${CMOR_VAR}_${SIM1}.nc
  cdo mergetime ${FILES2} ${CMOR_VAR}_${SIM2}.nc 
#ncks  -O -d time,10,160 ${CMOR_VAR}_${SIM2}.nc ${CMOR_VAR}_${SIM2}.nc
elif [ "${MODEL}" == "cesm" ] && [ "${SIM1}" == "crop"  ]; then
  cp ${DIR1}/${CMOR_VAR}_CROP_LAMACLIMA.e211.B2000cmip6.f09_g17.crop-i308.CROP.clm2.h0.nc ${CMOR_VAR}_${SIM1}.nc
  cp ${DIR2}/${CMOR_VAR}_CTL_LAMACLIMA.e211.B2000cmip6.f09_g17.control-i196.clm2.h0.nc ${CMOR_VAR}_${SIM2}.nc

elif [ "${MODEL}" == "cesm" ] && [ "${SIM1}" == "frst"  ]; then
  cp ${DIR1}/${CMOR_VAR}_FRST_LAMACLIMA.e211.B2000cmip6.f09_g17.forest-i308.clm2.h0.nc ${CMOR_VAR}_${SIM1}.nc
  cp ${DIR2}/${CMOR_VAR}_CTL_LAMACLIMA.e211.B2000cmip6.f09_g17.control-i196.clm2.h0.nc ${CMOR_VAR}_${SIM2}.nc
  
fi    

if [ "$SCENARIO_COMBNATION" == "frst-ctl" ] && [ "${CMOR_TABLE}" == "Lmon" -o "${CMOR_TABLE}" == "Amon" ]; then
    ncks -O -d time,0,1931 ${CMOR_VAR}_${SIM2}.nc ${CMOR_VAR}_${SIM2}.nc  #Lmon
elif [ "$SCENARIO_COMBNATION" == "frst-ctl" ] && [ "${CMOR_TABLE}" == "Eyr" ]; then
    ncks -O -d time,0,160 ${CMOR_VAR}_${SIM2}.nc ${CMOR_VAR}_${SIM2}.nc   #Eyr
elif [ "$SCENARIO_COMBNATION" == "harv-ctl" ] && [ "${CMOR_TABLE}" == "Lmon" -o "${CMOR_TABLE}" == "Amon" ]; then
    ncks  -O -d time,480,1919 ${CMOR_VAR}_${SIM2}.nc ${CMOR_VAR}_${SIM2}.nc
elif [ "$SCENARIO_COMBNATION" == "harv-ctl" ] && [ "${CMOR_TABLE}" == "Eyr" ]; then
    ncks  -O -d time,40,159 ${CMOR_VAR}_${SIM2}.nc ${CMOR_VAR}_${SIM2}.nc
elif [ "$SCENARIO_COMBNATION" == "harv-frst" ] && [ "${CMOR_TABLE}" == "Lmon" -o "${CMOR_TABLE}" == "Amon" ]; then
    ncks  -O -d time,360,1799 ${CMOR_VAR}_${SIM2}.nc ${CMOR_VAR}_${SIM2}.nc
elif [ "$SCENARIO_COMBNATION" == "harv-frst" ] && [ "${CMOR_TABLE}" == "Eyr" ]; then
    ncks  -O -d time,30,149 ${CMOR_VAR}_${SIM2}.nc ${CMOR_VAR}_${SIM2}.nc
elif [ "$SIM1" == "irri" ] && [ "${CMOR_TABLE}" == "Lmon" -o "${CMOR_TABLE}" == "Amon" ]; then
    ncks -O -d time,0,1919 ${CMOR_VAR}_${SIM2}.nc ${CMOR_VAR}_${SIM2}.nc  #Lmon
elif [ "$SIM1" == "irri" ] && [ "${CMOR_TABLE}" == "Eyr" ]; then
    ncks -O -d time,0,159 ${CMOR_VAR}_${SIM2}.nc ${CMOR_VAR}_${SIM2}.nc   #Eyr
  fi


cdo sub ${CMOR_VAR}_${SIM1}.nc ${CMOR_VAR}_${SIM2}.nc ${CMOR_VAR}_${SIM1}-${SIM2}_${MODEL}.nc

# JUST OVERWRITE CALCULATED FILES THAT EXIST, OTHERWISE ONE WOULD NOT CALL THIS SCRIPT....
# if [ ! -f ${CMOR_VAR}_${SIM1}-${SIM2}_${MODEL}_signal-separated.nc ]; then
#   cp ${CMOR_VAR}_${SIM1}-${SIM2}_${MODEL}.nc ${CMOR_VAR}_${SIM1}-${SIM2}_${MODEL}_signal-separated.nc
# else
#   echo "File ${CMOR_VAR}_${SIM1}-${SIM2}_${MODEL}_signal-separated.nc already exists...."
# fi
# if [ "${MODEL}" ==  ]
cp ${CMOR_VAR}_${SIM1}-${SIM2}_${MODEL}.nc ${CMOR_VAR}_${SIM1}-${SIM2}_${MODEL}_signal-separated.nc
#NC_FILE=${SCRATCH}/${SCENARIO_COMBNATION}/${CMOR_TABLE}/${CMOR_VAR}/${CMOR_VAR}_${SIM1}-${SIM2}_${MODEL}.nc
NC_FILE=${SCRATCH}/${SCENARIO_COMBNATION}/${CMOR_TABLE}/${CMOR_VAR}/${CMOR_VAR}_${SIM1}-${SIM2}_${MODEL}_signal-separated.nc

echo "Calculate signal separation for variable ${CMOR_VAR} in file ${NC_FILE}"
python -i /pf/b/b380949/separation/orig_sep/prepare_local_nonlocal_total_effects2.py ${NC_FILE} ${CMOR_VAR}

done

