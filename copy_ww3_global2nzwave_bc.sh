# Script to convert WW3 nc to nc4
#Usage
#set-roms-intel # loading netcdf libraries
#./convert_ww3_nc_to_nc4.sh 20000101 20001231 NZWAVE-ERA5

#set -x
export LANG=en_US.UTF-8
export HDF5_DISABLE_VERSION_CHECK=1

############# Parametros de entrada ##########
inicio=$1  #{YYYYMMDD}
limite=$2  #{YYYYMMDD}  
current=${inicio}

#expt=$3 # NZWAVE-HR

PATH_IN=/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/hindcast/GLOBALWAVE-GFDL-CCAM
PATH_OUT=/scale_wlg_persistent/filesets/project/niwa03150/santanarc/hindcast/NZWAVE-GFDL-CCAM
PREFIX="nest_globalwave-17_to_nzwave-4" # qwqg00_2018012600_10mwind.gz and qwqg00.2018013100.fields.gz

#cp ~/hindcast/GLOBALWAVE-GFDL-CCAM/2030/*/*/00/nest_globalwave-17_to_nzwave-4*ww3 ~/waves/hindcast/NZWAVE-GFDL-CCAM/boundary_conditions/

while [ ${current} -le ${limite} ]; do

  YEAR=`date -u -d "${current}" +%Y`
  MONTH=`date -u -d "${current}" +%m`
  DAY=`date -u -d "${current}" +%d`

  #PATH_OUTPUT=$PATH_OUT/${YEAR}/${MONTH}/${DAY}/00/
  #mkdir -p $PATH_OUTPUT # /${YEAR}/${MONTH}/${DAY}/00/

  echo "Copying nesting file day: ${current}"

  SOURCE=$PATH_IN/${YEAR}/${MONTH}/${DAY}/00/${PREFIX}*ww3 # ${current}00*.gz #"{_10mwind,.fields}.gz"
  S_OUT=$PATH_OUT/${YEAR}/${MONTH}/${DAY}/00/ # ${current}00*.gz #"{_10mwind,.fields}.gz"
  mkdir -p $S_OUT 
  
  echo "${SOURCE}"

  cp $SOURCE $S_OUT/
  

  #FILES=`ls $SOURCE`
  #for f in ${FILES[@]}; do
  #
  # echo "Processing file: ${f}"
  #
  # base=`basename ${f} .nc`
  # #base=small_$base
  #
  # #nccopy -d1 classic.nc compressed.nc
  # nccopy -d9 -s ${f} ~/${base}.nc4
  # mv  ~/${base}.nc4 ${f}
  #
  #done

  wait

  current=`date -u -d "${current} +1 day" +%Y%m%d`

done


