# Script to download wind data from CCMP

#set -x
export LANG=en_US.UTF-8
export HDF5_DISABLE_VERSION_CHECK=1

############# Parametros de entrada ##########
inicio=$1  #{YYYYMMDD}
limite=$2  #{YYYYMMDD}  
current=${inicio}

expt=NZWAVE-HR

PREFIX="ww3g" # qwqg00_2018012600_10mwind.gz and qwqg00.2018013100.fields.gz

PATH_OUT=/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/hindcast/$expt/

while [ ${current} -le ${limite} ]; do

  YEAR=`date -u -d "${current}" +%Y`
  MONTH=`date -u -d "${current}" +%m`
  DAY=`date -u -d "${current}" +%d`

  PATH_OUTPUT=$PATH_OUT/${YEAR}/${MONTH}/${DAY}/00/
  mkdir -p $PATH_OUTPUT # /${YEAR}/${MONTH}/${DAY}/00/

  echo "Copying day: ${current}"

  SOURCE=/scale_wlg_devoper/filesets/archive/ecoconnect/EFS/$expt/${YEAR}/${MONTH}/${DAY}/00/${PREFIX}*.nc # ${current}00*.gz #"{_10mwind,.fields}.gz"

  echo "${SOURCE}"

  #exit 0

  cp -r ${SOURCE} ${PATH_OUTPUT}.

  wait

  current=`date -u -d "${current} +1 day" +%Y%m%d`

done


