# Script to download wind data from CCMP

#set -x
export LANG=en_US.UTF-8
export HDF5_DISABLE_VERSION_CHECK=1

############# Parametros de entrada ##########
inicio=$1  #{YYYYMMDD}
limite=$2  #{YYYYMMDD}  
current=${inicio}

PREFIX="qwqg00" # qwqg00_2018012600_10mwind.gz and qwqg00.2018013100.fields.gz
PATH_OUT=/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/hindcast/UKMO/

while [ ${current} -le ${limite} ]; do

  YEAR=`date -u -d "${current}" +%Y`
  MONTH=`date -u -d "${current}" +%m`
  DAY=`date -u -d "${current}" +%d`

  PATH_OUTPUT=$PATH_OUT/${YEAR}/${MONTH}/${DAY}/00/
  mkdir -p $PATH_OUTPUT # /${YEAR}/${MONTH}/${DAY}/00/

  echo "Copying day: ${current}"

  SOURCE=/scale_wlg_devoper/filesets/archive/ecoconnect/EFS/UKMO/${YEAR}/${MONTH}/${DAY}/00/${PREFIX}*.gz # ${current}00*.gz #"{_10mwind,.fields}.gz"

  echo "${SOURCE}"

  #exit 0

  cp -r ${SOURCE} ${PATH_OUTPUT}.

  wait

  current=`date -u -d "${current} +1 day" +%Y%m%d`

done


