# Script to download wind data from CCMP

#set -x
export LANG=en_US.UTF-8
export HDF5_DISABLE_VERSION_CHECK=1

############# Parametros de entrada ##########
inicio=$1  #{YYYYMMDD}
limite=$2  #{YYYYMMDD}  
current=${inicio}

#PREFIX="qwqg00" # qwqg00_2018012600_10mwind.gz and qwqg00.2018013100.fields.gz
PATH_OUT=/home/santanarc/ww3/data/UKMO/

while [ ${current} -le ${limite} ]; do

  YEAR=`date -u -d "${current}" +%Y`
  MONTH=`date -u -d "${current}" +%m`
  DAY=`date -u -d "${current}" +%d`

  PATH_OUTPUT=$PATH_OUT/${YEAR}/${MONTH}/${DAY}/00/
  echo ${PATH_OUTPUT}
  mkdir -p ${PATH_OUTPUT}

  echo "Copying day: ${current}"

  SOURCE=/niwa/fitzroy_archive/ecoconnect_oper/EFS/UKMO/${YEAR}/${MONTH}/${DAY}/00/qwqg00_${current}00_10mwind.gz # ${current}00*.gz #"{_10mwind,.fields}.gz"
  #SOURCE=/niwa/fitzroy_archive/ecoconnect_oper/EFS/UKMO/${YEAR}/${MONTH}/${DAY}/00/q*10mwind*.gz # ${current}00*.gz #"{_10mwind,.fields}.gz"
  echo "${SOURCE}"
  rsync -avz --progress santanarc@login.kupe.niwa.co.nz:${SOURCE} ${PATH_OUTPUT}.

  SOURCE=/niwa/fitzroy_archive/ecoconnect_oper/EFS/UKMO/${YEAR}/${MONTH}/${DAY}/00/qwgl_daily_${current}00_ice.gz # ${current}00*.gz #"{_10mwind,.fields}.gz"
  #SOURCE=/niwa/fitzroy_archive/ecoconnect_oper/EFS/UKMO/${YEAR}/${MONTH}/${DAY}/00/q*ice.gz # ${current}00*.gz #"{_10mwind,.fields}.gz"
  echo "${SOURCE}"
  rsync -avz --progress santanarc@login.kupe.niwa.co.nz:${SOURCE} ${PATH_OUTPUT}.

  #exit 0

  #mv *{wind,ice}*.gz ${PATH_OUTPUT}.

  wait

  current=`date -u -d "${current} +1 day" +%Y%m%d`

done


