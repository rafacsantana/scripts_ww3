# Script to download wind data from CCMP

#set -x
export LANG=en_US.UTF-8
export HDF5_DISABLE_VERSION_CHECK=1

############# Parametros de entrada ##########
inicio=$1  #{YYYYMMDD}
limite=$2  #{YYYYMMDD}  
current=${inicio}

#expt=NZWAVE-HR
expt=NZCSM
#expt=GLOBALWAVE

PREFIX="nwpsfc_"
SUFIX="00-utc_nzcsm_" # 000 to 048
PATH_OUT=/niwa/archive/ecoconnect/EFS/$expt 

#/niwa/archive/ecoconnect/EFS/NZCSM/2021/06/24/00/nwpsfc_2021062400-utc_nzcsm_010.um

rm -rf ~/tape_files.txt

while [ ${current} -le ${limite} ]; do

  YEAR=`date -u -d "${current}" +%Y`
  MONTH=`date -u -d "${current}" +%m`
  DAY=`date -u -d "${current}" +%d`

  #PATH_OUTPUT=$PATH_OUT/${YEAR}/${MONTH}/${DAY}/00/
  #mkdir -p $PATH_OUTPUT # /${YEAR}/${MONTH}/${DAY}/00/

  echo "Appending day: ${current}"

  for i in {0..48}; do

    ii=$(printf "%03d" ${i})

    SOURCE=$PATH_OUT/${YEAR}/${MONTH}/${DAY}/00/${PREFIX}${current}${SUFIX}${ii}.um 

    #echo "${SOURCE}"

    #exit 0

    echo -e ${SOURCE} >> ~/tape_files.txt

    wait

  done

  current=`date -u -d "${current} +1 day" +%Y%m%d`

done


