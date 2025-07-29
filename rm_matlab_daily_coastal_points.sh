# Script to download wind data from CCMP

#set -x
export LANG=en_US.UTF-8
export HDF5_DISABLE_VERSION_CHECK=1

############# Parametros de entrada ##########
expt=$1 # NZWAVE-GFDL-CCAM
YEAR=$2  #{YYYY}
#current=${inicio}

#expt=NZCSM
#expt=GLOBALWAVE

PREFIX="i"
SUFIX="grid_point__" # 000 to 048
PATH_OUT=~/hindcast/$expt/matlab/ 

#/niwa/archive/ecoconnect/EFS/NZCSM/2021/06/24/00/nwpsfc_2021062400-utc_nzcsm_010.um

#rm -rf ~/tape_files.txt

#for i 1..12
for i in {1..12}; do
#while [ ${current} -le ${limite} ]; do

  month=$(printf "%02d" ${i})

  echo "Removing month: ${YEAR}${month}"

  #for i in {0..48}; do


  SOURCE=$PATH_OUT${PREFIX}*${SUFIX}${YEAR}${month}*.mat 

  echo "rm ${SOURCE}"
  rm $SOURCE

    #exit 0

#    echo -e ${SOURCE} >> ~/tape_files.txt

    wait

#  done

#  current=`date -u -d "${current} +1 day" +%Y%m%d`

done


