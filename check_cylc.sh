#!/bin/bash

#!/bin/bash
# script to push UKMO data to Maui

stime=$1

# infite loop  
a=1
while [ $a -lt 200000 ]; do

  # check cylc code
  sc=`cylc scan`
  k=1
  kk=0
  for s in $sc; do
    k=$((k+1))
  
    if [ `expr $k % 2` == 0 ]; then
      kk=$((kk+1))
      echo $kk
      echo $s
      #cylc dump $s | grep state
      cylc dump $s | grep wave_run,
      cylc dump $s | grep failed
    fi
  done

  echo ''
  echo 'waiting for '$stime' seconds ...'
  echo ''
  echo 'If the model crashes remember to RERUN the incomplete days below:'
  echo '1997060700'
  echo ''
  date
  sleep $stime

done

