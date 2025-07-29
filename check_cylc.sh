#!/bin/bash

#!/bin/bash

stime=$1

# infite loop  
a=1
while [ $a -lt 200000 ]; do

  # check cylc code
  sc=`cylc scan`
  k=1
  kk=0
  date
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
  date
  #echo 'If the model crashes remember to RERUN the incomplete days below:'
  #echo 'cylc run nzwave_gfd 20410221T00 -- to run from a certain datetime'
  #echo 'cylc stop --kill nzwave_gfd -- to stop a nzwave_gfd suite'
  echo ''
  sleep $stime

done

