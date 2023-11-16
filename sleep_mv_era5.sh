#!/bin/bash
# script to push UKMO data to Maui

prefix=~/uoa03669/rsan613/get_ERA5/ERA5

a=1
while [ $a -lt 200000 ]; do

  for v in ${prefix}*.nc; do

    ab_size=$(du -c $v | tail -n 1 | awk '{print $1}')
    echo $v 'size' $ab_size

    if [ $ab_size -gt 17500000 ]; then
    #while [ $ab_size -lt 17500000 ]; do
            #sleep 30
            ab_size=$(du -c $v | tail -n 1 | awk '{print $1}')
            echo $v 'size' $ab_size
            echo 'mv '$v' ~/hindcast/ERA5/'
            mv $v ~/hindcast/ERA5/
    fi

  done
  echo 'waiting for five minutes ...'
  sleep 300

done



