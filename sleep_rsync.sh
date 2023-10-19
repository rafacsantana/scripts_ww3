#!/bin/bash
# script to push UKMO data to Maui

a=1
while [ $a -lt 200000 ]; do 

  rsync -avz --progress /home/santanarc/ww3/data/UKMO/* santanarc@maui:~/hindcast/UKMO/
  sleep 3600

done



