#!/bin/bash

for r in 04 05 06 07 08 09 10 11 12 13 15
do
    echo "RCM$r present running"
    bash /prj/aquacat/CodeABG/BASH/456_netcdf_script.sh $r future > /prj/aquacat/CodeABG/logs/"$r"fut_new.out 2> /prj/aquacat/CodeABG/logs/"$r"fut_new.err

done
echo "all done"