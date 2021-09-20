#!/bin/bash

for r in 04 05 06
do
    echo "RCM$r present running"
    bash /prj/aquacat/CodeABG/BASH/7_13_netcdf_script.sh $r present > /prj/aquacat/CodeABG/logs/"$r"pres_new.out 2> /prj/aquacat/CodeABG/logs/"$r"pres_new.err

    echo "RCM$r future running"
    bash /prj/aquacat/CodeABG/BASH/7_13_netcdf_script.sh $r future > /prj/aquacat/CodeABG/logs/"$r"fut_new.out 2> /prj/aquacat/CodeABG/logs/"$r"fut_new.err

done
echo "all done"