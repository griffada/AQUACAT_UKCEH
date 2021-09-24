#!/bin/bash

for p in present future
do
    

    echo "RCM$r present running"
    bash /prj/aquacat/CodeABG/BASH/7_13_netcdf_script.sh 04 $p > /prj/aquacat/CodeABG/logs/"$r"pres_new.out 2> /prj/aquacat/CodeABG/logs/"$r"pres_new.err

    echo "RCM$r future running"
    bash /prj/aquacat/CodeABG/BASH/7_13_netcdf_script.sh 04 $p > /prj/aquacat/CodeABG/logs/"$r"fut_new.out 2> /prj/aquacat/CodeABG/logs/"$r"fut_new.err

done
echo "all done"