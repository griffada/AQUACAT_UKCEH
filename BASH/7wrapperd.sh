#!/bin/bash

    echo "RCM01 future running"
    bash /prj/aquacat/CodeABG/BASH/7_13_netcdf_script.sh 01 future > /prj/aquacat/CodeABG/logs/01fut_new.out 2> /prj/aquacat/CodeABG/logs/01fut_new.err

for r in 13 15
do
    echo "RCM$r present running"
    bash /prj/aquacat/CodeABG/BASH/7_13_netcdf_script.sh $r present > /prj/aquacat/CodeABG/logs/"$r"pres_new.out 2> /prj/aquacat/CodeABG/logs/"$r"pres_new.err

    echo "RCM$r future running"
    bash /prj/aquacat/CodeABG/BASH/7_13_netcdf_script.sh $r future > /prj/aquacat/CodeABG/logs/"$r"fut_new.out 2> /prj/aquacat/CodeABG/logs/"$r"fut_new.err

done
echo "all done"