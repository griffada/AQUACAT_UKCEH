#!/bin/bash

    echo "RCM01 future running"
    bash /prj/aquacat/CodeABG/BASH/890_netcdf_script.sh 01 future > /prj/aquacat/CodeABG/logs/01fut_new.out 2> /prj/aquacat/CodeABG/logs/01fut_new.err

for r in 10 11 12
do
    echo "RCM$r present running"
    bash /prj/aquacat/CodeABG/BASH/890_netcdf_script.sh $r present > /prj/aquacat/CodeABG/logs/"$r"pres_new.out 2> /prj/aquacat/CodeABG/logs/"$r"pres_new.err

    echo "RCM$r future running"
    bash /prj/aquacat/CodeABG/BASH/890_netcdf_script.sh $r future > /prj/aquacat/CodeABG/logs/"$r"fut_new.out 2> /prj/aquacat/CodeABG/logs/"$r"fut_new.err

done
echo "all done"