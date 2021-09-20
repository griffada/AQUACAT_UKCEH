#!/bin/bash

for r in 12 13 15
do
    echo "RCM$r present running"
    bash /prj/aquacat/CodeABG/BASH/456_netcdf_script.sh $r present > /prj/aquacat/CodeABG/logs/"$r"pres_new.out 2> /prj/aquacat/CodeABG/logs/"$r"pres_new.err

done
echo "all done"