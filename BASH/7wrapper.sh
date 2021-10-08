#!/bin/bash

for r in 06 07 08
    do
    for p in present future
    do

        echo "RCM$r $p running"
        bash /prj/aquacat/CodeABG/BASH/456g_netcdf_script.sh $r $p > /prj/aquacat/CodeABG/logs/"$r"_"$p"_new.out 2> /prj/aquacat/CodeABG/logs/"$r"_"$p"_new.err

    done
echo "all done"