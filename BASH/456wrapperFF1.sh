#!/bin/bash

echo "RCM01 future running"
bash /prj/aquacat/CodeABG/BASH/456FF_netcdf_script.sh 05 future FF > /prj/aquacat/CodeABG/logs/z05ff_new.out 2> /prj/aquacat/CodeABG/logs/z05ff_new.err &

echo "RCM01 future running"
bash /prj/aquacat/CodeABG/BASH/456FF_netcdf_script.sh 06 future FF > /prj/aquacat/CodeABG/logs/z06ff_new.out 2> /prj/aquacat/CodeABG/logs/z06ff_new.err &

echo "RCM01 future running"
bash /prj/aquacat/CodeABG/BASH/456FF_netcdf_script.sh 07 future FF > /prj/aquacat/CodeABG/logs/z07ff_new.out 2> /prj/aquacat/CodeABG/logs/z07ff_new.err &

wait

echo "RCM01 future running"
bash /prj/aquacat/CodeABG/BASH/456FF_netcdf_script.sh 08 future FF > /prj/aquacat/CodeABG/logs/z08ff_new.out 2> /prj/aquacat/CodeABG/logs/z08ff_new.err &

echo "RCM01 future running"
bash /prj/aquacat/CodeABG/BASH/456FF_netcdf_script.sh 09 future FF > /prj/aquacat/CodeABG/logs/z09ff_new.out 2> /prj/aquacat/CodeABG/logs/z09ff_new.err &

echo "RCM01 future running"
bash /prj/aquacat/CodeABG/BASH/456FF_netcdf_script.sh 10 future FF > /prj/aquacat/CodeABG/logs/z10ff_new.out 2> /prj/aquacat/CodeABG/logs/z10ff_new.err &

wait

date -R
echo "Done."