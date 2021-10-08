#!/bin/bash

echo "RCM01 present running"
bash /prj/aquacat/CodeABG/BASH/456FF_netcdf_script.sh 15 future FF > /prj/aquacat/CodeABG/logs/z08ff_new.out 2> /prj/aquacat/CodeABG/logs/z08ff_new.err &

echo "RCM01 present running"
bash /prj/aquacat/CodeABG/BASH/456FF_netcdf_script.sh 04 future FF > /prj/aquacat/CodeABG/logs/z10ff_new.out 2> /prj/aquacat/CodeABG/logs/z10ff_new.err &

wait

echo "RCM01 present running"
bash /prj/aquacat/CodeABG/BASH/456FF_netcdf_script.sh 11 future FF > /prj/aquacat/CodeABG/logs/z05ff_new.out 2> /prj/aquacat/CodeABG/logs/z05ff_new.err &

echo "RCM01 present running"
bash /prj/aquacat/CodeABG/BASH/456FF_netcdf_script.sh 12 future FF > /prj/aquacat/CodeABG/logs/z06ff_new.out 2> /prj/aquacat/CodeABG/logs/z06ff_new.err &

echo "RCM01 present running"
bash /prj/aquacat/CodeABG/BASH/456FF_netcdf_script.sh 13 future FF > /prj/aquacat/CodeABG/logs/z07ff_new.out 2> /prj/aquacat/CodeABG/logs/z07ff_new.err &

wait

date -R
echo "Done."