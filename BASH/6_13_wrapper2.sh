#!/bin/bash

echo "RCM05 future running"
bash /prj/aquacat/CodeABG/BASH/456FF_netcdf_script.sh 05 future > /prj/aquacat/CodeABG/logs/z05f_new.out 2> /prj/aquacat/CodeABG/logs/z05f_new.err &

echo "RCM06 future running"
bash /prj/aquacat/CodeABG/BASH/456FF_netcdf_script.sh 06 future > /prj/aquacat/CodeABG/logs/z06f_new.out 2> /prj/aquacat/CodeABG/logs/z06f_new.err &

echo "RCM07 future running"
bash /prj/aquacat/CodeABG/BASH/456FF_netcdf_script.sh 07 future > /prj/aquacat/CodeABG/logs/z07f_new.out 2> /prj/aquacat/CodeABG/logs/z07f_new.err &

wait

echo "RCM08 future running"
bash /prj/aquacat/CodeABG/BASH/456FF_netcdf_script.sh 08 future > /prj/aquacat/CodeABG/logs/z08f_new.out 2> /prj/aquacat/CodeABG/logs/z08f_new.err &

echo "RCM09 future running"
bash /prj/aquacat/CodeABG/BASH/456FF_netcdf_script.sh 09 future > /prj/aquacat/CodeABG/logs/z09f_new.out 2> /prj/aquacat/CodeABG/logs/z09f_new.err &

echo "RCM10 future running"
bash /prj/aquacat/CodeABG/BASH/456FF_netcdf_script.sh 10 future > /prj/aquacat/CodeABG/logs/z10f_new.out 2> /prj/aquacat/CodeABG/logs/z10f_new.err &

wait

date -R
echo "Done."