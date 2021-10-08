#!/bin/bash

echo "RCM05 present running"
bash /prj/aquacat/CodeABG/BASH/456FF_netcdf_script.sh 05 present > /prj/aquacat/CodeABG/logs/z05p_new.out 2> /prj/aquacat/CodeABG/logs/z05p_new.err &

echo "RCM06 present running"
bash /prj/aquacat/CodeABG/BASH/456FF_netcdf_script.sh 06 present > /prj/aquacat/CodeABG/logs/z06p_new.out 2> /prj/aquacat/CodeABG/logs/z06p_new.err &

echo "RCM07 present running"
bash /prj/aquacat/CodeABG/BASH/456FF_netcdf_script.sh 07 present > /prj/aquacat/CodeABG/logs/z07p_new.out 2> /prj/aquacat/CodeABG/logs/z07p_new.err &

wait

echo "RCM08 present running"
bash /prj/aquacat/CodeABG/BASH/456FF_netcdf_script.sh 08 present > /prj/aquacat/CodeABG/logs/z08p_new.out 2> /prj/aquacat/CodeABG/logs/z08p_new.err &

echo "RCM09 present running"
bash /prj/aquacat/CodeABG/BASH/456FF_netcdf_script.sh 09 present > /prj/aquacat/CodeABG/logs/z09p_new.out 2> /prj/aquacat/CodeABG/logs/z09p_new.err &

echo "RCM10 present running"
bash /prj/aquacat/CodeABG/BASH/456FF_netcdf_script.sh 10 present > /prj/aquacat/CodeABG/logs/z10p_new.out 2> /prj/aquacat/CodeABG/logs/z10p_new.err &

wait

date -R
echo "Done."