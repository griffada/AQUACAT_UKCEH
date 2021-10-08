#!/bin/bash

echo "RCM01 present running"
bash /prj/aquacat/CodeABG/BASH/456FF_netcdf_script.sh 01 present > /prj/aquacat/CodeABG/logs/z01p_new.out 2> /prj/aquacat/CodeABG/logs/z01p_new.err &

echo "RCM15 present running"
bash /prj/aquacat/CodeABG/BASH/456FF_netcdf_script.sh 15 present > /prj/aquacat/CodeABG/logs/z15p_new.out 2> /prj/aquacat/CodeABG/logs/z15p_new.err &

echo "RCM04 present running"
bash /prj/aquacat/CodeABG/BASH/456FF_netcdf_script.sh 04 present > /prj/aquacat/CodeABG/logs/z04p_new.out 2> /prj/aquacat/CodeABG/logs/z04p_new.err &

wait

echo "RCM11 present running"
bash /prj/aquacat/CodeABG/BASH/456FF_netcdf_script.sh 11 present > /prj/aquacat/CodeABG/logs/z11p_new.out 2> /prj/aquacat/CodeABG/logs/z11p_new.err &

echo "RCM12 present running"
bash /prj/aquacat/CodeABG/BASH/456FF_netcdf_script.sh 12 present > /prj/aquacat/CodeABG/logs/z12p_new.out 2> /prj/aquacat/CodeABG/logs/z12p_new.err &

echo "RCM13 present running"
bash /prj/aquacat/CodeABG/BASH/456FF_netcdf_script.sh 13 present > /prj/aquacat/CodeABG/logs/z13p_new.out 2> /prj/aquacat/CodeABG/logs/z13p_new.err &

wait

date -R
echo "Done."