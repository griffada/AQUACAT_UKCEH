#!/bin/bash

echo "RCM01 present running"
bash /prj/aquacat/CodeABG/BASH/456FF_netcdf_script.sh 01 future > /prj/aquacat/CodeABG/logs/z01f_new.out 2> /prj/aquacat/CodeABG/logs/z01f_new.err &

echo "RCM15 present running"
bash /prj/aquacat/CodeABG/BASH/456FF_netcdf_script.sh 15 future > /prj/aquacat/CodeABG/logs/z15f_new.out 2> /prj/aquacat/CodeABG/logs/z15f_new.err &

echo "RCM04 present running"
bash /prj/aquacat/CodeABG/BASH/456FF_netcdf_script.sh 04 future > /prj/aquacat/CodeABG/logs/z04f_new.out 2> /prj/aquacat/CodeABG/logs/z04f_new.err &

wait

echo "RCM11 present running"
bash /prj/aquacat/CodeABG/BASH/456FF_netcdf_script.sh 11 future > /prj/aquacat/CodeABG/logs/z11f_new.out 2> /prj/aquacat/CodeABG/logs/z11f_new.err &

echo "RCM12 present running"
bash /prj/aquacat/CodeABG/BASH/456FF_netcdf_script.sh 12 future > /prj/aquacat/CodeABG/logs/z12f_new.out 2> /prj/aquacat/CodeABG/logs/z12f_new.err &

echo "RCM13 present running"
bash /prj/aquacat/CodeABG/BASH/456FF_netcdf_script.sh 13 future > /prj/aquacat/CodeABG/logs/z13f_new.out 2> /prj/aquacat/CodeABG/logs/z13f_new.err &

wait

date -R
echo "Done."