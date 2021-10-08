#!/bin/bash

echo "RCM01 future running"
bash /prj/aquacat/CodeABG/BASH/106_netcdf_script.sh 01 future > /prj/aquacat/CodeABG/logs/z01f_new.out 2> /prj/aquacat/CodeABG/logs/z01f_new.err &

echo "RCM04 future  running"
bash /prj/aquacat/CodeABG/BASH/106_netcdf_script.sh 04 future > /prj/aquacat/CodeABG/logs/z04f_new.out 2> /prj/aquacat/CodeABG/logs/z04f_new.err &

echo "RCM15 future  running"
bash /prj/aquacat/CodeABG/BASH/106_netcdf_script.sh 15 future > /prj/aquacat/CodeABG/logs/z15f_new.out 2> /prj/aquacat/CodeABG/logs/z15f_new.err &

wait

echo "RCM08 future  running"
bash /prj/aquacat/CodeABG/BASH/106_netcdf_script.sh 08 future > /prj/aquacat/CodeABG/logs/z08f_new.out 2> /prj/aquacat/CodeABG/logs/z08f_new.err &

echo "RCM09 future  running"
bash /prj/aquacat/CodeABG/BASH/106_netcdf_script.sh 09 future > /prj/aquacat/CodeABG/logs/z09f_new.out 2> /prj/aquacat/CodeABG/logs/z09f_new.err &

echo "RCM10 future  running"
bash /prj/aquacat/CodeABG/BASH/106_netcdf_script.sh 10 future > /prj/aquacat/CodeABG/logs/z10f_new.out 2> /prj/aquacat/CodeABG/logs/z10f_new.err &

wait


echo "RCM05 future running"
bash /prj/aquacat/CodeABG/BASH/106_netcdf_script.sh 05 future > /prj/aquacat/CodeABG/logs/z05f_new.out 2> /prj/aquacat/CodeABG/logs/z05f_new.err &

echo "RCM05 future  running"
bash /prj/aquacat/CodeABG/BASH/106_netcdf_script.sh 06 future > /prj/aquacat/CodeABG/logs/z06f_new.out 2> /prj/aquacat/CodeABG/logs/z06f_new.err &

echo "RCM05 future  running"
bash /prj/aquacat/CodeABG/BASH/106_netcdf_script.sh 07 future > /prj/aquacat/CodeABG/logs/z07f_new.out 2> /prj/aquacat/CodeABG/logs/z07f_new.err &

wait

echo "RCM11 future  running"
bash /prj/aquacat/CodeABG/BASH/106_netcdf_script.sh 11 future > /prj/aquacat/CodeABG/logs/z11f_new.out 2> /prj/aquacat/CodeABG/logs/z11f_new.err &

echo "RCM12 future  running"
bash /prj/aquacat/CodeABG/BASH/106_netcdf_script.sh 12 future > /prj/aquacat/CodeABG/logs/z12f_new.out 2> /prj/aquacat/CodeABG/logs/z12f_new.err &

echo "RCM13 future  running"
bash /prj/aquacat/CodeABG/BASH/106_netcdf_script.sh 13 future > /prj/aquacat/CodeABG/logs/z13f_new.out 2> /prj/aquacat/CodeABG/logs/z13f_new.err &

wait

date -R
echo "Done."