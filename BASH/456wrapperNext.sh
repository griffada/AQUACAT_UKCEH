#!/bin/bash

echo "RCM01 present running"
bash /prj/aquacat/CodeABG/BASH/456i_netcdf_script.sh 05 present > /prj/aquacat/CodeABG/logs/z05p_new.out 2> /prj/aquacat/CodeABG/logs/z05p_new.err &

echo "RCM01 present running"
bash /prj/aquacat/CodeABG/BASH/456i_netcdf_script.sh 06 present > /prj/aquacat/CodeABG/logs/z06p_new.out 2> /prj/aquacat/CodeABG/logs/z06p_new.err &

echo "RCM01 present running"
bash /prj/aquacat/CodeABG/BASH/456i_netcdf_script.sh 07 present > /prj/aquacat/CodeABG/logs/z07p_new.out 2> /prj/aquacat/CodeABG/logs/z07p_new.err &

wait

echo "RCM01 present running"
bash /prj/aquacat/CodeABG/BASH/456i_netcdf_script.sh 05 future > /prj/aquacat/CodeABG/logs/z05f_new.out 2> /prj/aquacat/CodeABG/logs/z05f_new.err &

echo "RCM01 present running"
bash /prj/aquacat/CodeABG/BASH/456i_netcdf_script.sh 06 future > /prj/aquacat/CodeABG/logs/z06f_new.out 2> /prj/aquacat/CodeABG/logs/z06f_new.err &

echo "RCM01 present running"
bash /prj/aquacat/CodeABG/BASH/456i_netcdf_script.sh 07 future > /prj/aquacat/CodeABG/logs/z07f_new.out 2> /prj/aquacat/CodeABG/logs/z07f_new.err &

wait

echo "RCM01 present running"
bash /prj/aquacat/CodeABG/BASH/456i_netcdf_script.sh 08 present > /prj/aquacat/CodeABG/logs/z08p_new.out 2> /prj/aquacat/CodeABG/logs/z08p_new.err &

echo "RCM01 present running"
bash /prj/aquacat/CodeABG/BASH/456i_netcdf_script.sh 09 present > /prj/aquacat/CodeABG/logs/z09p_new.out 2> /prj/aquacat/CodeABG/logs/z09p_new.err &

echo "RCM01 present running"
bash /prj/aquacat/CodeABG/BASH/456i_netcdf_script.sh 10 present > /prj/aquacat/CodeABG/logs/z10p_new.out 2> /prj/aquacat/CodeABG/logs/z10p_new.err &

wait

echo "RCM01 present running"
bash /prj/aquacat/CodeABG/BASH/456i_netcdf_script.sh 08 future > /prj/aquacat/CodeABG/logs/z08f_new.out 2> /prj/aquacat/CodeABG/logs/z08f_new.err &

echo "RCM01 present running"
bash /prj/aquacat/CodeABG/BASH/456i_netcdf_script.sh 09 future > /prj/aquacat/CodeABG/logs/z09f_new.out 2> /prj/aquacat/CodeABG/logs/z09f_new.err &

echo "RCM01 present running"
bash /prj/aquacat/CodeABG/BASH/456i_netcdf_script.sh 10 future > /prj/aquacat/CodeABG/logs/z10f_new.out 2> /prj/aquacat/CodeABG/logs/z10f_new.err &

wait

date -R
echo "Done."