#! /bin/bash

for i in 01 04 05 06 07 08 09 10 11 12 13 15
do

    sh BASH/101_119_netcdf_script.sh $i future NW > "$i"future.out 2> "$i"future.err

done