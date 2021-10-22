#!/bin/bash

for r in 01 04 05 06 07 08 09 10 11 12 13 15
do
for p in "198012_201011" "205012_208011"
    do
	
	if [[ "$r" == "04" && "$p" == "198012_201011" ]]
	then
		continue
	else
	
    echo $r
    echo $p
	
	mv /prj/aquacat/Data/RCM"$r"_"$p"/NW/eventHT_region_NW_RCM"$r"_"$p$".nc /prj/aquacat/Data/RCM"$r"_"$p"/NW/eventHT_region_NW_RCM"$r"_"$p".nc
	fi
	
	done
done