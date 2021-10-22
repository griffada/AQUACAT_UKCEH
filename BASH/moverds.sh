#!/bin/bash

# This script should rerun the EC PoE, HT PoE, and summary scripts

for r in 04 05 06 07 08 09 10 11 12 13 15 01
do
for p in "198012_201011" "205012_208011"
    do

    echo $r
    echo $p
	
	mkdir -vp /prj/aquacat/InterimData/RCM"$r"_"$p"/NW
	
	
	mv -v /prj/aquacat/Data/RCM"$r"_"$p"/NW/marginal_models.rds /prj/aquacat/InterimData/RCM"$r"_"$p"/NW/marginal_models.rds 
	mv -v /prj/aquacat/Data/RCM"$r"_"$p"/NW/coefficients_NW_POT2_pc01_RCM"$r"_"$p".rds /prj/aquacat/InterimData/RCM"$r"_"$p"/NW/coefficients_NW_POT2_pc01_RCM"$r"_"$p".rds
	mv -v /prj/aquacat/Data/RCM"$r"_"$p"/NW/depStructNW_POT2_pc01_RCM"$r"_"$p".rds /prj/aquacat/InterimData/RCM"$r"_"$p"/NW/depStructNW_POT2_pc01_RCM"$r"_"$p".rds
	mv -v /prj/aquacat/Data/RCM"$r"_"$p"/NW/zScoresNW_POT2_pc01_RCM"$r"_"$p".rds /prj/aquacat/InterimData/RCM"$r"_"$p"/NW/zScoresNW_POT2_pc01_RCM"$r"_"$p".rds
	mv -v /prj/aquacat/Data/RCM"$r"_"$p"/NW/transformedData.rds /prj/aquacat/InterimData/RCM"$r"_"$p"/NW/transformedData.rds 
	
	done
done
echo "Complete."