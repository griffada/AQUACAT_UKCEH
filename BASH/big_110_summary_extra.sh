#!/bin/bash

# This script should rerun the summary scripts

for r in 01 04 05 06 07 08 09 10 11 12 13 15
do
for p in present future
    do

    echo $r
    echo $p

	echo "110NOBS running"
	Rscript --vanilla 110N_OBS_PoEEstimation.R $r $p NW > ./logs/"$r"_"$p"_NW_110O.out 2> ./logs/"$r"_"$p"_NW_110O.err
	tail ./logs/"$r"_"$p"_NW_110O.out
	tail ./logs/"$r"_"$p"_NW_110O.err
	
	echo "110NEC running"
	Rscript --vanilla 110N_OBS_PoEEstimation.R $r $p NW > ./logs/"$r"_"$p"_NW_110E.out 2> ./logs/"$r"_"$p"_NW_110E.err
	tail ./logs/"$r"_"$p"_NW_110E.out
	tail ./logs/"$r"_"$p"_NW_110E.err

    done
done
echo "Complete"