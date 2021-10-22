#!/bin/bash

# This script should rerun the summary scripts

for r in 01 04 05 06 07 08 09 10 11 12 13 15
do
for p in present
    do

    echo $r
    echo $p

	echo "110N running"
	Rscript --vanilla 110N_HT_PoEEstimation.R $r $p NW > ./logs/"$r"_"$p"_NW_110N.out 2> ./logs/"$r"_"$p"_NW_110N.err
	tail ./logs/"$r"_"$p"_NW_110N.out
	tail ./logs/"$r"_"$p"_NW_110N.err

    done
done
echo "Complete"