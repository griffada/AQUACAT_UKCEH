#!/bin/bash

# This script should rerun the EC PoE, HT PoE, and summary scripts

for r in 01 04 05 06 07 08 09 10 11 12 13 15
do
for p in present
    do

    echo $r
    echo $p

	echo "106cN running"
    Rscript --vanilla 106cN_PoE_Estimation.R $r $p > ./logs/"$r"_"$p"_106cN_rerun.out 2> ./logs/"$r"_"$p"_106cN_rerun.err
    tail ./logs/"$r"_"$p"_106cN_rerun.err

    echo "117Nrunning"
    Rscript --vanilla 117N_EmpiricalBeta_PoE.R $r $p > ./logs/"$r"_"$p"_117N_rerun.out 2> ./logs/"$r"_"$p"_117N_rerun.err
    tail ./logs/"$r"_"$p"_117N_rerun.err
	
    echo "113 OBS event summary running"
    python3 -u 113bN_proper_event_summary.py $r $p > ./logs/"$r"_"$p"_113.out 2> ./logs/"$r"_"$p"_113.err
    tail ./logs/"$r"_"$p"_113.err

    echo "113 EC event summary running"
    python3 -u 113bN_EC_proper_event_summary.py $r $p > ./logs/"$r"_"$p"_113ec.out 2> ./logs/"$r"_"$p"_113ec.err
    tail ./logs/"$r"_"$p"_113ec.err
	
    done
done
echo "Complete"