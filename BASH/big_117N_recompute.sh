#!/bin/bash

# This script should rerun the EC PoE, HT PoE, and summary scripts

for r in 04 05 06 07 08 09 10 11 12 13 15
do
for p in present future
    do

    echo $r
    echo $p

    # echo "117N regional event summary running"
    # Rscript --vanilla 117N_EmpiricalBeta_PoE.R $r $p > ./logs/"$r"_"$p"_117N_rerun.out 2> ./logs/"$r"_"$p"_117N_rerun.err
    # tail ./logs/"$r"_"$p"_114_rerun.err

    # echo "110N running"
    # Rscript --vanilla 110N_HT_PoEEstimation.R $r $p NW > ./logs/"$1"_"$2"_NW_110N_rerun.out 2> ./logs/"$1"_"$2"_NW_110N_rerun.err
    # tail ./logs/"$1"_"$2"_NW_110N_rerun.err

    # echo "113N running"
    # Rscript --vanilla 113N_proper_event_summary.R $1 $2 > ./logs/"$1"_"$2"_113N.out 2> ./logs/"$1"_"$2"_113N.err
    # tail ./logs/"$1"_"$2"_113N.out
    # tail ./logs/"$1"_"$2"_113N.err

    echo "113BN running"
    python3 ecsummpy.py $r $p > ./logs/"$r"_"$p"_113BN.out 2> ./logs/"$r"_"$p"_113BN.err
    tail ./logs/"$r"_"$p"_113BN.out
    tail ./logs/"$r"_"$p"_113BN.err

    # echo "114N running"
    # Rscript --vanilla 114N_regional_proper_event_summary.R $1 $2 NW > ./logs/"$1"_"$2"_NW_114N.out 2> ./logs/"$1"_"$2"_NW_114N.err
    # tail ./logs/"$1"_"$2"_114N.out
    # tail ./logs/"$1"_"$2"_114N.err

    done
done
echo "Complete"