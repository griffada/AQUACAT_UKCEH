#!/bin/bash

# This script should rerun the summary scripts

for r in 01 04 05 06 07 08 09 10 11 12 13 15
do
for p in present future
    do

    echo $r
    echo $p
    
    # echo "108 event splitting"
    # Rscript --vanilla 108_Event_Splitting.R $r $p > ./logs/"$r"_"$p"_108_rerun.out 2> ./logs/"$r"_"$p"_108_rerun.err
    # tail ./logs/"$r"_"$p"_108_rerun.err

    echo "113 OBS event summary running"
    Rscript --vanilla 113N_proper_event_summary.R $r $p > ./logs/"$r"_"$p"_113_rerun.out 2> ./logs/"$r"_"$p"_113_rerun.err
    tail ./logs/"$r"_"$p"_113_rerun.err

    # echo "113 EC event summary running"
    # Rscript --vanilla 113b_EC_proper_event_summary.R $r $p > ./logs/"$r"_"$p"_113b_rerun.out 2> ./logs/"$r"_"$p"_113b_rerun.err
    # tail ./logs/"$r"_"$p"_113b_rerun.err

    # echo "114 regional event summary running"
    # Rscript --vanilla 114_regional_proper_event_summary.R $r $p NW > ./logs/"$r"_"$p"_114_rerun.out 2> ./logs/"$r"_"$p"_114_rerun.err
    # tail ./logs/"$r"_"$p"_114_rerun.err

    done
done
echo "Complete"