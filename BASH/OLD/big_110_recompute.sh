#!/bin/bash

# This script should rerun the summary scripts

for r in 01 04 05 06 07 08 09 10 11 12 13 15
do
for p in present future
    do

    echo $r
    echo $p

    # echo "108 eventNumbering"
    # Rscript --vanilla 108bN_eventNumbering.R $r $p > ./logs/"$r"_"$p"_108b_again.out 2> ./logs/"$r"_"$p"_108b_again.err

     echo "110 event PoE again"
     Rscript --vanilla 110N_HT_PoEEstimation.R $r $p NW > ./logs/"$r"_"$p"_110_rerun.out

    echo "114 regional event summary running"
    Rscript --vanilla 114N_regional_proper_event_summary.R $r $p NW > ./logs/"$r"_"$p"_114_rerun.out 2> ./logs/"$r"_"$p"_114_rerun.err
    tail ./logs/"$r"_"$p"_114_rerun.err

    done
done
echo "Complete"
