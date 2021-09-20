#!/bin/bash

# This script should rerun the summary scripts

for r in 01 04 05 06 07 08 09 10 11 12 13 15
do
for p in present future
    do

    echo $r
    echo $p

    echo "114 regional event summary running"
    Rscript --vanilla 114_regional_proper_event_summary.R $r $p NW > ./logs/"$r"_"$p"_114_rerun.out 2> ./logs/"$r"_"$p"_114_rerun.err
    tail ./logs/"$r"_"$p"_114_rerun.err

    done
done
echo "Complete"