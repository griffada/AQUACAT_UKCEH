#!/bin/bash

# This script should rerun the summary scripts

for r in 01 04 05 06 07 08 09 10 11 12 13 15
do
for p in present future 
    do

    echo $r
    echo $p

    echo "114 event summary running"
    python3 -u 114N_regional_proper_event_summary.py $r $p > ./logs/"$r"_"$p"_114.out 2> ./logs/"$r"_"$p"_114.err
    tail ./logs/"$r"_"$p"_114.out

    done
done
echo "Complete"