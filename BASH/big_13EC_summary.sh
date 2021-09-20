#!/bin/bash

# This script should rerun the summary scripts

for r in 01 04 05 06 07 08 09 10 11 12 13 15
do
for p in present
    do

    echo $r
    echo $p

    echo "113 OBS event summary running"
    python3 -u 113bN_EC_proper_event_summary.py $r $p > ./logs/"$r"_"$p"_113.out 2> ./logs/"$r"_"$p"_113.err
    tail ./logs/"$r"_"$p"_113.out

    # Rscript --vanilla 113bN_EC_proper_event_summary_after_python.R $r $p > ./logs/"$r"_"$p"_113b.out 2> ./logs/"$r"_"$p"_113b.err
    # tail ./logs/"$r"_"$p"_113b.out

    done
done
echo "Complete"