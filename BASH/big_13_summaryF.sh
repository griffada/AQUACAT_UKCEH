#!/bin/bash

# This script should rerun the summary scripts

for r in 04 05 06 07 08 09 10 11 12 13
do
for p in present future
    do

    echo $r
    echo $p

    echo "113 OBS event summary running"
    python3 -u 113bN_proper_event_summary.py $r $p > ./logs/"$r"_"$p"_113.out 2> ./logs/"$r"_"$p"_113.err
    tail ./logs/"$r"_"$p"_113.out

    # Rscript --vanilla 113aN_proper_event_summary_after_python.R $r $p > ./logs/"$r"_"$p"_113a.out 2> ./logs/"$r"_"$p"_113a.err
    # tail ./logs/"$r"_"$p"_113a.out

    done
done
echo "Complete"