#!/bin/bash

# This script should rerun the summary scripts

for r in 09 10 11 12 13 15
do
for p in future 
    do

    echo $r
    echo $p

    echo "113 OBS event summary running"
    python3 -u 113bN_proper_event_summary.py $r $p > ./logs/"$r"_"$p"_113.out 2> ./logs/"$r"_"$p"_113.err
    tail ./logs/"$r"_"$p"_113.out

    echo "113 OBS event summary running"
    python3 -u 113bN_proper_event_summary.py $r $p FF > ./logs/"$r"_"$p"_113ff.out 2> ./logs/"$r"_"$p"_113ff.err
    tail ./logs/"$r"_"$p"_113.out

        echo "113 OBS event summary running"
    python3 -u 113bN_EC_proper_event_summary.py $r $p > ./logs/"$r"_"$p"_113ec.out 2> ./logs/"$r"_"$p"_113ec.err
    tail ./logs/"$r"_"$p"_113.out

    echo "113 OBS event summary running"
    python3 -u 113bN_EC_proper_event_summary.py $r $p FF > ./logs/"$r"_"$p"_113ecff.out 2> ./logs/"$r"_"$p"_113ecff.err
    tail ./logs/"$r"_"$p"_113.out

    done
done
echo "Complete"