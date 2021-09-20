#!/bin/bash

# This script should rerun the EC PoE, HT PoE, and summary scripts

for r in 01 04 05 06 07 08 09 10 11 12 13 15
do
for p in present future
    do

    echo $r
    echo $p

    python3 -u 108N_Event_Splitting.py $r $p> ./logs/"$r"_"$p"_108.out 2> ./logs/"$r"_"$p"_108N.out

    done
done
echo "Complete"