#!/bin/bash

# This script should rerun the summary scripts

for r in 01 04 05 06 07 08 09 10 11 12 13 15
do
for p in present future
    do

    echo $r
    echo $p

	echo "108N running"
	python3 -u 108N_Event_Splitting.py $r $p > ./logs/"$r"_"$p"_108N.out 2> ./logs/"$r"_"$p"_108N.err
	tail ./logs/"$r"_"$p"_108N.out
	tail ./logs/"$r"_"$p"_108N.err

    done
done
echo "Complete"