#!/bin/bash

# This bash script reruns just the event splitting 108

for period in present future
do
    for rcm in 01 04 05 06 07 08 09 10 11 12 13 15
    do
    echo $rcm
    echo $period
    Rscript --vanilla 111_SQL_Compilation.R $rcm $period > ./logs/"$rcm"_"$period"_111.out 2> ./logs/"$rcm"_"$period"_111.err
    tail ./logs/"$rcm"_"$period"_111.out
    done
done