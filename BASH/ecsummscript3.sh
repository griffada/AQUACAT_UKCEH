#!/bin/bash

# This script should rerun the summary scripts

for r in 01 04 05 06 07 08 09 10 11 12 13 15
do

echo $r
echo $1

python3 -u ./113bN_EC_proper_event_summary.py $r $1

done
echo "Complete"
