#!/bin/bash

# This bash script runs the bare essentials to get to 108 for RCM01 present
# In theory we have all these already, except for 107, where new 

echo "EC2 running"

for i in 1 2 3 4 5
do

    Rscript --vanilla EMPCOP3.R $1 $2 > ./logs/"$1"_"$2"_EC2.out 2> ./logs/"$1"_"$2"_EC2.err
    tail ./logs/"$1"_"$2"_EC2.out
    tail ./logs/"$1"_"$2"_EC2.err

done

echo "Complete"
