#!/bin/bash

# This bash script runs the bare essentials to get to 108 for RCM01 present
# In theory we have all these already, except for 107, where new 

#SBATCH -p short-serial
#SBATCH -o present_%J.out
#SBATCH -e present_%J.err
#SBATCH --time 23:59:00
#SBATCH --mem=12000
#SBATCH --chdir=/home/users/adagri/AQUACAT/CodeABG

echo "Running 102-115 not using 103"
module load jaspy/3.7/r20200606

echo $EN
echo $PE
echo $RE

echo "108 running"
Rscript --vanilla 108_Event_Splitting.R $EN $PE $RE > ./logs/"$EN"_"$PE"_108.out 2> ./logs/"$EN"_"$PE"_108.err
tail ./logs/"$EN"_"$PE"_108.out
tail ./logs/"$EN"_"$PE"_108.err

echo " "
echo "**********"
echo " "

echo "109 running"
Rscript --vanilla 109_HeffTawn_Modelling.R $EN $PE $RE > ./logs/"$EN"_"$RE"_"$PE"_109.out 2> ./logs/"$EN"_"$RE"_"$PE"_109.err
tail ./logs/"$EN"_"$RE"_"$PE"_109.out
tail ./logs/"$EN"_"$RE"_"$PE"_109.err

echo "110 running"
Rscript --vanilla 110_HT_PoEEstimation.R $EN $PE $RE > ./logs/"$EN"_"$RE"_"$PE"_110.out 2> ./logs/"$EN"_"$RE"_"$PE"_110.err
tail ./logs/"$EN"_"$RE"_"$PE"_110.out
tail ./logs/"$EN"_"$RE"_"$PE"_110.err

echo "114 running"
Rscript --vanilla 114_regional_proper_event_summary.R $EN $PE $RE> ./logs/"$EN"_"$RE"_"$PE"_114.out 2> ./logs/"$EN"_"$RE"_"$PE"_114.err
tail ./logs/"$EN"_"$RE"_"$PE"_114.out
tail ./logs/"$EN"_"$RE"_"$PE"_114.err

echo "Complete"
