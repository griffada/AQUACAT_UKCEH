#!/bin/bash

# This bash script runs the bare essentials to get from 109 to 111 for a specified region, RCM and period.

#SBATCH -p short-serial
#SBATCH -o present_%J.out
#SBATCH -e present_%J.err
#SBATCH --time 23:59:00
#SBATCH --mem=12000
#SBATCH --chdir=/home/users/adagri/AQUACAT/CodeABG

echo "Running 108-111"
module load jaspy/3.7/r20200606


echo "108 01 present"
Rscript --vanilla 108_Event_Splitting.R $EN $PE > ./logs/108_"$EN"_"$PE".out 2> ./logs/108_"$EN"_"$PE".err
tail ./logs/108_"$EN"_"$PE".out

echo "109 01 present"
Rscript --vanilla 109_HeffTawn_Modelling.R $EN $PE $RE > ./logs/109_"$EN"_"$PE"_"$RE".out 2> ./logs/109_"$EN"_"$PE"_"$RE".err
tail ./logs/109_"$EN"_"$PE"_"$RE".out

echo "110 01 present"
Rscript --vanilla 110_HT_PoEEstimation.R $EN $PE $RE > ./logs/110_"$EN"_"$PE"_"$RE".out 2> ./logs/110_"$EN"_"$PE"_"$RE".err
tail ./logs/110_"$EN"_"$PE"_"$RE".out

echo "111 01 present"
Rscript --vanilla 111_SQL_Compilation.R $EN $PE $RE > ./logs/111_"$EN"_"$PE"_"$RE".out 2> ./logs/111_"$EN"_"$PE"_"$RE".err
tail ./logs/111_"$EN"_"$PE"_"$RE".out

echo "Complete"
