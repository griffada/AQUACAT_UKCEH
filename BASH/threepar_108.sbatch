#!/bin/bash

# This bash script runs the bare essentials to get to 108 for RCM01 present
# In theory we have all these already, except for 107, where new 

#SBATCH -p short-serial
#SBATCH -o present_%J.out
#SBATCH -e present_%J.err
#SBATCH --time 23:59:00
#SBATCH --mem=12000
#SBATCH --chdir=/home/users/adagri/AQUACAT/CodeABG

echo "Running slimline_testing6"
module load jaspy/3.7/r20200606

echo "102 01 present"
Rscript --vanilla 102_Threshold_Extract.R $EN $PE > ./logs/102_"$EN"_"$PE".out 2> ./logs/102_"$EN"_"$PE".err
tail ./logs/102_"$EN"_"$PE".out

echo " "
echo "**********"
echo " "

echo "103 01 present"
Rscript --vanilla 103_Regional_Splitting.R > ./logs/103.out 2> ./logs/pipe103.err 
tail ./logs/103.out

echo 
echo "**********"
echo 

echo "104 01 present"
Rscript --vanilla 104_Event_Extract.R $EN $PE > ./logs/104_"$EN"_"$PE".out 2> ./logs/104_"$EN"_"$PE".err
tail ./logs/104_"$EN"_"$PE".out

echo ""
echo "**********"
echo ""

echo "105 01 present"
Rscript --vanilla 105_Event_Summary.R $EN $PE > ./logs/105_"$EN"_"$PE".out 2> ./logs/105_"$EN"_"$PE".err
tail ./logs/105_"$EN"_"$PE".out

echo " "
echo "**********"
echo " "

echo "106 01 present"
Rscript --vanilla 106_PoE_Estimation.R $EN $PE > ./logs/106_"$EN"_"$PE".out 2> ./logs/106_"$EN"_"$PE".err
tail ./logs/106_"$EN"_"$PE".out

echo " "
echo "**********"
echo " "

echo "107 01 present"
Rscript --vanilla 107_Empirical_Copula.R $EN $PE > ./logs/107_"$EN"_"$PE".out 2> ./logs/107_"$EN"_"$PE".err 
tail ./logs/107_"$EN"_"$PE".out

echo " "
echo "**********"
echo " "

echo "108 01 present"
Rscript --vanilla 108_Event_Splitting.R $EN $PE > ./logs/108_"$EN"_"$PE".out 2> ./logs/108_"$EN"_"$PE".err
tail ./logs/108_"$EN"_"$PE".out

echo "Complete"
