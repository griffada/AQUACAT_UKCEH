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
Rscript --vanilla 102_Threshold_Extract.R $EN $PE > pipe102_pr.out 2> pipe102_pr.out
tail pipe102_pr.out

echo "\n"
echo "**********"
echo "\n"

echo "103 01 present"
Rscript --vanilla 103_Regional_Splitting.R > pipe103.out 2> pipe103.out 
tail pipe103.out

echo "\n"
echo "**********"
echo "\n"

echo "104 01 present"
Rscript --vanilla 104_Event_Extract.R $EN $PE > pipe104_pr.out 2> pipe104_pr.out 
tail pipe104_pr.out

echo "\n"
echo "**********"
echo "\n"

echo "105 01 present"
Rscript --vanilla 105_Event_Summary.R $EN $PE > pipe105_pr.out 2> pipe105_pr.out 
tail pipe105_pr.out

echo "\n"
echo "**********"
echo "\n"

echo "106 01 present"
Rscript --vanilla 106_PoE_Estimation.R $EN $PE > pipe106_pr.out 2> pipe106_pr.out 
tail pipe106_pr.out

echo "\n"
echo "**********"
echo "\n"

echo "107 01 present"
Rscript --vanilla 107_Empirical_Copula.R $EN $PE > pipe107_pr.out 2> pipe107_pr.out 
tail pipe107_pr.out

echo "\n"
echo "**********"
echo "\n"

echo "108 01 present"
Rscript --vanilla 108_Event_Splitting.R $EN $PE $RE > pipe108_pr.out 2> pipe108_pr.out 
tail pipe108_pr.out

echo "\n"
echo "**********"
echo "\n"

echo "Complete"
