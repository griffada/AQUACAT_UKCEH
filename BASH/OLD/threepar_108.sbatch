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

echo $EN
echo $PE

echo "102 Threshold running"
Rscript --vanilla 102_Threshold_Extract.R $EN $PE > ./logs/"$EN"_"$PE"_102.out 2> ./logs/"$EN"_"$PE"_102.err
tail ./logs/"$EN"_"$PE"_102.out

echo " "
echo "**********"
echo " "

echo "103 Regions running"
Rscript --vanilla 103_Regional_Splitting.R > ./logs/"$EN"_"$PE"_103.out 2> ./logs/"$EN"_"$PE"_103.err 
tail ./logs/103.out
echo ""
echo "**********"
echo ""

echo "104 Event Extract running"
Rscript --vanilla 104b_Event_Extract.R $EN $PE > ./logs/"$EN"_"$PE"_104b.out 2> ./logs/"$EN"_"$PE"_104b.err
tail ./logs/"$EN"_"$PE"_104b.out

echo ""
echo "**********"
echo ""

echo "105 OBS summary running"
Rscript --vanilla 105_Event_Summary.R $EN $PE > ./logs/"$EN"_"$PE"_105.out 2> ./logs/"$EN"_"$PE"_105.err
tail ./logs/"$EN"_"$PE"_105.out

echo " "
echo "**********"
echo " "

echo "106 PoE obs running"
Rscript --vanilla 106_PoE_Estimation.R $EN $PE > ./logs/"$EN"_"$PE"_106.out 2> ./logs/"$EN"_"$PE"_106.err
tail ./logs/"$EN"_"$PE"_106.out

echo " "
echo "**********"
echo " "

echo "107 EC dpe running"
Rscript --vanilla 107_Empirical_Copula.R $EN $PE > ./logs/"$EN"_"$PE"_107.out 2> ./logs/"$EN"_"$PE"_107.err 
tail ./logs/"$EN"_"$PE"_107.out

echo " "
echo "**********"
echo " "

echo "107 EC obs running"
Rscript --vanilla 107b_Empirical_Copula.R $EN $PE > ./logs/"$EN"_"$PE"_107b.out 2> ./logs/"$EN"_"$PE"_107b.err 
tail ./logs/"$EN"_"$PE"_107.out

echo " "
echo "**********"
echo " "

echo "108 OBS EC event splitting running"
Rscript --vanilla 108_Event_Splitting.R $EN $PE > ./logs/"$EN"_"$PE"_108.out 2> ./logs/"$EN"_"$PE"_108.err
tail ./logs/"$EN"_"$PE"_108.out

echo " "
echo "**********"
echo " "

echo "113 OBS event summary running"
Rscript --vanilla 113_proper_event_summary.R $EN $PE > ./logs/"$EN"_"$PE"_113.out 2> ./logs/"$EN"_"$PE"_113.err
tail ./logs/"$EN"_"$PE"_113.out

echo "Complete"
