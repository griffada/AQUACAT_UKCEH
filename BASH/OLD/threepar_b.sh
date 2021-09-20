#!/bin/bash

# This bash script runs the bare essentials to get to 108 for RCM01 present
# In theory we have all these already, except for 107, where new 

echo $1
echo $2

# echo "102 Threshold running"
# Rscript --vanilla 102_Threshold_Extract.R $1 $2 > ./logs/"$1"_"$2"_102.out 2> ./logs/"$1"_"$2"_102.err
# tail ./logs/"$1"_"$2"_102.out
# tail ./logs/"$1"_"$2"_102.err

# echo " "
# echo "**********"
# echo " "

echo "103 Regions running"
Rscript --vanilla 103_Regional_Splitting.R > ./logs/"$1"_"$2"_103.out 2> ./logs/"$1"_"$2"_103.err 
tail ./logs/103.out
echo ""
echo "**********"
echo ""

# echo "104 Event Extract running"
# Rscript --vanilla 104b_Event_Extract.R $1 $2 > ./logs/"$1"_"$2"_104b.out 2> ./logs/"$1"_"$2"_104b.err
# tail ./logs/"$1"_"$2"_104b.out
# tail ./logs/"$1"_"$2"_104b.err

# echo ""
# echo "**********"
# echo ""

# echo "105 OBS summary running"
# Rscript --vanilla 105_Event_Summary.R $1 $2 > ./logs/"$1"_"$2"_105.out 2> ./logs/"$1"_"$2"_105.err
# tail ./logs/"$1"_"$2"_105.out
# tail ./logs/"$1"_"$2"_105.err

# echo " "
# echo "**********"
# echo " "

# echo "106 PoE obs running"
# Rscript --vanilla 106b_PoE_Estimation.R $1 $2 > ./logs/"$1"_"$2"_106.out 2> ./logs/"$1"_"$2"_106.err
# tail ./logs/"$1"_"$2"_106.out
# tail ./logs/"$1"_"$2"_106.err

# echo " "
# echo "**********"
# echo " "

# echo "107 EC dpe running"
# Rscript --vanilla 107_Empirical_Copula.R $1 $2 > ./logs/"$1"_"$2"_107.out 2> ./logs/"$1"_"$2"_107.err 
# tail ./logs/"$1"_"$2"_107.out
# tail ./logs/"$1"_"$2"_107.err

# echo " "
# echo "**********"
# echo " "

echo "107 EC obs running"
Rscript --vanilla 107b_Empirical_Copula.R $1 $2 > ./logs/"$1"_"$2"_107b.out 2> ./logs/"$1"_"$2"_107b.err 
tail ./logs/"$1"_"$2"_107b.out
tail ./logs/"$1"_"$2"_107b.err

echo " "
echo "**********"
echo " "

echo "108 OBS EC event splitting running"
Rscript --vanilla 108_Event_Splitting.R $1 $2 > ./logs/"$1"_"$2"_108.out 2> ./logs/"$1"_"$2"_108.err
tail ./logs/"$1"_"$2"_108.out

echo " "
echo "**********"
echo " "

echo "113 OBS event summary running"
Rscript --vanilla 113_proper_event_summary.R $1 $2 > ./logs/"$1"_"$2"_113.out 2> ./logs/"$1"_"$2"_113.err
tail ./logs/"$1"_"$2"_113.out
tail ./logs/"$1"_"$2"_113.err

echo "113 EC event summary running"
Rscript --vanilla 113b_EC_proper_event_summary.R $1 $2 > ./logs/"$1"_"$2"_113b.out 2> ./logs/"$1"_"$2"_113b.err
tail ./logs/"$1"_"$2"_113b.out
tail ./logs/"$1"_"$2"_113b.err

echo "Complete"
