#!/bin/bash

# This bash script runs the bare essentials to get to 108 for RCM01 present
# In theory we have all these already, except for 107, where new 


echo "102 01 present"
Rscript --vanilla 102_Threshold_Extract.R $1 $2 > ./logs/102_"$1"_"$2".out 2> ./logs/102_"$1"_"$2".err
tail ./logs/102_"$1"_"$2".out

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
Rscript --vanilla 104_Event_Extract.R $1 $2 > ./logs/104_"$1"_"$2".out 2> ./logs/104_"$1"_"$2".err
tail ./logs/104_"$1"_"$2".out

echo ""
echo "**********"
echo ""

echo "105 01 present"
Rscript --vanilla 105_Event_Summary.R $1 $2 > ./logs/105_"$1"_"$2".out 2> ./logs/105_"$1"_"$2".err
tail ./logs/105_"$1"_"$2".out

echo " "
echo "**********"
echo " "

echo "106 01 present"
Rscript --vanilla 106_PoE_Estimation.R $1 $2 > ./logs/106_"$1"_"$2".out 2> ./logs/106_"$1"_"$2".err
tail ./logs/106_"$1"_"$2".out

echo " "
echo "**********"
echo " "

echo "107 01 present"
Rscript --vanilla 107_Empirical_Copula.R $1 $2 > ./logs/107_"$1"_"$2".out 2> ./logs/107_"$1"_"$2".err 
tail ./logs/107_"$1"_"$2".out

echo " "
echo "**********"
echo " "

echo "108 01 present"
Rscript --vanilla 108_Event_Splitting.R $1 $2 > ./logs/108_"$1"_"$2".out 2> ./logs/108_"$1"_"$2".err
tail ./logs/108_"$1"_"$2".out

echo "Complete"
