#!/bin/bash

# This bash script runs the bare essentials to get to 108 for RCM01 present
# In theory we have all these already, except for 107, where new 


echo "102 01 present"
Rscript --vanilla 102_Threshold_Extract.R $1 $2 > pipe102_"$1"_"$2"_pr.out 2> pipe102_"$1"_"$2"_pr.out
tail pipe102_"$1"_"$2"_pr.out

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
Rscript --vanilla 104_Event_Extract.R $1 $2 > pipe104_"$1"_"$2"_pr.out 2> pipe104_"$1"_"$2"_pr.out 
tail pipe104_"$1"_"$2"_pr.out

echo "\n"
echo "**********"
echo "\n"

echo "105 01 present"
Rscript --vanilla 105_Event_Summary.R $1 $2 > pipe105_"$1"_"$2"_pr.out 2> pipe105_"$1"_"$2"_pr.out 
tail pipe105_"$1"_"$2"_pr.out

echo "\n"
echo "**********"
echo "\n"

echo "106 01 present"
Rscript --vanilla 106_PoE_Estimation.R $1 $2 > pipe106_"$1"_"$2"_pr.out 2> pipe106_"$1"_"$2"_pr.out 
tail pipe106_"$1"_"$2"_pr.out

echo "\n"
echo "**********"
echo "\n"

echo "107 01 present"
Rscript --vanilla 107_Empirical_Copula.R $1 $2 > pipe107_"$1"_"$2"_pr.out 2> pipe107_"$1"_"$2"_pr.out 
tail pipe107_"$1"_"$2"_pr.out

echo "\n"
echo "**********"
echo "\n"

echo "108 01 present"
Rscript --vanilla 108_Event_Splitting.R $1 $2 > pipe108_"$1"_"$2"_pr.out 2> pipe108_"$1"_"$2"_pr.out 
tail pipe108_"$1"_"$2"_pr.out

echo "Complete"
