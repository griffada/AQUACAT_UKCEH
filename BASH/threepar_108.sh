#!/bin/bash

# This bash script runs the bare essentials to get to 108 for RCM01 present
# In theory we have all these already, except for 107, where new 

echo $1
echo $2

echo "102 running"
Rscript --vanilla 102_Threshold_Extract.R $1 $2 > ./logs/"$1"_"$2"_102.out 2> ./logs/"$1"_"$2"_102.err
tail ./logs/"$1"_"$2"_102.out

echo " "
echo "**********"
echo " "

#echo "103 running"
#Rscript --vanilla 103_Regional_Splitting.R > ./logs/"$1"_"$2"_103.out 2> ./logs/"$1"_"$2"_103.err 
#tail ./logs/103.out

#echo ""
#echo "**********"
#echo ""

echo "104 running"
Rscript --vanilla 104_Event_Extract.R $1 $2 > ./logs/"$1"_"$2"_104.out 2> ./logs/"$1"_"$2"_104.err
tail ./logs/"$1"_"$2"_104.out

echo ""
echo "**********"
echo ""

echo "105 running"
Rscript --vanilla 105_Event_Summary.R $1 $2 > ./logs/"$1"_"$2"_105.out 2> ./logs/"$1"_"$2"_105.err
tail ./logs/"$1"_"$2"_105.out

echo " "
echo "**********"
echo " "

echo "106 running"
Rscript --vanilla 106_PoE_Estimation.R $1 $2 > ./logs/"$1"_"$2"_106.out 2> ./logs/"$1"_"$2"_106.err
tail ./logs/"$1"_"$2"_106.out

echo " "
echo "**********"
echo " "

echo "107 running"
Rscript --vanilla 107_Empirical_Copula.R $1 $2 > ./logs/"$1"_"$2"_107.out 2> ./logs/"$1"_"$2"_107.err 
tail ./logs/"$1"_"$2"_107.out

echo " "
echo "**********"
echo " "

echo "108 running"
Rscript --vanilla 108_Event_Splitting.R $1 $2 > ./logs/"$1"_"$2"_108.out 2> ./logs/"$1"_"$2"_108.err
tail ./logs/"$1"_"$2"_108.out

echo " "
echo "**********"
echo " "

echo "113 running"
Rscript --vanilla 113_proper_event_summary.R $1 $2 > ./logs/"$1"_"$2"_113.out 2> ./logs/"$1"_"$2"_113.err
tail ./logs/"$1"_"$2"_113.out

echo "Complete"
