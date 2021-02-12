#!/bin/bash

# This bash script runs the bare essentials to get to 108 for RCM01 present
# In theory we have all these already, except for 107, where new 

echo "109 01 present"
Rscript --vanilla 109_HeffTawn_Modelling.R $1 $2 $3> ./logs/"$1"_"$3"_"$2"_109.out 2> ./logs/"$1"_"$3"_"$2"_109.err
tail ./logs/"$1"_"$3"_"$2"_109.out

echo "110 01 present"
Rscript --vanilla 110_HT_PoEEstimation.R $1 $2 > ./logs/"$1"_"$3"_"$2"_110.out 2> ./logs/"$1"_"$3"_"$2"_110.err
tail ./logs/"$1"_"$3"_"$2"_110.out

echo "111 01 present"
Rscript --vanilla 111_SQL_Compilation.R $1 $2 > ./logs/"$1"_"$3"_"$2"_111.out 2> ./logs/"$1"_"$3"_"$2"_111.err
tail ./logs/"$1"_"$3"_"$2"_111.out

echo "Complete"
