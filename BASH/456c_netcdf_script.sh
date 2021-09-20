#!/bin/bash

# This bash script runs the bare essentials to get to 108 for RCM01 present
# In theory we have all these already, except for 107, where new 

echo "101 running"
Rscript --vanilla 101_YAML_creation.R $1 $2 > ./logs/"$1"_"$2"_101.out 2> ./logs/"$1"_"$2"_101.err
tail ./logs/"$1"_"$2"_101.out
tail ./logs/"$1"_"$2"_101.err

echo "102 running"
Rscript --vanilla 102e_Threshold_Extract.R $1 $2 > ./logs/"$1"_"$2"_102.out 2> ./logs/"$1"_"$2"_102.err
tail ./logs/"$1"_"$2"_102.out
tail ./logs/"$1"_"$2"_102.err

echo "104 running"
Rscript --vanilla 104b_Event_Extract.R $1 $2 > ./logs/"$1"_"$2"_104.out 2> ./logs/"$1"_"$2"_104.err
tail ./logs/"$1"_"$2"_104.out
tail ./logs/"$1"_"$2"_104.err

echo "105N running"
Rscript --vanilla 105N_Event_Summary.R $1 $2 > ./logs/"$1"_"$2"_105N.out 2> ./logs/"$1"_"$2"_105N.err
tail ./logs/"$1"_"$2"_105N.out
tail ./logs/"$1"_"$2"_105N.err

echo "106N running"
Rscript --vanilla 106N_PoE_EstimationB.R $1 $2 > ./logs/"$1"_"$2"_106N.out 2> ./logs/"$1"_"$2"_106N.err
tail ./logs/"$1"_"$2"_106N.out
tail ./logs/"$1"_"$2"_106N.err 

echo "107N running"
Rscript --vanilla 107N_EmpiricalBeta_CopulaB.R $1 $2 > ./logs/"$1"_"$2"_107N.out 2> ./logs/"$1"_"$2"_107N.err
tail ./logs/"$1"_"$2"_107N.out
tail ./logs/"$1"_"$2"_107N.err

echo "117N running"
Rscript --vanilla 117N_EmpiricalBeta_PoEB.R $1 $2 > ./logs/"$1"_"$2"_117N.out 2> ./logs/"$1"_"$2"_117N.err
tail ./logs/"$1"_"$2"_117N.out
tail ./logs/"$1"_"$2"_117N.err

# echo "108N running"
# Rscript --vanilla 108N_Event_Splitting.R $1 $2 > ./logs/"$1"_"$2"_108N.out 2> ./logs/"$1"_"$2"_108N.err
# tail ./logs/"$1"_"$2"_108N.out
# tail ./logs/"$1"_"$2"_108N.err

# echo "113N running"
# python3 -u 113N_proper_event_summary.py $1 $2 > ./logs/"$1"_"$2"_113N.out 2> ./logs/"$1"_"$2"_113N.err
# tail ./logs/"$1"_"$2"_113N.out
# tail ./logs/"$1"_"$2"_113N.err

 #echo "113BN running"
 #python3 -u 113bN_EC_proper_event_summary.py $1 $2 > ./logs/"$1"_"$2"_113BN.out 2> ./logs/"$1"_"$2"_113BN.err
# tail ./logs/"$1"_"$2"_113BN.out
# tail ./logs/"$1"_"$2"_113BN.err

# echo "109N running"
# Rscript --vanilla 109N_HeffTawn_Modelling.R $1 $2 > ./logs/"$1"_"$2"_109N.out 2> ./logs/"$1"_"$2"_109N.err
# tail ./logs/"$1"_"$2"_109N.out
# tail ./logs/"$1"_"$2"_109N.err

# echo "110N running"
# Rscript --vanilla 110N_HT_PoEEstimation.R $1 $2 NW > ./logs/"$1"_"$2"_NW_110N.out 2> ./logs/"$1"_"$2"_NW_110N.err
# tail ./logs/"$1"_"$2"_NW_110N.out
# tail ./logs/"$1"_"$2"_NW_110N.err

# echo "114N running"
# python3 -u 114N_regional_proper_event_summary.py $1 $2 > ./logs/"$1"_"$2"_NW_114N.out 2> ./logs/"$1"_"$2"_NW_114N.err
# tail ./logs/"$1"_"$2"_NW_114N.out
# tail ./logs/"$1"_"$2"_NW_114N.err

# echo "119N running"
# Rscript --vanilla 119_HT_tidying.R $1 $2 NW > ./logs/"$1"_"$2"_NW_119N.out 2> ./logs/"$1"_"$2"_NW_119N.err
# tail ./logs/"$1"_"$2"_NW_119N.out
# tail ./logs/"$1"_"$2"_NW_119N.err

echo "Complete"
