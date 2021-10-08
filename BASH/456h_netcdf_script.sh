#!/bin/bash

# This bash script runs the bare essentials to get to 108 for RCM01 present
# In theory we have all these already, except for 107, where new 

if [ $2 = "present" ]; then
    exit 1
fi

echo "101 running"
Rscript --vanilla 101_YAML_creation.R $1 $2 > ./logs/"$1"_"$2"_101.out 2> ./logs/"$1"_"$2"_101.err
tail ./logs/"$1"_"$2"_101.out
tail ./logs/"$1"_"$2"_101.err

echo "102 running"
Rscript --vanilla 102g_Threshold_Extract.R $1 $2 > ./logs/"$1"_"$2"_102.out 2> ./logs/"$1"_"$2"_102.err
tail ./logs/"$1"_"$2"_102.out
tail ./logs/"$1"_"$2"_102.err

echo "104 running"
Rscript --vanilla 104b_Event_Extract.R $1 $2 > ./logs/"$1"_"$2"_104.out 2> ./logs/"$1"_"$2"_104.err
tail ./logs/"$1"_"$2"_104.out
tail ./logs/"$1"_"$2"_104.err

echo "105N running"
Rscript --vanilla 105N_Event_Summary.R $1 $2 > ./logs/"$1"_"$2"_105.out 2> ./logs/"$1"_"$2"_105.err
tail ./logs/"$1"_"$2"_105.out
tail ./logs/"$1"_"$2"_105.err

echo "106N running"
Rscript --vanilla 106cN_PoE_Estimation.R $1 $2 > ./logs/"$1"_"$2"_106.out 2> ./logs/"$1"_"$2"_106.err
tail ./logs/"$1"_"$2"_106.out
tail ./logs/"$1"_"$2"_106.err 

echo "107N running"
Rscript --vanilla 107cN_EmpiricalBeta_Copula.R $1 $2 > ./logs/"$1"_"$2"_107.out 2> ./logs/"$1"_"$2"_107.err
tail ./logs/"$1"_"$2"_107.out
tail ./logs/"$1"_"$2"_107.err

echo "117N running"
Rscript --vanilla 117cN_EmpiricalBeta_PoE.R $1 $2 > ./logs/"$1"_"$2"_117.out 2> ./logs/"$1"_"$2"_117.err
tail ./logs/"$1"_"$2"_117.out
tail ./logs/"$1"_"$2"_117.err

# echo "108N running"
# Rscript --vanilla 108N_Event_Splitting.R $1 $2 > ./logs/"$1"_"$2"_108N.out 2> ./logs/"$1"_"$2"_108N.err
# tail ./logs/"$1"_"$2"_108N.out
# tail ./logs/"$1"_"$2"_108N.err

echo "113N running"
python3 -u 113bN_proper_event_summary.py $1 $2 > ./logs/"$1"_"$2"_113.out 2> ./logs/"$1"_"$2"_113.err
tail ./logs/"$1"_"$2"_113.out
tail ./logs/"$1"_"$2"_113.err

echo "113BN running"
python3 -u 113bN_EC_proper_event_summary.py $1 $2 > ./logs/"$1"_"$2"_113B.out 2> ./logs/"$1"_"$2"_113B.err
tail ./logs/"$1"_"$2"_113B.out
tail ./logs/"$1"_"$2"_113B.err

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
