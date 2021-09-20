#!/bin/bash

# This slurm sbatch script runs all the steps of the pipeline saving and reading
# to netCDF files.

#SBATCH -p short-serial
#SBATCH -o present_%J.out
#SBATCH -e present_%J.err
#SBATCH --time 23:59:00
#SBATCH --mem=12000
#SBATCH --chdir=/home/users/adagri/AQUACAT/CodeABG

echo "Running 102_113_netcdf_script"
module load jaspy

echo $EN
echo $PE

echo "101 running"
Rscript --vanilla 101_YAML_creation.R $EN $PE > ./logs/"$EN"_"$PE"_101.out 2> ./logs/"$EN"_"$PE"_101.err
tail ./logs/"$EN"_"$PE"_101.out
tail ./logs/"$EN"_"$PE"_101.err

echo "102 running"
Rscript --vanilla 102b_Threshold_Extract.R $EN $PE > ./logs/"$EN"_"$PE"_102.out 2> ./logs/"$EN"_"$PE"_102.err
tail ./logs/"$EN"_"$PE"_102.out
tail ./logs/"$EN"_"$PE"_102.err

echo "104 running"
Rscript --vanilla 104b_Event_Extract.R $EN $PE > ./logs/"$EN"_"$PE"_104.out 2> ./logs/"$EN"_"$PE"_104.err
tail ./logs/"$EN"_"$PE"_104.out
tail ./logs/"$EN"_"$PE"_104.err

echo "105N running"
Rscript --vanilla 105N_Event_Summary.R $EN $PE > ./logs/"$EN"_"$PE"_105N.out 2> ./logs/"$EN"_"$PE"_105N.err
tail ./logs/"$EN"_"$PE"_105N.out
tail ./logs/"$EN"_"$PE"_105N.err

echo "106N running"
Rscript --vanilla 106N_PoE_Estimation.R $EN $PE > ./logs/"$EN"_"$PE"_106N.out 2> ./logs/"$EN"_"$PE"_106N.err
tail ./logs/"$EN"_"$PE"_106N.out
tail ./logs/"$EN"_"$PE"_106N.err

if 

echo "107N running"
Rscript --vanilla 107N_EmpiricalBeta_Copula.R $EN $PE > ./logs/"$EN"_"$PE"_107N.out 2> ./logs/"$EN"_"$PE"_107N.err
tail ./logs/"$EN"_"$PE"_107N.out
tail ./logs/"$EN"_"$PE"_107N.err

echo "117N running"
Rscript --vanilla 117N_EmpiricalBeta_PoE.R $EN $PE > ./logs/"$EN"_"$PE"_117N.out 2> ./logs/"$EN"_"$PE"_117N.err
tail ./logs/"$EN"_"$PE"_117N.out
tail ./logs/"$EN"_"$PE"_117N.err

echo "108N running"
Rscript --vanilla 108N_Event_Splitting.R $EN $PE > ./logs/"$EN"_"$PE"_108N.out 2> ./logs/"$EN"_"$PE"_108N.err
tail ./logs/"$EN"_"$PE"_108N.out
tail ./logs/"$EN"_"$PE"_108N.err

echo "113N running"
Rscript --vanilla 113N_proper_event_summary.R $EN $PE> ./logs/"$EN"_"$PE"_113N.out 2> ./logs/"$EN"_"$PE"_113N.err
tail ./logs/"$EN"_"$PE"_113N.out
tail ./logs/"$EN"_"$PE"_113N.err

echo "113BN running"
Rscript --vanilla 113bN_EC_proper_event_summary.R $EN $PE> ./logs/"$EN"_"$PE"_113BN.out 2> ./logs/"$EN"_"$PE"_113BN.err
tail ./logs/"$EN"_"$PE"_113BN.out
tail ./logs/"$EN"_"$PE"_113BN.err


echo "Complete"
