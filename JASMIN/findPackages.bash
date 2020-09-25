#!/bin/bash

#SBATCH -partition=short-serial
#SBATCH -o fp_%J.out
#SBATCH -e fp_%J.err
#SBATCH --time 00:30

echo "Running slimline_testing6"
module load jaspy/3.7/r20200606
Rscript --vanilla ~/AQUACAT/CodeABG/findPackages.R