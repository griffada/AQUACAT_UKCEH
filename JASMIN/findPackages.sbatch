#!/bin/bash

#SBATCH -p short-serial
#SBATCH -o speed_%J.out
#SBATCH -e speed_%J.err
#SBATCH --time 00:20
#SBATCH --mem=8G
#SBATCH --chdir=/home/users/adagri/AQUACAT/CodeABG

echo "Running slimline_testing6"
module load jaspy/3.7/r20200606
Rscript --vanilla findPackages.R