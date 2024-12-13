#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=1-00:00:00
#SBATCH -J admixture

# this environment contains plink, admixture
# conda activate env_others

chr=$1
k_value=$2
FILE=chr$chr"_"filterQUAL30.LD0.3

admixture --cv $FILE.bed $k_value > log${k_value}.$FILE.out


