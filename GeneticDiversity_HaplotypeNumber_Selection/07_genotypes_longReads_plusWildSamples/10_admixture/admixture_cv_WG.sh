#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=2-00:00:00
#SBATCH -J admixture

# this environment contains plink, admixture
# conda activate env_others

k_value=$1
FILE=$2

admixture --cv $FILE.bed $k_value > log${k_value}.$FILE.out


