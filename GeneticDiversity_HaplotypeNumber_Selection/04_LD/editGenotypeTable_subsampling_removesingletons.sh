#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=1-00:00:00
#SBATCH -J subSamplingGT

file=$1
output=$2
min_dis=$3

python3 editGenotypeTable_subsampling_removesingletons.py $file $output $min_dis 


