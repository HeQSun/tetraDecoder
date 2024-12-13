#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00
#SBATCH -J calculateSS

genotype_file=$1
win_size=$2
output_name=$3

python calculateSS.py $genotype_file $win_size $output_name





