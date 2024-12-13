#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00
#SBATCH -J meanLD_cal

chrID=$1
win_size=$2
chr_size=$3
output_name=$4

python3 meanLDvalues_subsampling.py $chrID $win_size $chr_size $output_name


