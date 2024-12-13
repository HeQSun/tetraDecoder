#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=12:00:00
#SBATCH -J Hap_pi

input=$1
output=$2
chr=$3

# calculate population parameters:
python3 calculatePi_haplotypes.py $input $output $chr

