#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=12:00:00
#SBATCH -J Pop_Par

list_samples=$1
chr_used=$2
win_size=$3
prefix=$4

# calculate population parameters:
python3 geneticDiversity_genotyping.py merged_genotypes.$prefix.$chr_used.txt  list_haplotypesNames.$prefix.$chr_used.txt $win_size $prefix"_"Win$win_size"_"$chr_used



