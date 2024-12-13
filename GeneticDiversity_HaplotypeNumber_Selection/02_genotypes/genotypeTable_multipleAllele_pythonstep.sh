#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=12:00:00
#SBATCH -J Genotype_table

list_samples=$1
chr_used=$2
prefix=$3


# calculate population parameters:
python3 genotypeTable_multipleAlleles.py merged_genotypes.$prefix.$chr_used.txt  list_haplotypesNames.$prefix.$chr_used.txt $prefix"_"$chr_used





