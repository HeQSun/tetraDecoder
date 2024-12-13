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

# python3 /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr01_missing_genotypeTable.txt 10000 NSS_chr01
python calculateSS_removingNAN.py $genotype_file $win_size $output_name


