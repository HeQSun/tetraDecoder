#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=12:00:00
#SBATCH -J HeAr

vcf=$1
output=$2
samples=$3

python ArtificialHeterozygosity.py $vcf $output $samples


