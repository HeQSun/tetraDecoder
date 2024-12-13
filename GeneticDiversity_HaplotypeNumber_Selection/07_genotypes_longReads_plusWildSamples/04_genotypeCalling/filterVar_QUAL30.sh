#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=4700mb
#SBATCH --time=4:00:00
#SBATCH -J QUAL_Val

#run=03
run=$1
chr=$2

mkdir -p merged_vcf/chr

bcftools view -O z --targets $chr -o merged_vcf/chr"/"$chr"_"filterQUAL30.vcf.gz -e 'QUAL<=30' merged_vcf/$run"/"ChrGroup_$run.vcf.gz

tabix -p vcf merged_vcf/chr"/"$chr"_"filterQUAL30.vcf.gz


