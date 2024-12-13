#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=3:00:00

# mamba activate env_others
# sbatch ../ConvertFromVCF.sh $vcf allele_count_chr$chr.txt SampleA_Hap_1 

vcf=$1
output=$2
sample_information=$3

python3 ConvertFromVCF.py $vcf $output $sample_information


