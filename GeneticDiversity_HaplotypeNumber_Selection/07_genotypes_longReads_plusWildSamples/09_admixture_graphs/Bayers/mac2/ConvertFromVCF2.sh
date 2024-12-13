#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --cpus-per-task=1
#SBATCH --mem=4500mb
#SBATCH --time=2:00:00
#SBATCH -J VCFtoInput



# mamba activate env_others
# sbatch ../ConvertFromVCF.sh $vcf allele_count_chr$chr.txt SampleA_Hap_1 

vcf=$1
output=$2
sample_information=$3

python3 ConvertFromVCF.py $vcf $output $sample_information


