#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00

chr=$1

mkdir -p splitVCF
mkdir -p splitVCF/$chr

ls -1 /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/splitVCF/$chr"/"*.vcf.gz > list_vcf_$chr.txt


for vcf in $(cat list_vcf_$chr.txt); do echo $win
vcf_file=${vcf##*/}
vcf_file_base=${vcf_file%.vcf.gz}
vcftools --gzvcf $vcf --recode --recode-INFO-all --keep list_WildSamples.txt  --out splitVCF/$chr"/"$vcf_file_base
gzip -c splitVCF/$chr"/"$vcf_file_base.recode.vcf > splitVCF/$chr"/"$vcf_file_base.vcf.gz
rm splitVCF/$chr"/"$vcf_file_base.recode.vcf
done



