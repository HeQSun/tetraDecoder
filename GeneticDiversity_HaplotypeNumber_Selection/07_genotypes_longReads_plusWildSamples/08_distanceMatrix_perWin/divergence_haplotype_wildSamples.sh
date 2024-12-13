#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00

chr=$1
output_prefix=$2

mkdir -p $output_prefix
mkdir -p $output_prefix/chr$chr

ls -1 /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/splitVCF/$chr"/"*.vcf.gz > list_vcf_Chr$chr.txt

for vcf in $( cat list_vcf_Chr$chr.txt ); do echo $vcf
# Extract the filename
filename="${vcf##*/}"
# Extract the number using parameter expansion
number="${filename#split.}"
window="${number%%.*}"
python3 divergence_haplotype_wildSamples.py $vcf $chr $output_prefix/chr$chr"/"dis_chr$chr"_"Win$window
done



