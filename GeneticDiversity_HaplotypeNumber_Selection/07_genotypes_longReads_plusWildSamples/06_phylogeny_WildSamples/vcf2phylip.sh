#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00
#SBATCH -J vcf2phylip


chr=$1

mkdir -p splitVCF_phy
mkdir -p splitVCF_phy/$chr

ls -1 /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/splitVCF/$chr"/"*.vcf.gz > list_vcf_$chr.txt

for vcf in $(cat list_vcf_$chr.txt); do echo $win
vcf_file=${vcf##*/}
vcf_file_base=${vcf_file%.vcf.gz}
python /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/vcf2phylip/vcf2phylip.py --input ./splitVCF/$chr"/"$vcf_file --fasta --nexus --nexus-binary --min-samples-locus 10 --output-folder ./splitVCF_phy/$chr 
done


"done finished"


