#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=10:00:00
#SBATCH -J Plink


VCF=$1

FILE=WG_filterQUAL30.LD0.3
outFile=$2

#vcftools --gzvcf $VCF --remove-indels --recode --recode-INFO-all --stdout | gzip -c > tem_$FILE.vcf.gz

plink --vcf tem_$FILE.vcf.gz --double-id --allow-extra-chr \
        --snps-only \
        --geno 0.1 \
        --mac 2 \
        --set-missing-var-ids @:#:'$1':'$2' \
        --indep-pairwise 10kb 3000 0.5 \
        --out $outFile --vcf-half-call 'missing'

plink --vcf tem_$FILE.vcf.gz --double-id --allow-extra-chr \
        --set-missing-var-ids @:#:'$1':'$2' \
        --snps-only \
        --geno 0.1 \
        --mac 2 \
        --extract $outFile.prune.in \
        --vcf-half-call 'missing' \
        --make-bed --out $outFile

# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
cp $outFile.bim $outFile.bim_original
awk '{$1=0;print $0}' $outFile.bim > $outFile.bim.tmp
mv $outFile.bim.tmp $outFile.bim

