#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=3:00:00
#SBATCH -J Plink


VCF=$1

FILE=WG_filterQUAL30.LD0.3

#vcftools --gzvcf $VCF --remove-indels --recode --recode-INFO-all --stdout | gzip -c > tem_$FILE.vcf.gz

plink --vcf tem_$FILE.vcf.gz --double-id --allow-extra-chr \
        --snps-only \
        --set-missing-var-ids @:#:'$1':'$2' \
        --indep-pairwise 50 25 0.3 \
        --out $FILE --vcf-half-call 'missing'

plink --vcf tem_$FILE.vcf.gz --double-id --allow-extra-chr \
        --set-missing-var-ids @:#:'$1':'$2' \
        --snps-only \
        --extract $FILE.prune.in \
        --vcf-half-call 'missing' \
        --make-bed --out $FILE

# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
cp $FILE.bim $FILE.bim_original
awk '{$1=0;print $0}' $FILE.bim > $FILE.bim.tmp
mv $FILE.bim.tmp $FILE.bim

