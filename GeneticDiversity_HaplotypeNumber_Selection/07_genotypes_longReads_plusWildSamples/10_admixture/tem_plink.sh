#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=3:00:00
#SBATCH -J Plink


chr=$1

VCF=/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/merged_vcf/chr/chr$chr"_"filterQUAL30.vcf.gz

FILE=chr$chr"_"filterQUAL30.LD0.3

#zcat $VCF | grep -v "##contig=<ID=scaffold" | awk -F"\t" 'BEGIN{ FS=OFS="\t" }{if($1 ~ /^chr/) $3="." ; print $0}' > tem_$FILE.vcf

plink --vcf tem_$FILE.vcf.gz --double-id --allow-extra-chr \
	--snps-only \
	--set-missing-var-ids @:# \
	--indep-pairwise 50 25 0.3 \
	--out $FILE --vcf-half-call 'missing' 
	
plink --vcf tem_$FILE.vcf.gz --double-id --allow-extra-chr \
	--set-missing-var-ids @:# \
	--snps-only \
	--extract $FILE.prune.in \
	--vcf-half-call 'missing' \
	--make-bed --out $FILE	

# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
cp $FILE.bim $FILE.bim_original
awk '{$1=0;print $0}' $FILE.bim > $FILE.bim.tmp
mv $FILE.bim.tmp $FILE.bim


