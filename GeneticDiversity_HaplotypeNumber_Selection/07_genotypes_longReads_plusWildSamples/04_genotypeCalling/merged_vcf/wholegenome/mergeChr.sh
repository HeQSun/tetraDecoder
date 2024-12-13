#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00
#SBATCH -J MergeVCF

vcf=$(ls ../chr/*vcf.gz)
bcftools concat --threads 5 $vcf | grep -v "##contig=<ID=scaffold" | awk -F"\t" 'BEGIN{ FS=OFS="\t" }{if($1 ~ /^chr/) $3="." ; print $0}' | gzip -c > WholeGenome_filterQUAL30.vcf.gz


