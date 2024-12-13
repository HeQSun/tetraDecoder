#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --cpus-per-task=1
#SBATCH --mem=4700mb
#SBATCH --time=8:00:00
#SBATCH -J QUAL_Val

#run=03
run=$1
output=$2

echo "chr,QUAL,NVar" > $output
zcat merged_vcf/$run"/"ChrGroup_$run.vcf.gz | grep -v "^#" | awk '{print $1"\t"$6}' | sort | uniq -c | awk '{print $2","$3","$1}' >> $output



