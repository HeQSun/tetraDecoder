#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=8:00:00
#SBATCH -J splitvcf

vcf=$1
output=$2
winSize=$3
nSnps=$4

mkdir -p splitVCF
mkdir -p splitVCF/$output

zcat $vcf | grep -v "##contig=<ID=scaffold" | awk -F"\t" 'BEGIN{ FS=OFS="\t" }{if($1 ~ /^chr/) $3="." ; print $0}' > splitVCF/tem_$output.vcf

java -jar /dss/dsshome1/lxc03/di36guz2/miniconda3/envs/myjvarkit/share/jvarkit-2024.04.20-0/jvarkit.jar biostar497922 -n $nSnps -o splitVCF/$output -D $winSize ./splitVCF/tem_$output.vcf

rm ./splitVCF/tem_$output.vcf


