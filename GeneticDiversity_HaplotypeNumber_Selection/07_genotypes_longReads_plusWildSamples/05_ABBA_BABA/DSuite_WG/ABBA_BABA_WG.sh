#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=24:00:00

sample1=$1
sample2=$2
sample3=$3
sample4=$4
output_prefix=$5
folder=$6

mkdir -p $folder

vcf=/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/merged_vcf/wholegenome/WholeGenome_filterQUAL30.vcf.gz

bcftools query -l $vcf | awk '{print $1"\txxx"}' | \
    sed "s/${sample1}.*xx/${sample1}\t${sample1}/g" | \
    sed "s/${sample2}.*xx/${sample2}\t${sample2}/g" | \
    sed "s/${sample3}.*xx/${sample3}\t${sample3}/g" | \
    sed "s/${sample4}.*xx/${sample4}\tOutgroup/g" > species_sets_WholeGenome_${output_prefix}.txt

echo -e "P1\tP2\tP3\tDstatistic\tZ-score\tp-value\tf4-ratio\tBBAA\tABBA\tBABA" > Values_WholeGenome_$output_prefix.txt

# delete files from previous window
if test -f ./species_sets_WholeGenome_${output_prefix}_WholeGenome_$output_prefix"_tree.txt" ; then
  rm species_sets_WholeGenome_${output_prefix}_WholeGenome_$output_prefix"_"*.txt
fi

# calculate D values
/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/Dsuite/Build/Dsuite Dtrios -c -n WholeGenome_$output_prefix -t astral.merged_trees.WGenome.min25.nwk $vcf species_sets_WholeGenome_${output_prefix}.txt

# add values to output file:
grep "" species_sets_WholeGenome_${output_prefix}_WholeGenome_$output_prefix"_tree.txt" | grep -v "ABBA" >> Values_WholeGenome_$output_prefix.txt

rm species_sets_WholeGenome_${output_prefix}.txt

mv Values_WholeGenome_$output_prefix.txt $folder"/"




