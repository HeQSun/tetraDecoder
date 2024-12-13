#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00


#chr=12
#sample1=SampleE_Hap_4
#sample2=SRR15458993
#sample3=SRR15458986
#sample4=SRR15458989
#output_prefix=bre_vem_Cul_cho_E4

chr=$1
sample1=$2
sample2=$3
sample3=$4
sample4=$5
output_prefix=$6
folder=$7

mkdir -p $folder

ls -1 /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/splitVCF/$chr"/"*.vcf.gz > list_vcf_Chr$chr"_"$output_prefix.txt

vcf=$(head -1 list_vcf_Chr$chr"_"$output_prefix.txt)

bcftools query -l $vcf | awk '{print $1"\txxx"}' | \
    sed "s/${sample1}.*xx/${sample1}\t${sample1}/g" | \
    sed "s/${sample2}.*xx/${sample2}\t${sample2}/g" | \
    sed "s/${sample3}.*xx/${sample3}\t${sample3}/g" | \
    sed "s/${sample4}.*xx/${sample4}\tOutgroup/g" > species_sets_${chr}_${output_prefix}.txt

# delete pre existing output file:
#if test -f ./Values_Chr$chr"_"$output_prefix.txt ; then
#  rm Values_Chr$chr"_"$output_prefix.txt
#fi

echo -e "chr\twin\tP1\tP2\tP3\tDstatistic\tZ-score\tp-value\tf4-ratio\tBBAA\tABBA\tBABA" > Values_Chr$chr"_"$output_prefix.txt


for vcf in $( cat list_vcf_Chr$chr"_"$output_prefix.txt ); do echo $vcf
# Extract the filename
filename="${vcf##*/}"
# Extract the number using parameter expansion
number="${filename#split.}"
window="${number%%.*}"
# delete files from previous window
if test -f ./species_sets_${chr}_${output_prefix}_chr$chr"_"$output_prefix"_tree.txt" ; then
  rm species_sets_${chr}_${output_prefix}_chr$chr"_"$output_prefix"_"*.txt
fi
# calculate D values
/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/Dsuite/Build/Dsuite Dtrios -c -n chr$chr"_"$output_prefix -t astral.merged_trees.WGenome.min25.nwk $vcf species_sets_${chr}_${output_prefix}.txt
# add values to output file:
grep "" species_sets_${chr}_${output_prefix}_chr$chr"_"$output_prefix"_tree.txt" | grep -v "ABBA" | awk -v win="$window" -v chrN="$chr" 'BEGIN {OFS="\t"} {print chrN, win, $0}' >> Values_Chr$chr"_"$output_prefix.txt
rm species_sets_${chr}_${output_prefix}_chr$chr"_"$output_prefix"_"*.txt
done

rm species_sets_${chr}_${output_prefix}.txt list_vcf_Chr$chr"_"$output_prefix.txt

mv Values_Chr$chr"_"$output_prefix.txt $folder"/"


