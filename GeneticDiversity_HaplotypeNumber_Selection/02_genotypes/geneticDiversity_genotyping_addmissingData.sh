#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00
#SBATCH -J addMissingGT

#python geneticDiversity_genotyping_missingData.sh AllSamplesPopPar_Win10000_chr01_genotypeTable.txt list_haplotypes_chr01.txt AllSamplesPopPar_Win10000_chr01_missing chr01
#genotype_table=AllSamplesPopPar_Win10000_chr01_genotypeTable.txt
#list_samples=list_haplotypes_chr01.txt
#prefix=test_missing_AllSamplesPopPar_Win10000_chr01
#chr=chr01

genotype_table=$1
list_samples=$2
prefix=$3
chr=$4

head -1 $genotype_table | sed 's/\t/\n/'g | grep hap | cat -n | sed 's/ //g' > sample_order.$prefix.txt
cp $genotype_table $prefix"_"genotypeTable.txt
awk '{if($1!="CHR") print $1"\t"$2"\t"$3}' $genotype_table > $prefix"_"genotype.bed 

for haplotype_file in $(cat $list_samples )
do echo $haplotype_file

haplotype_name=$(echo $haplotype_file | sed 's@.*/01_nucmer_syri_vs_DM/nucmer_syri_@@g' | sed 's@vs_DM/refDM_chr.*_vs_qry_@@g' | sed 's@/syri.out@@g' )
echo $haplotype_name
echo $chr
index=$(grep -w $haplotype_name sample_order.$prefix.txt | awk '{print $1}')
echo $index

awk '{if($11=="DEL" && $3-($2-1)>100) print $1"\t"$2-1"\t"$3}' $haplotype_file > $prefix"_"del_$chr"_"$haplotype_name.txt

awk '{if($11=="HDR" && $3-($2-1)>100) print $1"\t"$2-1"\t"$3}' $haplotype_file > $prefix"_"hdr_$chr"_"$haplotype_name.txt

awk '{if($11=="CPL" && $3-($2-1)>100) print $1"\t"$2-1"\t"$3}' $haplotype_file > $prefix"_"cpl_$chr"_"$haplotype_name.txt

awk '{if($1!="-" && $11=="NOTAL" && $3-($2-1)>100) print $1"\t"$2-1"\t"$3}' $haplotype_file > $prefix"_"NOTAL_$chr"_"$haplotype_name.txt

cat $prefix"_"del_$chr"_"$haplotype_name.txt $prefix"_"NOTAL_$chr"_"$haplotype_name.txt $prefix"_"hdr_$chr"_"$haplotype_name.txt $prefix"_"cpl_$chr"_"$haplotype_name.txt  | sort -k1,1 -k2,2n > $prefix"_"ALLMissing_$chr"_"$haplotype_name.txt

bedtools merge -i $prefix"_"ALLMissing_$chr"_"$haplotype_name.txt > $prefix"_"ALLMissing_merged_$chr"_"$haplotype_name.txt

#bedtools intersect  -wa -wb -a $prefix"_"genotype.bed -b $prefix"_"del_$chr"_"$haplotype_name.txt | awk '{varSize=$3-$2 ; delSize=$6-$5 ; if (varSize!=delSize) print $2"\t"$3}' > list_$prefix"_"del_$chr"_"$haplotype_name.txt
#del_file=list_$prefix"_"del_$chr"_"$haplotype_name.txt

#python geneticDiversity_genotyping_addmissingData.py $prefix"_"genotypeTable.txt $index $del_file $haplotype_name

bedtools intersect  -wa -wb -a $prefix"_"genotype.bed -b $prefix"_"ALLMissing_merged_$chr"_"$haplotype_name.txt | awk '{varSize=$3-$2 ; delSize=$6-$5 ; if (varSize!=delSize) print $2"\t"$3}' > list_$prefix"_"ALLMissing_merged_$chr"_"$haplotype_name.txt
ALLMissing_file=list_$prefix"_"ALLMissing_merged_$chr"_"$haplotype_name.txt

python geneticDiversity_genotyping_addmissingData.py $prefix"_"genotypeTable.txt $index $ALLMissing_file $haplotype_name

mv tem_$prefix"_"genotypeTable.txt $prefix"_"genotypeTable.txt

#rm $prefix"_"del_$chr"_"$haplotype_name.txt list_$prefix"_"del_$chr"_"$haplotype_name.txt

rm $prefix"_"del_$chr"_"$haplotype_name.txt $prefix"_"NOTAL_$chr"_"$haplotype_name.txt $prefix"_"ALLMissing_$chr"_"$haplotype_name.txt $prefix"_"ALLMissing_merged_$chr"_"$haplotype_name.txt list_$prefix"_"ALLMissing_merged_$chr"_"$haplotype_name.txt $prefix"_"hdr_$chr"_"$haplotype_name.txt $prefix"_"cpl_$chr"_"$haplotype_name.txt

done

rm sample_order.$prefix.txt $prefix"_"genotype.bed 


