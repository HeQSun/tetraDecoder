#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=12:00:00
#SBATCH -J Pop_Par

list_samples=$1
chr_used=$2
win_size=$3
prefix=$4

#ls -1 ../01_nucmer_syri_vs_DM/nucmer_syri_*_vs_DM/refDM_$chr_used""*/syri.out > list_haplotypes_$chr_used.txt

# produce genotype table per haplotype/sample
for haplotype_file in $(cat $list_samples )
do echo $haplotype_file
haplotype_name=$(echo $haplotype_file | sed 's@.*/01_nucmer_syri_vs_DM/nucmer_syri_@@g' | sed 's@vs_DM/refDM_chr.*_vs_qry_@@g' | sed 's@/syri.out@@g' )
echo $haplotype_name
echo $chr_used
#produce bed file per haplotype with chr, start-1, end, ref, alt, sample
awk -v del="DEL" -v ins="INS" -v snp="SNP" -v haplotype=$haplotype_name '{if($11==del || $11==ins || $11==snp) print $1"\t"$2-1"\t"$3"\t"$4"\t"$5"\t"haplotype}' $haplotype_file | sort | uniq > genotypes.$prefix.$chr_used.$haplotype_name.txt
# produce file with list og haplotype's names:
echo $haplotype_name >> list_haplotypesNames.$prefix.$chr_used.txt
done

# merged file with genotypes with all samples:
cat genotypes.$prefix.$chr_used.* | sort -nk2,2 -k4,4 -k5,5 | awk -v chr=$chr_used 'BEGIN {sample_list=""; startpos=0; endpos=0; ref=""; alt=""} {if($2==startpos && $3==endpos && $4==ref && $5==alt) sample_list=sample_list","$6 ; else {{if(sample_list!="") print chr"\t"startpos"\t"endpos"\t"ref"\t"alt"\t"sample_list}; sample_list=$6; startpos=$2; endpos=$3; ref=$4; alt=$5 }} END {print chr"\t"startpos"\t"endpos"\t"ref"\t"alt"\t"sample_list}' > merged_genotypes.$prefix.$chr_used.txt

# remove files:
rm genotypes.$prefix.$chr_used.*

# calculate population parameters:
python3 geneticDiversity_genotyping.py merged_genotypes.$prefix.$chr_used.txt  list_haplotypesNames.$prefix.$chr_used.txt $win_size $prefix"_"Win$win_size"_"$chr_used



