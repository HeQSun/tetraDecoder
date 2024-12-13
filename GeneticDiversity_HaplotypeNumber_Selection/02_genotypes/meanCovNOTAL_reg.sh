#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00
#SBATCH -J covNOTAL_reg

winSize=$1
list_samples=$2
prefix=$3
chr=$4

for haplotype_file in $(cat $list_samples )
do echo $haplotype_file

haplotype_name=$(echo $haplotype_file | sed 's@.*/01_nucmer_syri_vs_DM/nucmer_syri_@@g' | sed 's@vs_DM/refDM_chr.*_vs_qry_@@g' | sed 's@/syri.out@@g' )
echo $haplotype_name
echo $chr

awk '{if($11=="HDR" && $3-($2-1)>100) print $1"\t"$2-1"\t"$3}' $haplotype_file > $prefix"_"hdr_$chr"_"$haplotype_name.txt

awk '{if($11=="CPL" && $3-($2-1)>100) print $1"\t"$2-1"\t"$3}' $haplotype_file > $prefix"_"cpl_$chr"_"$haplotype_name.txt

awk '{if($1!="-" && $11=="NOTAL" && $3-($2-1)>100) print $1"\t"$2-1"\t"$3}' $haplotype_file > $prefix"_"NOTAL_$chr"_"$haplotype_name.txt

cat $prefix"_"NOTAL_$chr"_"$haplotype_name.txt $prefix"_"hdr_$chr"_"$haplotype_name.txt $prefix"_"cpl_$chr"_"$haplotype_name.txt  | sort -k1,1 -k2,2n > $prefix"_"ALLMissing_$chr"_"$haplotype_name.txt

bedtools merge -i $prefix"_"ALLMissing_$chr"_"$haplotype_name.txt > $prefix"_"$haplotype_name.bed

rm $prefix"_"hdr_$chr"_"$haplotype_name.txt $prefix"_"cpl_$chr"_"$haplotype_name.txt $prefix"_"NOTAL_$chr"_"$haplotype_name.txt $prefix"_"ALLMissing_$chr"_"$haplotype_name.txt

#awk '{if($1!="-" && $11=="NOTAL" && $3-($2-1)>100) print $1"\t"$2-1"\t"$3}' $haplotype_file | sort -k1,1 -k2,2n > $prefix"_"$haplotype_name.ded

done

cat $prefix"_"*hap*.bed | sort -k1,1 -k2,2n > $prefix"_All".bed

chr_size=$(awk 'BEGIN{maxVal=0} {if($2>$3) {line_max=$2} else {line_max=$3} ; if (line_max>maxVal) maxVal=line_max} END {print maxVal}' $prefix"_All".bed )

rm Win_$chr.bed
for start in $(seq 0 $winSize $chr_size )
do end=$(echo $start | awk -v WS=$winSize '{print $1+WS-1}' )
echo -e $chr'\t'$start'\t'$end >> Win_$chr.bed
done


bedtools coverage -a Win_$chr.bed -b $prefix"_All".bed -d | awk -v CHR=$chr 'BEGIN {win=0 ; SumCov=0 ; NLines=0} {if($2==win) {SumCov+=$5 ; NLines+=1} else {print CHR"\t"win"\t"SumCov/NLines ; win=$2 ; SumCov=$5 ; NLines=1} } END {print CHR"\t"win"\t"SumCov/NLines}' > meanCoverage_$prefix.txt


rm $prefix"_"*hap*.bed $prefix"_All".bed Win_$chr.bed


