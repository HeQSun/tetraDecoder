#!/bin/bash

#usage
if [ "$#" -le 7 ]
  then
echo "hic_binning.sh <1.bam> <2.min_hic_contact> <3.dm_res_grouping_details> <4.tig_win_marker> <5.min_hapctg_size> <6.allele_table> <7.ctg_alignment_to_dm> <8.lg>"
exit 1
fi

bam=$1
min_hic_contact=$2
dm_res_grouping_details=$3
tig_win_marker=$4
min_hapctg_size=$5
allele_table=$6
ctg_alignment_to_dm=$7
lg=$8

samtools view -F 3840 ${bam} | hic_binning_vt - ${min_hic_contact} ${dm_res_grouping_details} ${tig_win_marker} ${min_hapctg_size} ${allele_table} ${ctg_alignment_to_dm} ${lg} > hic_binning.log
