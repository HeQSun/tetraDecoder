#!/bin/bash

#usage
if [ "$#" -le 2 ]
  then
echo "blast_contigs_summary.sh <work_path> <lg> <real_hapi>"
exit 1
fi

work_path=$1
lg=$2
real_hapi=$3

cd ${work_path}

# collect hits
>all_LessThan350kb_blast_top1_stat.txt    
while read ctg; do
    cd ${work_path}
    cd ${ctg}_blast_result
#    cat ${ctg}_against_ncbi_nt_nt_*.oblast > ${ctg}_against_ncbi_nt_all.oblast
#    ls -l ${ctg}_against_ncbi_nt_all.oblast
#    ctgsize=../../../../homLG_${lg}_LG_${real_hapi}_LessThan350kb.ctgsizes
#    get_organelle_ctg_stat ${ctg}_against_ncbi_nt_all.oblast ${ctgsize} > ${ctg}_blast_top1_stat.txt            
    grep -A 1 'Report: ' ${ctg}_blast_top1_stat.txt >> ../all_LessThan350kb_blast_top1_stat.txt
    cd ..
done < ../../../homLG_${lg}_LG_${real_hapi}_LessThan350kb.ctgids # to confirm directory
# 
# get and plot organelle-hitting contigs
egrep 'chlor|mito|pla' all_LessThan350kb_blast_top1_stat.txt | sed 's/           //g' > all_LessThan350kb_blast_top1_stat_subset_organelle.txt
# work_path="./"
Rscript /biodata/dep_mercier/grp_schneeberger/projects/methods/src_shq/z_GameteBinning_potato_multivar/src_bash_scripts/viz_blast_ratio.R all_LessThan350kb_blast_top1_stat_subset_organelle.txt ${work_path} homLG_${lg}_LG_${real_hapi}
