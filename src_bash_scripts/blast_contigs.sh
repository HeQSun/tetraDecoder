#!/bin/bash
#Updated by Craig for most recent NCBI database in /opt/share/blastdb/ncbi/ (nt.000 - nt.109) (previously was nt.00 - nt.76)


#usage
if [ "$#" -le 4 ]
  then
echo "blast_contigs.sh <work_path> <ctg_path> <ctg_id> <lg> <real_hapi>"
exit 1
fi

work_path=$1
ctg_path=$2
ctg_id=$3
lg=$4
real_hapi=$5

cd ${work_path}

# 1. blast to get all hits for this ctg
for i in {0..9}; do 
    blastn -query ${ctg_path}/"chr"${ctg_id}.fa -db /opt/share/blastdb/ncbi/nt.00${i} -out ${ctg_id}against_nt_nt_0${i}.oblast -outfmt 0 -max_target_seqs 5
done
for i in {10..99}; do 
    blastn -query ${ctg_path}/"chr"${ctg_id}.fa -db /opt/share/blastdb/ncbi/nt.0${i} -out ${ctg_id}against_nt_nt_${i}.oblast -outfmt 0 -max_target_seqs 5
done
for i in {100..109}; do 
    blastn -query ${ctg_path}/"chr"${ctg_id}.fa -db /opt/share/blastdb/ncbi/nt.${i} -out ${ctg_id}against_nt_nt_${i}.oblast -outfmt 0 -max_target_seqs 5
done    
#
mkdir ${ctg_id}_blast_result
mv ${ctg_id}against_nt_nt_*.oblast ${ctg_id}_blast_result
#
# 2. collect to get top hits
#
cd ${ctg_id}_blast_result
cat ${ctg_id}against_nt_nt_*.oblast > ${ctg_id}against_nt_all.oblast
ls -l ${ctg_id}against_nt_all.oblast
#
ctgsize=../../../../homLG_${lg}_LG_${real_hapi}_LessThan350kb.ctgsizes
SCRIPT_DIR=/biodata/dep_mercier/grp_schneeberger/projects/methods/src_shq/Useful_pipeline/tetraploid_genome_assembly/z_aux_customized_tools_scripts
${SCRIPT_DIR}/get_organelle_ctg_stat ${ctg_id}against_nt_all.oblast ${ctgsize} > ${ctg_id}_blast_top1_stat.txt            
#
# 3. clean
rm ${ctg_id}*.oblast
#
# final result: ${ctg_id}_blast_top1_stat.txt to be further visualized

