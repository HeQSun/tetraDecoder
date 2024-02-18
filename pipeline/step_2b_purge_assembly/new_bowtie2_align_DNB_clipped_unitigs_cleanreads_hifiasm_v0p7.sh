#!/bin/bash
#    
# 1. purge overlaps among unitigs using purge_overlaps
# 2. align illumina reads to the clipped unitig assembly
#
# 11 AUG 2022
# Hequan Sun
#
#
# TASK TODO: align short reads to both the raw assembly and the clipped assembly, to check window-based coverage!
#
echo "##############################################################################################################" 
# set up parameters: 
subset_id=1
sample="XXXX"
gfa_file=${sample}.p_utg.gfa
wd="/this/is/your/working/path" # where gfa_file is
readpath="/this/is/the/path/to/your/wgs/short/reads/"

# step 0. clip overlap between unitigs: purge_overlaps
# 
cd ${wd}
mkdir subset${subset_id}_illu_re_align_purge_ovl
cd subset${subset_id}_illu_re_align_purge_ovl
mkdir DNB_align
cd DNB_align

# update name of read files below
ln -s ${readpath}/clean_${sample}_read_R1.fq.gz ln_${sample}_read_R1.fq.gz
ln -s ${readpath}/clean_${sample}_read_R2.fq.gz ln_${sample}_read_R2.fq.gz
#

#<<COMMENT
# step 0. clip unitigs => useful result: clipped_${sample}.p_utg.gfa.fa and clipped_${sample}.p_utg.gfa.fa.ctgsizes
min_ctg_len=15000
purge_overlaps ${gfa_file} ${min_ctg_len} > purge_overlaps.log
fasta_length clipped_${sample}.p_utg.gfa.fa | grep '>' | sed 's/>//g' > clipped_${sample}.p_utg.gfa.fa.ctgsizes
calc_CN50.pl clipped_${sample}.p_utg.gfa.fa 1000000000 1 > clipped_${sample}.p_utg.gfa.N50.calc.result.txt
#COMMENT

# step 1. index reference

current_date=`date`
echo "        index reference ${current_date}." 
bowtie2-build --threads ${nthread} clipped_${sample}.p_utg.gfa.fa clipped_${sample}.p_utg.gfa.fa
current_date=`date`
echo "        index reference done on ${current_date}." 

# step 2. align reads 

current_date=`date`
echo "        align reads on ${current_date}." 
#
R1=ln_${sample}_read_R1.fq.gz
R2=ln_${sample}_read_R2.fq.gz
#
bowtie2 -p ${nthread} -x clipped_${sample}.p_utg.gfa.fa -1 ${R1} -2 ${R2} | samtools view -@ ${nthread} -bS - | samtools sort -m 2G -@ ${nthread} -o ${sample}_PE_sorted.bam -
#
current_date=`date`
echo "        align reads done on ${current_date}."
#
samtools index -@ ${nthread} ${sample}_PE_sorted.bam
#
current_date=`date`
echo "        bam index done on ${current_date}." 
#
#
########################################################################################################################
echo "Step 2: find depth for each bam; and take the sum of depths from all bams at each position"
#
samtools depth -@ ${nthread} -aa ${sample}_PE_sorted.bam > ${sample}_PE_sorted_only_depth.txt
# note, you can delete this file, once CNV_HQ_v3 is successfully done: ${sample}_PE_sorted_only_depth.txt
#
current_date=`date`
echo "Step 2 done on ${current_date}."
#
#
########################################################################################################################
echo "Step 4: find distribution of average depth at non-overlapping windows: winstep = winsize"
#
chrsizes=clipped_${sample}.p_utg.gfa.fa.ctgsizes
samplecol=1
# Update depth value below: do you know roughly you the depth of wgs short reads?
avgdepth=100
# You can submit CNV_HQ_v3 with your job scheduling system and run them in parallel.
for winsize in 1000 10000 20000 50000; do
    CNV_HQ_v3 ${chrsizes} ${sample}_PE_sorted_only_depth.txt ${winsize} ${winsize} ${samplecol} ${avgdepth} > slidingwin_depth_${winsize}.log
done
#
current_date=`date`
echo "Step 4. to check CNV_HQ_v3.log. ${current_date}"






