#!/bin/bash
#
# align illumina reads to the clipped unitig assembly with locally assembled contigs
#     output to: /a2_initial_assembly/sample_*/hifiasm_asm_v0p7/subset3_illu_re_align_purge_ovl/DNB_align/
# 09 SEP 2022
# Hequan Sun
#
#
#
#
echo "##############################################################################################################"
echo "get threads: $1"
nthread=$1
echo "get sample: $2"
sample=$2
echo "get subset id: $3"
subset_id=$3
#
subset_id_last=$((subset_id-1))
#
# subset_id=1
# sample=C

####COMMENT
#
# step 0. select sequences to build a new fasta
#
cd /your/working/path//
mkdir subset${subset_id}_illu_re_align_purge_ovl
cd subset${subset_id}_illu_re_align_purge_ovl
mkdir DNB_align
cd DNB_align
# build assembly with local assembly
awk 'FNR==1{print ""}1' ../../subset${subset_id_last}_illu_re_align_purge_ovl/DNB_align/clipped${subset_id_last}_${sample}_hifiasm.p_utg.gfa.fa ../../subset${subset_id}_local_asm/${sample}_assembly_local.fasta > clipped${subset_id}_${sample}_hifiasm.p_utg.gfa.fa
sed -i 's/contig_/l/g' clipped${subset_id}_${sample}_hifiasm.p_utg.gfa.fa
#
fasta_length clipped${subset_id}_${sample}_hifiasm.p_utg.gfa.fa | grep '>' | sed 's/>//g' > clipped${subset_id}_${sample}_hifiasm.p_utg.gfa.fa.ctgsizes
#
# calcuate N50
bsub -o calc_CN50.log -e calc_CN50.err "/path/to/GitHub-schneeberger-group-jac/toolbox/Analysis/Assembly/calc_CN50.pl clipped${subset_id}_${sample}_hifiasm.p_utg.gfa.fa 844000000 1 > clipped${subset_id}_${sample}_hifiasm.p_utg.gfa.N50.calc.result.txt"

# step 1. index reference

current_date=`date`
echo "        index reference ${current_date}."

#cd /your/working/path//subset${subset_id}_illu_re_align_purge_ovl/DNB_align

bsub -o fa_index.log -e fa_index.err "samtools faidx clipped${subset_id}_${sample}_hifiasm.p_utg.gfa.fa"
bowtie2-build --threads ${nthread} clipped${subset_id}_${sample}_hifiasm.p_utg.gfa.fa clipped${subset_id}_${sample}_hifiasm.p_utg.gfa.fa

current_date=`date`
echo "        index reference done on ${current_date}."
#### COMMENT


# step 2. align reads

current_date=`date`
echo "        align reads on ${current_date}."

#cd /your/working/path//subset${subset_id}_illu_re_align_purge_ovl/DNB_align
#
readpath=/your/read/path/to/clean_${sample}
ln -s ${readpath}/clean_${sample}_BGI_DNB_R1.fq.gz ln_${sample}_BGI_DNB_R1.fq.gz
ln -s ${readpath}/clean_${sample}_BGI_DNB_R2.fq.gz ln_${sample}_BGI_DNB_R2.fq.gz
R1=ln_${sample}_BGI_DNB_R1.fq.gz
R2=ln_${sample}_BGI_DNB_R2.fq.gz
#
bowtie2 -p ${nthread} -x clipped${subset_id}_${sample}_hifiasm.p_utg.gfa.fa -1 ${R1} -2 ${R2} | samtools view -@ ${nthread} -bS - | samtools sort -m 2G -@ ${nthread} -o ${sample}_PE_sorted.bam -
#
current_date=`date`
echo "        align reads done on ${current_date}."
#
samtools index -@ ${nthread} ${sample}_PE_sorted.bam
#
current_date=`date`
echo "        bam index done on ${current_date}."

#


########################################################################################################################
echo "Step 2: find depth for each bam"
#
samtools depth -@ ${nthread} -aa ${sample}_PE_sorted.bam > ${sample}_PE_sorted_only_depth.txt
#
current_date=`date`
echo "Step 2 done on ${current_date}."
#
#
########################################################################################################################
echo "Step 4: find distribution of average depth at non-overlapping windows: winstep = winsize"
#
chrsizes=clipped${subset_id}_${sample}_hifiasm.p_utg.gfa.fa.ctgsizes
samplecol=1
avgdepth=130
for winsize in 1000 10000 20000 50000; do
    bsub -o CNV_HQ_v3.log -e CNV_HQ_v3.err -q normal -R "rusage[mem=2000]" -R "span[hosts=1]" -M 2000 "CNV_HQ_v3 ${chrsizes} ${sample}_PE_sorted_only_depth.txt ${winsize} ${winsize} ${samplecol} ${avgdepth} >slidingwin_depth_${winsize}.log"
done
#
current_date=`date`
echo "Step 4 submitted to LSF, and you need to check CNV_HQ_v3.log later. ${current_date}"
#
########################################################################################################################
echo "Step 5: convert bam to cram"
nthread=8
bsub -o bam2cram.log -e bam2cram.err -q ioheavy -n ${nthread} -R "span[hosts=1] rusage[mem=2000]" -M 2000 "samtools view -@ ${nthread} --write-index -T clipped${subset_id}_${sample}_hifiasm.p_utg.gfa.fa -C -o ${sample}_PE_sorted.cram ${sample}_PE_sorted.bam"
#
# COMMENT
