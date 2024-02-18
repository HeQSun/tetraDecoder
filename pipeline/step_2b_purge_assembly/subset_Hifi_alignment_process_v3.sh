#!/bin/bash
# subset_Hifi_alignment_process_v3.sh (LSF-oriented)
# align HiFi reads to the (initial/purged) assembly to get a coverage view of contigs.
# 23 April 2022
# Hequan Sun

# usage

if [ $# -eq 0 ]
  then
echo  
echo "Usage: subset_Hifi_alignment_process_v3.sh <1.sample-id:A/B/..> <2.alignmen-folder> <3.subset-id:0/1/...> <4.thread> <5.genome> <6.reads>"
echo
exit 1
fi

echo 
echo
echo "##############################################################################################################" 
echo "subset_Hifi_alignment_process_v3 (LSF-oriented) starting ..."
current_date=`date`
echo ${current_date}

# step 0. get variables

echo "##############################################################################################################" 
echo "get sample id: $1" #  wd=$(pwd)
sample=$1
echo "get alignment folder: $2"
folder=$2
echo "get subset id: $3"
subset_id=$3
echo "get thread for read alignment: $4"
thread=$4
echo "get genome: $5"
genome=$5
echo "get reads: $6"
ptthifi=$6

# sample="B"
# folder=/netscratch/dep_mercier/grp_schneeberger/projects/Potato_multipleCultivars/s2_10_Cultivars_PacBio_HiFi/a2_initial_assembly/sample_${sample}/hifiasm_asm/XXX/
cd ${folder}

# genome=../${sample}_hifiasm.bp.p_utg.gfa.fasta
# genome=../${sample}_hifiasm.p_utg.gfa.fasta
# ptthifi=../sample_${sample}.fastq

ls -l ${genome} ${genome}.ctgsizes
ls -l ${ptthifi}

########################################################################################################################
echo "Step 1: read alignment and sort sam to bam: \"-N 0 allow 0 secondary alignment\" -- however, minimap2 does not this, but automatically recommending \"-N 1 --secondary=no\" " # -- 20200725
#
minimap2 -ax map-pb -t ${thread} -N 1 --secondary=no ${genome} ${ptthifi} | samtools view -@ ${thread} -bS - | samtools sort -@ ${thread} -o ptt${sample}_hiasm_ref_pilon_subset${subset_id}.bam -
current_date=`date`
echo "        read alignment done on ${current_date}." 
#
#
########################################################################################################################
echo "Step 2: find depth for each bam; and take the sum of depths from all bams at each position"
#
samtools depth -aa ptt${sample}_hiasm_ref_pilon_subset${subset_id}.bam > ptt${sample}_hiasm_ref_pilon_subset${subset_id}_pb_only_depth.txt
#
current_date=`date`
echo "Step 2 done on ${current_date}."
#
#
########################################################################################################################
echo "Step 3: find distribution of depth at genomic positions"
#
bsub -o samtoolsDepthHisto.log -e samtoolsDepthHisto.err -q normal -R "rusage[mem=2000]" -R "span[hosts=1]" -M 2000 "samtoolsDepthHisto ptt${sample}_hiasm_ref_pilon_subset${subset_id}_pb_only_depth.txt >samtoolsDepthHisto.log"
#
current_date=`date`
echo "Step 3 submitted to LSF, and you need to check samtoolsDepthHisto.log later. ${current_date}"
#
#
########################################################################################################################
echo "Step 4: find distribution of average depth at non-overlapping windows: winstep = winsize"
#
#chrsizes=../DNB_align/clipped_${sample}_hifiasm.bp.p_utg.gfa.fa.ctgsizes
#chrsizes=../DNB_align/clipped_${sample}_hifiasm.p_utg.gfa.fa.ctgsizes
chrsizes=${genome}.ctgsizes
samplecol=1
avgdepth=30
for winsize in 1000 10000 20000 50000; do
    bsub -o CNV_HQ_v3.log -e CNV_HQ_v3.err -q normal -R "rusage[mem=2000]" -R "span[hosts=1]" -M 2000 "CNV_HQ_v3 ${chrsizes} ptt${sample}_hiasm_ref_pilon_subset${subset_id}_pb_only_depth.txt ${winsize} ${winsize} ${samplecol} ${avgdepth} >slidingwin_depth_${winsize}.log"
done
#
current_date=`date`
echo "Step 4 submitted to LSF, and you need to check CNV_HQ_v3.log later. ${current_date}"
#
#
########################################################################################################################
echo "Step 5: Check how well it aligned => HiFi_alignment_quality_Hifiasm_ref_pilon_long_contigs.pdf --"
#           Result:
#                  => Only xxx Gb clipped off xxx reads; among which xxx reads made the major clipping of xxx Gb.
#                  => Alignments without proper cigar: xxx (xxx Gb)
#           Conclusion: Hifi reads are very well aligned to the contig-selected reference genome.
#
bsub -o read_clip_length_finder.log -e read_clip_length_finder.err -q normal -R "rusage[mem=2000]" -R "span[hosts=1]" -M 2000 "samtools view ptt${sample}_hiasm_ref_pilon_subset${subset_id}.bam | read_clip_length_finder 0 2 ptt${sample}_hiasm_ref_pilon_subset${subset_id}_clip > read_clip_length_finder_subset${subset_id}.log"
#
current_date=`date`
echo "Step 5 submitted to LSF, and you need to check read_clip_length_finder.log later. ${current_date}"
#
#
######################################################################################################################## bsub -q ioheavy -n 20 -R "span[hosts=1] rusage[mem=48000]" -M 48000
echo "Step 6: check stats of reads and index the bam for later IGVing"
# stats
bsub -o samtools_flagstat.log -e samtools_flagstat.err -q ioheavy -n 4 -R "rusage[mem=2000]" -R "span[hosts=1]" -M 2000 "samtools flagstat -@ 4 ptt${sample}_hiasm_ref_pilon_subset${subset_id}.bam"
# index
bsub -o samtools_index.log -e samtools_index.err -q ioheavy -n 4 -R "rusage[mem=2000]" -R "span[hosts=1]" -M 2000 "samtools index -@ 4 ptt${sample}_hiasm_ref_pilon_subset${subset_id}.bam"
#

current_date=`date`
echo "Step 6 submitted to LSF, and you need to check samtools_flagstat.log and samtools_index.log later. ${current_date}"
#
#
########################################################################################################################
# 
current_date=`date`
echo "subset_Hifi_alignment_process_v3 all done/submitted on ${current_date}."
echo 
echo

