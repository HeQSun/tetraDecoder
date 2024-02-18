# Here is the pipeline for aligning hic reads to clipped assembly: https://github.com/tangerzhang/ALLHiC/wiki
#      and separat hic reads of four haplotypes to each linkage group.

wd=/your/working/path/
cd ${wd}

# step 1. index reference

wd=/your/working/path/
cd ${wd}
#
ref_path="/your/assembly/path//a2_initial_assembly/"
#
for sample in A B C D E F G H I J O; do
    cd ${ref_path}/sample_${sample}/hifiasm_asm_v0p7/subset4_illu_re_align_purge_ovl/DNB_align/
    echo ${sample}": "
    ll clipped4_${sample}_hifiasm.p_utg.gfa.fa
    bsub -o index.log -q normal -R "rusage[mem=10000]" -R "span[hosts=1]" -M 10000 "bwa index -a bwtsw clipped4_${sample}_hifiasm.p_utg.gfa.fa"
    bsub -o index.log -q normal -R "rusage[mem=10000]" -R "span[hosts=1]" -M 10000 "samtools faidx clipped4_${sample}_hifiasm.p_utg.gfa.fa"
    # bsub -o index.log -q ioheavy -n 8 -R "rusage[mem=10000]" -R "span[hosts=1]" -M 10000 "bowtie2-build --threads 8 clipped4_${sample}_hifiasm.p_utg.gfa.fa clipped4_${sample}_hifiasm.p_utg.gfa.fa"
    cd ${wd}
done

# step 2. align hi-c reads to the ungrouped contigs

# 2.1 align each end:

wd=/your/working/path/bgi_deep_seq/
cd ${wd}
for sample in A B C D E F G H I J O; do
   #
   wd=/your/working/path/bgi_deep_seq/
   cd ${wd}
   mkdir hic_${sample}
   cd hic_${sample}
   mkdir bwa_align
   cd bwa_align
   genome=../../../../a2_initial_assembly/sample_${sample}/hifiasm_asm_v0p7/subset4_illu_re_align_purge_ovl/DNB_align/clipped4_${sample}_hifiasm.p_utg.gfa.fa
   readpath=//your/omic_read/path/read_BGI/sample_${sample}/
   Rs=`ls ${readpath}*_*.fq.gz`
   for readset in ${sample}_L1 ${sample}_L2 ${sample}_L3; do
       R1="${readpath}/${readset}_1.fq.gz"
       R2="${readpath}/${readset}_2.fq.gz"
       #
       ll ${genome} ${R1} ${R2}
       #
       bsub -m "hpc002 hpc003 hpc004 hpc005 hpc006 hpc007" -o bowtie2.log -e bowtie2.err -q multicore40 -n 30 -R "rusage[mem=36000]" -R "span[hosts=1]" -M 36000 "bowtie2 -x ${genome} -1 ${R1} -2 ${R2} -p 30 --local | samtools view -@ 30 -bS - | samtools sort -n -@ 30 -o ${readset}_PE_RN_sorted_local.bam -; samtools view -@ 30 -C -T ${genome} -o ${readset}_PE_RN_sorted_local.cram ${readset}_PE_RN_sorted_local.bam; ls -l *.bam *.cram"
       #
   done
   #
done

# Important preparation steps: before step 3 and followings, run s2_a4_ref_based_grouping_hifiasm_v0p7.sh first to get ref-groupped list of contigs - you have done this.

# step 3. extract lg-wise hic reads: samtools view -q 20 -o ${readset}_${lg}_extract_q20.bam ${readset}_${lg}_extract.bam - not used 20221126
#         omnic_read_extracter: includes intra-hap and inter-hap read pairs within the same lg
#
wd=/your/working/path/bgi_deep_seq/
cd ${wd}
for sample in A B C D E F G H I J O; do
   #
   wd=/your/working/path/bgi_deep_seq/
   cd ${wd}
   cd hic_${sample}
   cd bwa_align
   for readset in ${sample}_L1 ${sample}_L2 ${sample}_L3; do
       for lg in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12; do
           group=../../../../a4_alignment_based_linkage_grouping/zhifiasm_v0p7_dm_${sample}_asm/dm_res_grouping_details_${sample}_${lg}.txt
           ll ${group} ${readset}_PE_RN_sorted_local.cram
           bsub -o omnic_read_extracter_r.log -e omnic_read_extracter_r.err -q ioheavy -R "rusage[mem=3000]" -R "span[hosts=1]" -M 3000 " samtools view -F 3840 ${readset}_PE_RN_sorted_local.cram | omnic_read_extracter - ${group} ${readset}_${lg}; samtools sort -n -@ 1 -o ${readset}_${lg}_extract.bam ${readset}_${lg}_extract.sam; rm ${readset}_${lg}_extract.sam"
       done
       #
   done
   #
done
#

# step 4. extract paired-end hic reads in bam to fastqs
#
wd=/your/working/path/bgi_deep_seq/
cd ${wd}
for sample in A B C D E F G H I J O; do
   #
   wd=/your/working/path/bgi_deep_seq/
   cd ${wd}
   cd hic_${sample}
   cd bwa_align
   for readset in ${sample}_L1 ${sample}_L2 ${sample}_L3; do
       for lg in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12; do
           ls -l ${readset}_${lg}_extract.bam
           bsub -o samtools_bam2fastq_r.log -e samtools_bam2fastq_r.err -q normal -R "rusage[mem=3000]" -R "span[hosts=1]" -M 3000 "samtools fastq -1 ${readset}_${lg}_extract_R1.fastq.gz -2 ${readset}_${lg}_extract_R2.fastq.gz -0 ${readset}_${lg}_OTHER1.fastq.gz -s ${readset}_${lg}_OTHER2.fastq.gz -n -F 0x900 ${readset}_${lg}_extract.bam"
       done
       #
   done
   #
done

# now go to ./s2_a3_hic_alignment_hic_binning_new_method_for_phasing_asm_substep1.sh
