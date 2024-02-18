# Here for locally assembled contigs from purge3, we remove their ends with low coverage, and define tig-markers!

# step 1. create fasta with sequence selection
# Note: cannot use multicore20/40 due to TERM_THREADLIMIT: job killed after reaching LSF thread limit.
# 1.1 align short reads

wd=/your/working/path//
cd ${wd}
# TODO J ioheavy - vs done A B C D E F G H
for sample in E I; do
#
wd=/your/working/path//
cd ${wd}/sample_${sample}
cd hifiasm_asm_v0p7
#
subset_id=4
mkdir subset${subset_id}_illu_re_align_purge_ovl
cd subset${subset_id}_illu_re_align_purge_ovl
mkdir DNB_align
cd DNB_align
#
nthread=40
bsub -q ioheavy -n ${nthread} -R "span[hosts=1] rusage[mem=100000]" -M 100000 -o subset_DNB_alignment_process.log -e subset_DNB_alignment_process.err "../../../../new_bowtie2_align_DNB_clipped_unitigs_cleanreads_4purge${subset_id}_hifiasm_v0p7.sh ${nthread} ${sample} ${subset_id} > new_bowtie2_align_DNB_clipped_unitigs_cleanreads_4purge${subset_id}_hifiasm_v0p7.log"
#
cd ${wd}
#
done

# 1.2 get new assembly stats
#
wd=/your/working/path//
cd ${wd}
head -n 100 sample_*/hifiasm_asm_v0p7/subset4_illu_re_align_purge_ovl/DNB_align/*.N50.calc.result.txt > summary_hifiasm_v0p7_N50_purged_4th.txt
#head -n 100 sample_*/hifiasm_asm_v0p7/*.N50.calc.result.txt                                           > summary_hifiasm_v0p7_N50_raw.txt

# 1.3 R visualize
#
cd /your/working/path//
sample_x_cnv_initial_contigs_winsize10kb_step10kb_depth_DNB_purged_unitig_level_hifiasm_v0p7.R
sed -e '$s/$/\n/' -s  sample_*/hifiasm_asm_v0p7/subset4_illu_re_align_purge_ovl/DNB_align/sample_*_tig_ratio > sample_10_tig_ratio_summary_4th_purge.txt


# 2. align hifi reads

# 2.1 align long reads: get bam and convert bam to cram
for sample in A B C D E F G H I J O; do
#
subset_id=4
wd=/your/working/path//
cd ${wd}/sample_${sample}
cd hifiasm_asm_v0p7
cd subset${subset_id}_illu_re_align_purge_ovl
#
mkdir HIFI_align
cd HIFI_align
#
this_path="./"
nthread=20
genome=../DNB_align/clipped${subset_id}_${sample}_hifiasm.p_utg.gfa.fa
reads=../../../../../a1_HiFi_stats/sample_${sample}/sample_${sample}_zcat_HiFi_trim.fastq.gz
ll ${genome} ${reads}
# run substep1. get alignments etc.
bsub -q ioheavy -n ${nthread} -R "span[hosts=1] rusage[mem=80000]" -M 80000 -o subset_Hifi_alignment_process_v3.log -e subset_Hifi_alignment_process_v3.err "../../../../subset_Hifi_alignment_process_v3.sh ${sample} ${this_path} ${subset_id} ${nthread} ${genome} ${reads} > subset_Hifi_alignment_process_details.log"
# run substep2. convert bam to cram after finishing alignment above: save ~45%
#
cd ${wd}/
#
done


# step 3: call variants with hifi reads:
#
# 3.1 split hifi bam
wd=/your/working/path//
cd ${wd}
#
for sample in A B C D E F G H I J; do
#
wd=/your/working/path//
cd ${wd}/sample_${sample}
cd hifiasm_asm_v0p7
#
subset_id=4
cd subset${subset_id}_illu_re_align_purge_ovl
cd HIFI_align
#
mkdir variant_call_sc
cd variant_call_sc
#
genome=../../DNB_align/clipped${subset_id}_${sample}_hifiasm.p_utg.gfa.fa
bam=../ptt${sample}_hiasm_ref_pilon_subset${subset_id}.bam
ll ${genome} ${bam}
#
sed 's/\t/\t1\t/g' ../../DNB_align/clipped${subset_id}_${sample}_hifiasm.p_utg.gfa.fa.ctgsizes | sort -k3,3 --random-sort > clipped${subset_id}_${sample}_interval
split -l 800 clipped${subset_id}_${sample}_interval clipped${subset_id}_${sample}_interval_split
for interval in clipped${subset_id}_${sample}_interval_split*; do mv ${interval} ${interval}.txt; done
#
for interval in clipped*_*_interval_split*.txt; do
    this_base=$(basename "$interval")
    this_base_no_suffix=${this_base%.*}
    cat ${interval} > ${this_base_no_suffix}.bed
    bsub -q normal -R "span[hosts=1] rusage[mem=1000]" -M 1000 -o split_bam.log -e split_bam.err "samtools view -h -b -L ${this_base_no_suffix}.bed ${bam} > ${this_base_no_suffix}.bam; samtools index ${this_base_no_suffix}.bam"
done
#
wd=/your/working/path//
cd ${wd}
#
done
#
# 3.2 - perform hifi variant calling: singularity version 3.6.4 plus deepvariant v1.4.0
#         otherwise error: Set an empty environment as LD_LIBRARY_PATH may mix dependencies, just rely only on the library cache or its own lookup mechanism.
#                see here: https://github.com/apptainer/singularity/pull/5669, export LD_LIBRARY_PATH=""
#         caution on settings for running "singularity run deepvariant"
#
wd=/your/working/path//
cd ${wd}
#
for sample in A B C D E F G H I J; do
#
wd=/your/working/path//
cd ${wd}/sample_${sample}
cd hifiasm_asm_v0p7
#
subset_id=4
cd subset${subset_id}_illu_re_align_purge_ovl
cd HIFI_align
#
cd variant_call_sc
# redundant - remove later after finishing deepvariant
genome=../../DNB_align/clipped${subset_id}_${sample}_hifiasm.p_utg.gfa.fa
genome_index=../../DNB_align/clipped${subset_id}_${sample}_hifiasm.p_utg.gfa.fa.fai
cp -a ${genome} ${genome_index} ./
ll *.fa *.fai
#
for interval in clipped*_*_interval_split*.txt; do
    this_base=$(basename "$interval")
    this_base_no_suffix=${this_base%.*}
    #
    # export TMPDIR="$PWD/zdv1_${this_base_no_suffix}_tmp/";export LC_CTYPE="";export LANG="";
    this_flag="zdv1"
    bsub -q ioheavy -n 4 -R "rusage[mem=16000]" -M 16000 -o ${this_flag}_${this_base_no_suffix}.log -e ${this_flag}_${this_base_no_suffix}.err "mkdir ${this_flag}_${this_base_no_suffix}_tmp; export TMPDIR=\"$PWD/${this_flag}_${this_base_no_suffix}_tmp/\";export LC_CTYPE=\"\";export LANG=\"\"; export LD_LIBRARY_PATH=\"\"; singularity run --bind ${PWD} /srv/netscratch/dep_mercier/grp_schneeberger/software/deepvariant/deepvariant.img /opt/deepvariant/bin/run_deepvariant --model_type PACBIO --ref clipped${subset_id}_${sample}_hifiasm.p_utg.gfa.fa --reads ${this_base_no_suffix}.bam --intermediate_results_dir ./${this_flag}_${this_base_no_suffix}_tmp --output_vcf ${this_flag}_${this_base_no_suffix}.vcf.gz --num_shards 4 "
    #
done
#
wd=/your/working/path//
cd ${wd}
#
done
#
#
# 3.3 merge hifi sub-vcfs for each sample and find raw hifi SNPs
#     result of raw hifi SNPs
# 2nd call
#	626386 z_merged_clipped2_A_interval_hifi.vcf
#	485958 z_merged_clipped2_B_interval_hifi.vcf
#	574118 z_merged_clipped2_C_interval_hifi.vcf
#	686542 z_merged_clipped2_D_interval_hifi.vcf
#	648359 z_merged_clipped2_E_interval_hifi.vcf
#	561776 z_merged_clipped2_F_interval_hifi.vcf
#	610965 z_merged_clipped2_G_interval_hifi.vcf
#	474294 z_merged_clipped2_H_interval_hifi.vcf
#	626522 z_merged_clipped2_I_interval_hifi.vcf
#	749426 z_merged_clipped2_J_interval_hifi.vcf
# 4th call - high quality snps
#	13300 z_merged_clipped4_A_interval_hifi_HQ_SNPs.vcf
#	 8682 z_merged_clipped4_B_interval_hifi_HQ_SNPs.vcf
#	 8078 z_merged_clipped4_C_interval_hifi_HQ_SNPs.vcf
#	 9580 z_merged_clipped4_D_interval_hifi_HQ_SNPs.vcf
#	 9943 z_merged_clipped4_E_interval_hifi_HQ_SNPs.vcf
#	 8817 z_merged_clipped4_F_interval_hifi_HQ_SNPs.vcf
#	12963 z_merged_clipped4_G_interval_hifi_HQ_SNPs.vcf
#	 6410 z_merged_clipped4_H_interval_hifi_HQ_SNPs.vcf
#	10154 z_merged_clipped4_I_interval_hifi_HQ_SNPs.vcf
#	13819 z_merged_clipped4_J_interval_hifi_HQ_SNPs.vcf
#
#
wd=/your/working/path//
cd ${wd}
#
for sample in A B C D E F G H I J; do
#
wd=/your/working/path//
cd ${wd}/sample_${sample}
cd hifiasm_asm_v0p7
#
subset_id=4
cd subset${subset_id}_illu_re_align_purge_ovl
cd HIFI_align
cd variant_call_sc
#
cat zdv1_clipped${subset_id}_*_interval_split*.vcf.gz | gunzip | grep -E -v '#|PASS' | awk '$5=="A" || $5=="C" || $5=="G" ||$5=="T" ' | awk '$4=="A" || $4=="C" || $4=="G" ||$4=="T" ' > z_merged_clipped${subset_id}_${sample}_interval_hifi.vcf
sed 's/:/\t/g' z_merged_clipped${subset_id}_${sample}_interval_hifi.vcf | awk '$16>15 && $17>=35 && $19>=0.38 && $19<=0.62' > z_merged_clipped${subset_id}_${sample}_interval_hifi_HQ_SNPs.vcf
wc -l z_merged_clipped${subset_id}_${sample}_interval_hifi.vcf z_merged_clipped${subset_id}_${sample}_interval_hifi_HQ_SNPs.vcf
#
wd=/your/working/path//
cd ${wd}
#
done

# step 4. define tig markers
#
#         new window marker generation => cnv_winsize10000_step10000_hq_merged_vs_hifi_markers_20221018_wsize50kb_final.txt considers hifi depth for marker type correction!
#             note: integrate correction with hifi depth to last column of hifi+illu:cnv_winsize10000_step10000_hq.txt
#
wd=/your/working/path//
cd ${wd}
#
for sample in A B C D E F G H I J O; do
#
wd=/your/working/path//
cd ${wd}/sample_${sample}
cd hifiasm_asm_v0p7
#
subset_id=4
cd subset${subset_id}_illu_re_align_purge_ovl
#
#### paste ./DNB_align/cnv_winsize10000_step10000_hq.txt ./HIFI_align/cnv_winsize10000_step10000_hq.txt | cut -f1-7,12,13 > ${sample}_cnv_winsize10000_step10000_hq_merged_vs_hifi.txt
# new version of tig_marker_finder considers hifi depth as well! 20200909
illu_median_hap_cov=$(cat ./DNB_align/sample_${sample}_avg_cov)
hifi_median_hap_cov=$(cat ./HIFI_align/sample_${sample}_avg_cov)
#### echo "tig_marker_finder_v2 ${sample}_cnv_winsize10000_step10000_hq_merged_vs_hifi.txt ${illu_median_hap_cov} ${hifi_median_hap_cov} 500000 20221102 > tig_marker_finder.log"
#### tig_marker_finder_v2 ${sample}_cnv_winsize10000_step10000_hq_merged_vs_hifi.txt ${illu_median_hap_cov} ${hifi_median_hap_cov} 500000 20221102 > tig_marker_finder.log
echo "tig_marker_finder_v2 ${sample}_cnv_winsize10000_step10000_hq_merged_vs_hifi.txt ${illu_median_hap_cov} ${hifi_median_hap_cov} 2000000 20221204 > tig_marker_finder.log"
tig_marker_finder_v2 ${sample}_cnv_winsize10000_step10000_hq_merged_vs_hifi.txt ${illu_median_hap_cov} ${hifi_median_hap_cov} 2000000 20221204 > tig_marker_finder.log
#
#
wd=/your/working/path//
cd ${wd}
#
done
