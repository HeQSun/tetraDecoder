# Here we 1. purge overlaps between unitigs with gfa file from hifiasm (v0p7), and 2. align reads to the cleaned assembly.
#    sub-script: 
#       purge_overlaps
#       fasta_length
#       calc_CN50.pl
#       new_bowtie2_align_DNB_clipped_unitigs_cleanreads.sh
# 01 SEP 2022
# Hequan Sun
#

# 1. align short reads - do this for both raw and purged assemblies - need some updates on scripts.

# where you have initial assemblies for the cultivars
wd=/your/working/path/ 
cd ${wd}
# C
for sample in R T; do
#
wd=/your/working/path/
cd ${wd}/sample_${sample}
cd hifiasm_asm_v0p7
#
subset_id=1
mkdir subset${subset_id}_illu_re_align_purge_ovl
cd subset${subset_id}_illu_re_align_purge_ovl
mkdir DNB_align
cd DNB_align
#
nthread=20
#
bsub -q ioheavy -n ${nthread} -R "span[hosts=1] rusage[mem=64000]" -M 64000 -o subset_DNB_alignment_process.log -e subset_DNB_alignment_process.err "/path/to/new_bowtie2_align_DNB_clipped_unitigs_cleanreads_hifiasm_v0p7.sh ${nthread} ${sample} ${subset_id} > new_bowtie2_align_DNB_clipped_unitigs_cleanreads_hifiasm_v0p7.log"
#nohup /path/to/new_bowtie2_align_DNB_clipped_unitigs.sh ${nthread} ${sample} ${subset_id} &> subset_DNB_alignment_process.log&
#
cd ${wd}
#
done

# R visualize: need the coverage file from CNV_HQ_v3 for both raw and purged assemblies!

cd /your/working/path/
sample_x_cnv_initial_contigs_winsize10kb_step10kb_depth_DNB_purged_unitig_level_hifiasm_v0p7.R
sed -e '$s/$/\n/' -s  sample_*/hifiasm_asm_v0p7/subset1_illu_re_align_purge_ovl/DNB_align/sample_*_tig_ratio > sample_10_tig_ratio_summary.txt

# 2. align hifi reads 

# R T
for sample in D; do
#
subset_id=1
wd=/your/working/path/
cd ${wd}/sample_${sample}
cd hifiasm_asm_v0p7
cd subset${subset_id}_illu_re_align_purge_ovl
#
mkdir HIFI_align
cd HIFI_align
#
this_path="./"
nthread=20
genome=../DNB_align/clipped_${sample}_hifiasm.p_utg.gfa.fa
reads=../../../../../a1_HiFi_stats/sample_${sample}/sample_${sample}_zcat_HiFi_trim.fastq.gz
# run from hpc
bsub -q ioheavy -n ${nthread} -R "span[hosts=1] rusage[mem=64000]" -M 64000 -o subset_Hifi_alignment_process_v3.log -e subset_Hifi_alignment_process_v3.err "../../../../subset_Hifi_alignment_process_v3.sh ${sample} ${this_path} ${subset_id} ${nthread} ${genome} ${reads} > subset_Hifi_alignment_process_details.log"
# run from dell - no bsub
# nohup ../../../../subset_Hifi_alignment_process_v2.sh ${sample} ${this_path} ${subset_id} ${nthread} ${genome} ${reads} &> subset_Hifi_alignment_process_details.log&
#
#
cd ${wd}/
# 
done





