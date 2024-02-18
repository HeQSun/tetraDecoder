# Here we use variants called in hifi to correct over-purging, or to recover haplotypes in regions
#    where we can observe variations in hifi data.

# you should have prepared bams of both illumina and hifi during coverage analysis

#
# step 1: call variants with short reads - to find non-haplotigs but with "sufficient" variations: bcftools 1.9
#         note, reads aligned in the last coverage analysis script
#
# 1.1 perform illumina variant calling
#
wd=/your/working/path/
cd ${wd}
#
for sample in A B C D E F G H I J; do
#
wd=/your/working/path/
cd ${wd}/sample_${sample}
cd hifiasm_asm_v0p7
#
subset_id=2
cd subset${subset_id}_illu_re_align_purge_ovl
cd DNB_align
#
mkdir variant_call_sc
cd variant_call_sc
# files: ploidy,fasta,bam
cp -a ../../../../../ploidy .
genome=../clipped${subset_id}_${sample}_hifiasm.p_utg.gfa.fa
bam=../${sample}_PE_sorted.cram
#
sed 's/\t/\t1\t/g' ../clipped${subset_id}_${sample}_hifiasm.p_utg.gfa.fa.ctgsizes | sort -k3,3 --random-sort > clipped${subset_id}_${sample}_interval
split -l 800 clipped${subset_id}_${sample}_interval clipped${subset_id}_${sample}_interval_split
for interval in clipped${subset_id}_${sample}_interval_split*; do mv ${interval} ${interval}.txt; done
#
for interval in *split*.txt; do
    bsub -q normal -R "span[hosts=1] rusage[mem=3000]" -M 4000 -o sc_variant_call.log -e sc_variant_call.err "bcftools mpileup --regions-file ${interval} -d 1000 -Ou -f ${genome} ${bam} | bcftools call -A -m -Ov -v --ploidy-file ./ploidy > ${interval}.vcf"
done
#
wd=/your/working/path/
cd ${wd}
#
done

# 1.2 merge illumina sub-vcfs for each sample and find raw illumina SNPs
#     result of raw illumina SNPs -- tool SHOREmap used: www.shoremap.org
#	48351 A_shoremap_converted/1_converted_variant.txt
#	36144 B_shoremap_converted/1_converted_variant.txt
#	40718 C_shoremap_converted/1_converted_variant.txt
#	58311 D_shoremap_converted/1_converted_variant.txt
#	40982 E_shoremap_converted/1_converted_variant.txt
#	44603 F_shoremap_converted/1_converted_variant.txt
#	46487 G_shoremap_converted/1_converted_variant.txt
#	35117 H_shoremap_converted/1_converted_variant.txt
#	48301 I_shoremap_converted/1_converted_variant.txt
#	58967 J_shoremap_converted/1_converted_variant.txt
#
wd=/your/working/path/
cd ${wd}
#
for sample in A B C D E F G H I J; do
#
wd=/your/working/path/
cd ${wd}/sample_${sample}
cd hifiasm_asm_v0p7
#
subset_id=2
cd subset${subset_id}_illu_re_align_purge_ovl
cd DNB_align
cd variant_call_sc
#
cat clipped*split*.vcf | grep -E -v '#|INDEL' > z_merged_clipped2_${sample}_interval.vcf
SHOREmap convert --marker z_merged_clipped2_${sample}_interval.vcf --folder ${sample}_shoremap_converted --indel-size 0 --min-AF 0.01 -no-c -no-r > SHOREmap_convert.log
wc -l ${sample}_shoremap_converted/*.txt
#
wd=/your/working/path/
cd ${wd}
#
done

# step 2: call variants with hifi reads:
#
# 2.1 split hifi bam
wd=/your/working/path/
cd ${wd}
#
for sample in A B C D E F G H I J; do
#
wd=/your/working/path/
cd ${wd}/sample_${sample}
cd hifiasm_asm_v0p7
#
subset_id=2
cd subset${subset_id}_illu_re_align_purge_ovl
cd HIFI_align
#
mkdir variant_call_sc
cd variant_call_sc
#
genome=../../DNB_align/clipped${subset_id}_${sample}_hifiasm.p_utg.gfa.fa
bam=../ptt${sample}_hiasm_ref_pilon_subset${subset_id}.bam
#
for interval in ../../DNB_align/variant_call_sc/clipped*_*_interval_split*.txt; do
    this_base=$(basename "$interval")
    this_base_no_suffix=${this_base%.*}
    cat ${interval} > ${this_base_no_suffix}.bed
    bsub -q normal -R "span[hosts=1] rusage[mem=1000]" -M 1000 -o split_bam.log -e split_bam.err "samtools view -h -b -L ${this_base_no_suffix}.bed ${bam} > ${this_base_no_suffix}.bam; samtools index ${this_base_no_suffix}.bam"
done
#
wd=/your/working/path/
cd ${wd}
#
done
#
# 2.2 - perform hifi variant calling: singularity version 3.6.4 plus deepvariant v1.4.0
#         otherwise error: Set an empty environment as LD_LIBRARY_PATH may mix dependencies, just rely only on the library cache or its own lookup mechanism.
#                see here: https://github.com/apptainer/singularity/pull/5669, export LD_LIBRARY_PATH=""
#         caution on settings for running "singularity run deepvariant"
#
wd=/your/working/path/
cd ${wd}
#
for sample in A B C D E F G H I J; do
#
wd=/your/working/path/
cd ${wd}/sample_${sample}
cd hifiasm_asm_v0p7
#
subset_id=2
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
for interval in ../../DNB_align/variant_call_sc/clipped*_*_interval_split*.txt; do
    this_base=$(basename "$interval")
    this_base_no_suffix=${this_base%.*}
    #
    # export TMPDIR="$PWD/zdv1_${this_base_no_suffix}_tmp/";export LC_CTYPE="";export LANG="";
    this_flag="zdv1"
    bsub -q ioheavy -n 4 -R "rusage[mem=16000]" -M 16000 -o ${this_flag}_${this_base_no_suffix}.log -e ${this_flag}_${this_base_no_suffix}.err "mkdir ${this_flag}_${this_base_no_suffix}_tmp; export TMPDIR=\"$PWD/${this_flag}_${this_base_no_suffix}_tmp/\";export LC_CTYPE=\"\";export LANG=\"\"; export LD_LIBRARY_PATH=\"\"; singularity run --bind ${PWD} /srv/netscratch/dep_mercier/grp_schneeberger/software/deepvariant/deepvariant.img /opt/deepvariant/bin/run_deepvariant --model_type PACBIO --ref clipped2_${sample}_hifiasm.p_utg.gfa.fa --reads ${this_base_no_suffix}.bam --intermediate_results_dir ./${this_flag}_${this_base_no_suffix}_tmp --output_vcf ${this_flag}_${this_base_no_suffix}.vcf.gz --num_shards 4 "
    #
done
#
wd=/your/working/path/
cd ${wd}
#
done
#
#
# 2.3 merge hifi sub-vcfs for each sample and find raw hifi SNPs
#     result of raw hifi SNPs
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
#
wd=/your/working/path/
cd ${wd}
#
for sample in A B C D E F G H I J; do
#
wd=/your/working/path/
cd ${wd}/sample_${sample}
cd hifiasm_asm_v0p7
#
subset_id=2
cd subset${subset_id}_illu_re_align_purge_ovl
cd HIFI_align
cd variant_call_sc
#
cat zdv1_clipped2_*_interval_split*.vcf.gz | gunzip | grep -E -v '#|PASS' | awk '$5=="A" || $5=="C" || $5=="G" ||$5=="T" ' | awk '$4=="A" || $4=="C" || $4=="G" ||$4=="T" ' > z_merged_clipped2_${sample}_interval_hifi.vcf
wc -l z_merged_clipped2_${sample}_interval_hifi.vcf
#
wd=/your/working/path/
cd ${wd}
#
done
#

# step 3. align raw contigs to purge2 contigs - if
# 		a purge2 contig is aligned by multiple raw contigs, and
# 		there are both DNB and HIFI SNPs,
#		re-collect the contig covering alternative alleles into purge2 contigs, to make a purge2_refined_ctgs.fa
#
# 3.1 align raw contigs to purge2 contigs

wd=/your/working/path/
cd ${wd}
#
for sample in A B C D E F G H I J; do
#
wd=/your/working/path/
cd ${wd}/sample_${sample}
cd hifiasm_asm_v0p7/subset2_illu_re_align_purge_ovl/
#
mkdir raw_asm_to_purge2_ctg
cd raw_asm_to_purge2_ctg
purge2_ref=../DNB_align/clipped2_${sample}_hifiasm.p_utg.gfa.fa
raw_asm=../../${sample}_hifiasm.p_utg.gfa.fasta
ll ${purge2_ref} ${raw_asm}
thread=10
#
bsub -q multicore20 -n 10 -R "span[hosts=1] rusage[mem=100000]" -M 100000 -o minimap2.log -e minimap2.err "minimap2 -ax asm20 -t ${thread} -N 1 --secondary=no ${purge2_ref} ${raw_asm} | samtools view -@ ${thread} -bS - | samtools sort -@ ${thread} -o ${sample}_raw_asm_aligned_to_purge2_ctg.bam -; samtools index -@ ${thread} ${sample}_raw_asm_aligned_to_purge2_ctg.bam"
#
wd=/your/working/path/
cd ${wd}
#
done

# step 4. extract long reads according to local variants.
#
wd=/your/working/path/
cd ${wd}
#
for sample in A B C D E F G H I J; do
#
wd=/your/working/path/
cd ${wd}/sample_${sample}
cd hifiasm_asm_v0p7/
#
mkdir subset3_local_asm
cd subset3_local_asm
#
#
#      marker format: contig_id   position .   maternal paternal  optional
#                     utg000103lc 727428   .   A        T         0.2 RefCall .....
marker=../subset2_illu_re_align_purge_ovl/HIFI_align/variant_call_sc/z_merged_clipped2_${sample}_interval_hifi.vcf
gap_size=20000 # bp
cram=../subset2_illu_re_align_purge_ovl/HIFI_align/ptt${sample}_hiasm_ref_pilon_subset2.cram
# run this first; 1 second each
/path/to/dev_bin/snp_interval_finder ${marker} ${gap_size} ${sample}_var; awk '$2>=10' ${sample}_var_snp_clusters.txt | sed 's/:/\t/g' | cut -f1 | sort | uniq > contigs_with_var_for_local_asm.list
#
mkdir z1_local_read_extract
cd z1_local_read_extract
# 4.1 prepare commands
ctg_cnt=0
>${sample}_bam_extracter_all.sh
while read ctg; do echo "$((ctg_cnt=ctg_cnt+1)): ${sample}-${ctg}"; echo "samtools view ../${cram} ${ctg} | /path/to/dev_bin//bam_extracter ../${marker} ${ctg}_local_hifi - > bam_extracter_${ctg}.log; gzip ${ctg}_local_hifi_extracted_reads.fq" >> ${sample}_bam_extracter_all.sh; done < ../contigs_with_var_for_local_asm.list
# 4.2 split cmds and submit
split -l 100 ${sample}_bam_extracter_all.sh ${sample}_bam_extracter_sub
for si in ${sample}_bam_extracter_sub*; do mv ${si} ${si}.sh; bsub -q normal -R "span[hosts=1] rusage[mem=6000]" -o bam_extracter.log -e bam_extracter.err "bash ./${si}.sh"; done
#
cd ..
#
wd=/your/working/path/
cd ${wd}
#
done
#

# step 5. local assembly with var-extracted hifi reads.
#
wd=/your/working/path/
cd ${wd}
#
for sample in A B C D E F G H I J; do
#
wd=/your/working/path/
cd ${wd}/sample_${sample}
cd hifiasm_asm_v0p7/
#
cd subset3_local_asm
mkdir z2_local_asm
cd z2_local_asm
ctg_cnt=0
>${sample}_local_assembly_all.sh
while read ctg; do echo "$((ctg_cnt=ctg_cnt+1)): ${sample}-${ctg}"; echo "cd /your/working/path/sample_${sample}/hifiasm_asm_v0p7/subset3_local_asm/z2_local_asm; mkdir local_${ctg}; cd local_${ctg}; flye --threads 1 -m 1000 --keep-haplotypes --scaffold --pacbio-hifi ../../z1_local_read_extract/${ctg}_local_hifi_extracted_reads.fq.gz --read-error 0.03 --out-dir . > flye.log" >> ${sample}_local_assembly_all.sh; done < ../contigs_with_var_for_local_asm.list
# 4.2 split cmds and submit
split -l 50 ${sample}_local_assembly_all.sh ${sample}_local_assembly_sub
for si in ${sample}_local_assembly_sub*; do mv ${si} ${si}.sh; bsub -q ioheavy -R "span[hosts=1] rusage[mem=6000]" -o local_assembly.log -e local_assembly.err "bash ./${si}.sh"; done
#
cd ..
#
wd=/your/working/path/
cd ${wd}
#
done
#
#

# step 6. combine locally assembled contigs
#
wd=/your/working/path/
cd ${wd}
# A B C D E F G H I J
for sample in A B C D E F G H I J; do
#
wd=/your/working/path/
cd ${wd}/sample_${sample}
cd hifiasm_asm_v0p7/
#
cd subset3_local_asm
# there are cases where no local contigs assembled - such cases are discarded.
echo ${sample}
while read ctg; do sed -i "s/>contig/>${ctg}_contig/g" z2_local_asm/local_${ctg}/assembly.fasta; done < contigs_with_var_for_local_asm.list
#
cat ./z2_local_asm/*/assembly.fasta > ${sample}_assembly_local.fasta
ll ${sample}_assembly_local.fasta
#
wd=/your/working/path/
cd ${wd}
#
done
#
#

# step 7. create fasta with sequence selection
# Note: cannot multicore20/40 as TERM_THREADLIMIT: job killed after reaching LSF thread limit.
# 1.1 align short reads

wd=/your/working/path/
cd ${wd}
# A B C D E F G ioheavy
for sample in H I J; do
#
wd=/your/working/path/
cd ${wd}/sample_${sample}
cd hifiasm_asm_v0p7
#
subset_id=3
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
