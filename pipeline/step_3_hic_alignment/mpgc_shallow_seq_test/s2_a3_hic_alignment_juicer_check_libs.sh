# Here is the juicer cmd for checking the initial omni-c libs from mpgc: if they are good with >20~30% unique reads, they would be sequenced deeply at BGI

wd=/working/path/
cd ${wd}

# step 1. index reference 

wd=/working/path/
cd ${wd}
#
ref_path="/path/to/assembly2_initial_assembly/"
# 
for sample in A B; do 
    cd ${ref_path}/sample_${sample}/hifiasm_asm/subset1_illu_re_align_purge_ovl/DNB_align/
    echo ${sample}": "
    ll clipped_${sample}_hifiasm.bp.p_utg.gfa.fa
    bsub -o index.log -q normal -R "rusage[mem=10000]" -R "span[hosts=1]" -M 10000 "bwa index -a bwtsw clipped_${sample}_hifiasm.bp.p_utg.gfa.fa"
    bsub -o index.log -q normal -R "rusage[mem=10000]" -R "span[hosts=1]" -M 10000 "samtools faidx clipped_${sample}_hifiasm.bp.p_utg.gfa.fa" 
    # bsub -o index.log -q normal -R "rusage[mem=10000]" -R "span[hosts=1]" -M 10000 "bowtie2-build clipped_${sample}_hifiasm.bp.p_utg.gfa.fa clipped_${sample}_hifiasm.bp.p_utg.gfa.fa"
    cd ${wd}
done

# step 2. prepare reads and configure juicer folders

wd=/working/path/
cd ${wd}

# A B C D E F G H I J
for sample in A B; do 
   #
   wd=/working/path/
   cd ${wd}
   mkdir hic_${sample}
   cd hic_${sample}
   mkdir 0_juicer_config
   cd 0_juicer_config
   #
   mkdir fastq chrsizes scripts references restriction_sites   
   #
   cd ./fastq
   #
   readpath=/path/to/reads/
   #
   R1=$(ls /path/to/reads/*/5424_${sample}_*R1*.fastq.gz)
   cat ${R1} > ${sample}_R1_mpgc.fastq.gz
   R2=$(ls /path/to/reads/*/5424_${sample}_*R2*.fastq.gz)
   cat ${R2} > ${sample}_R2_mpgc.fastq.gz   
   #
   cd ../chrsizes
   ln -s ../../../../../a2_initial_assembly/sample_${sample}/hifiasm_asm/subset1_illu_re_align_purge_ovl/DNB_align/clipped_${sample}_hifiasm.bp.p_utg.gfa.fa.ctgsizes ln_purged.ctgsizes
   #
   cd ../references
   ln -s ../../../../../a2_initial_assembly/sample_${sample}/hifiasm_asm/subset1_illu_re_align_purge_ovl/DNB_align/clipped_${sample}_hifiasm.bp.p_utg.gfa.fa ln_purged.fa
   ln -s ../../../../../a2_initial_assembly/sample_${sample}/hifiasm_asm/subset1_illu_re_align_purge_ovl/DNB_align/clipped_${sample}_hifiasm.bp.p_utg.gfa.fa.amb ln_purged.fa.amb
   ln -s ../../../../../a2_initial_assembly/sample_${sample}/hifiasm_asm/subset1_illu_re_align_purge_ovl/DNB_align/clipped_${sample}_hifiasm.bp.p_utg.gfa.fa.ann ln_purged.fa.ann
   ln -s ../../../../../a2_initial_assembly/sample_${sample}/hifiasm_asm/subset1_illu_re_align_purge_ovl/DNB_align/clipped_${sample}_hifiasm.bp.p_utg.gfa.fa.bwt ln_purged.fa.bwt
   ln -s ../../../../../a2_initial_assembly/sample_${sample}/hifiasm_asm/subset1_illu_re_align_purge_ovl/DNB_align/clipped_${sample}_hifiasm.bp.p_utg.gfa.fa.pac ln_purged.fa.pac
   ln -s ../../../../../a2_initial_assembly/sample_${sample}/hifiasm_asm/subset1_illu_re_align_purge_ovl/DNB_align/clipped_${sample}_hifiasm.bp.p_utg.gfa.fa.sa ln_purged.fa.sa
   ln -s ../../../../../a2_initial_assembly/sample_${sample}/hifiasm_asm/subset1_illu_re_align_purge_ovl/DNB_align/clipped_${sample}_hifiasm.bp.p_utg.gfa.fa.fai ln_purged.fa.fai   
   #
   cd ../scripts
   cp -a -r /netscratch/dep_mercier/grp_schneeberger/bin/juicer/juicer-1.6/CPU_from_Meng/* ./
   # cp -a -r /netscratch/dep_mercier/grp_schneeberger/bin/juicer/CPU/* ./   
   cd ..
   #
   wd=/working/path/
   cd ${wd}   
   #
done

# step 2. align hi-c reads to the ungrouped contigs 

wd=/working/path/
cd ${wd}

# A B C D E F G H I J
for sample in A B; do 
   #
   wd=/working/path/
   cd ${wd}
   cd hic_${sample}
   cd 0_juicer_config
   #
   genome=./references/ln_purged.fa
   Juicer_scripts_dir=/working/path/hic_${sample}/0_juicer_config
   Juicer_config_dir=/working/path/hic_${sample}/0_juicer_config
   chrsizes=./chrsizes/ln_purged.ctgsizes
   #  --cleanup not in v1.6
   bsub -q ioheavy -n 20 -R "span[hosts=1] rusage[mem=24000]" -M 24000 -o juicer_${sample}_mpgc.log -e juicer_${sample}_mpgc.err "./scripts/juicer.sh -g ${sample}_mpgc -D ${Juicer_scripts_dir} -d ${Juicer_config_dir} -s none -t 20 -z ${genome} -p ${chrsizes}"   
   #
   wd=/working/path/
   cd ${wd}   
   #
done

# Compare uniquely align reads, inter-chr reads with summary_OmniC.xlsx, to estimate complexlity of Omni-C lib (to decide going for deeper sequencing or not).


