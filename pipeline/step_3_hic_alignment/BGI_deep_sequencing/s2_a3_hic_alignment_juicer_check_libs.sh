# Here is the juicer cmd for checking the initial omni-c libs from mpgc: if they are good with >20~30% unique reads, they would be sequenced deeply at BGI

wd=/your/path/to//a3_hic_alignment/
cd ${wd}

# step 1. index reference

wd=/your/path/to//a3_hic_alignment/
cd ${wd}
#
ref_path="/your/path/to//a2_initial_assembly/"
#
for sample in A B C D E F G H I J; do
    cd ${ref_path}/sample_${sample}/hifiasm_asm/subset6_hifi_re_align_purge_dup/
    echo ${sample}": "
    ll purged.fa
    bsub -o index.log -q normal -R "rusage[mem=10000]" -R "span[hosts=1]" -M 10000 "bwa index -a bwtsw purged.fa"
    bsub -o index.log -q normal -R "rusage[mem=10000]" -R "span[hosts=1]" -M 10000 "samtools faidx purged.fa"
    bsub -o index.log -q normal -R "rusage[mem=10000]" -R "span[hosts=1]" -M 10000 "bowtie2-build purged.fa purged.fa"
    cd ${wd}
done

# step 2. prepare reads and configure juicer folders

wd=/your/path/to//a3_hic_alignment/
cd ${wd}

# A B C D E F G H I J
for sample in D; do
   #
   wd=/your/path/to//a3_hic_alignment/
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
   readpath=//your/path/to/omnic//reads/
   #
   R1=$(ls //your/path/to/omnic//reads/*/5424_${sample}_*R1*.fastq.gz)
   cat ${R1} > ${sample}_R1_mpgc.fastq.gz
   R2=$(ls //your/path/to/omnic//reads/*/5424_${sample}_*R2*.fastq.gz)
   cat ${R2} > ${sample}_R2_mpgc.fastq.gz
   #
   cd ../chrsizes
   ln -s ../../../../a2_initial_assembly/sample_${sample}/hifiasm_asm/subset6_hifi_re_align_purge_dup/purged.ctgsizes ln_purged.ctgsizes
   #
   cd ../references
   ln -s ../../../../a2_initial_assembly/sample_${sample}/hifiasm_asm/subset6_hifi_re_align_purge_dup/purged.fa ln_purged.fa
   ln -s ../../../../a2_initial_assembly/sample_${sample}/hifiasm_asm/subset6_hifi_re_align_purge_dup/purged.fa.amb ln_purged.fa.amb
   ln -s ../../../../a2_initial_assembly/sample_${sample}/hifiasm_asm/subset6_hifi_re_align_purge_dup/purged.fa.ann ln_purged.fa.ann
   ln -s ../../../../a2_initial_assembly/sample_${sample}/hifiasm_asm/subset6_hifi_re_align_purge_dup/purged.fa.bwt ln_purged.fa.bwt
   ln -s ../../../../a2_initial_assembly/sample_${sample}/hifiasm_asm/subset6_hifi_re_align_purge_dup/purged.fa.pac ln_purged.fa.pac
   ln -s ../../../../a2_initial_assembly/sample_${sample}/hifiasm_asm/subset6_hifi_re_align_purge_dup/purged.fa.sa ln_purged.fa.sa
   ln -s ../../../../a2_initial_assembly/sample_${sample}/hifiasm_asm/subset6_hifi_re_align_purge_dup/purged.fa.fai ln_purged.fa.fai
   #
   cd ../scripts
   cp -a -r /netscratch/dep_mercier/grp_schneeberger/bin/juicer/juicer-1.6/CPU_from_Meng/* ./
   cd ..
   #
   wd=/your/path/to//a3_hic_alignment/
   cd ${wd}
   #
done

# step 2. align hi-c reads to the ungrouped contigs

wd=/your/path/to//a3_hic_alignment/
cd ${wd}

# A B C D E F G H I J
for sample in D; do
   #
   wd=/your/path/to//a3_hic_alignment/
   cd ${wd}
   cd hic_${sample}
   cd 0_juicer_config
   #
   genome=./references/ln_purged.fa
   Juicer_scripts_dir=/your/path/to//a3_hic_alignment/hic_${sample}/0_juicer_config
   Juicer_config_dir=/your/path/to//a3_hic_alignment/hic_${sample}/0_juicer_config
   chrsizes=./chrsizes/ln_purged.ctgsizes
   #
   bsub -q ioheavy -n 20 -R "span[hosts=1] rusage[mem=24000]" -M 24000 -o juicer_${sample}_mpgc.log -e juicer_${sample}_mpgc.err "./scripts/juicer.sh -g ${sample}_mpgc -D ${Juicer_scripts_dir} -d ${Juicer_config_dir} -s none -t 20 -z ${genome} -p ${chrsizes}"
   #
   wd=/your/path/to//a3_hic_alignment/
   cd ${wd}
   #
done

# Conclusion 20220614: using C where shallow and deep sequencings are available and 75% uniquely aligned when shallow, >40% when deep: it is risky to deep sequence sample D, so we create a new lib for D.
# Rerun for a new omnic lib of D - 20220725:

# /path/to/java/jre1.8.0_181/bin/java -Xmx5000m -jar /netscratch/dep_mercier/grp_schneeberger/bin/Juicebox/juicebox_2.13.07.jar
