This is the pipeline explaining haplotyping in a tetraploid genome using Hi-C and long read sequencing.
=

#### Publicly available tools need to be installed. Note: zlib.h is required by some tools (for the purpose of zipping files), please install accordingly：

* [bwa 0.7.17](https://github.com/lh3/bwa/releases)
* [samtools 1.9](https://github.com/samtools/)
* [bowtie2 2.2.8](https://github.com/BenLangmead/bowtie2)
* [minimap2 2.20-r1061](https://github.com/lh3/minimap2)
* [ALLHiC, including filterBAM_forHiC.pl](https://github.com/tangerzhang/ALLHiC)
* [ALLHiC-dev](https://github.com/sc-zhang/ALLHiC_components)
* [gmap](http://research-pub.gene.com/gmap/src/gmap-gsnap-2012-07-20.v3.tar.gz)
* [hifiasm 0.7](https://github.com/chhylp123/hifiasm)

#### Developed tools/scripts in this work (Check [INSTALL](https://github.com/HeQSun/tetraDecoder/blob/main/INSTALL) for installation. Installation tested on linux distribution "Debian GNU/Linux 9 (stretch)" with x86_64 cpu architecture).

 * fasta_length-----------------# calculate length of sequences in fasta format
 * fasta_name_selecter--------# select a subset of sequences with sequence name from a fasta file
 * ref_linkage_grouper---------# based on alignment of contigs to a reference genomes, separate contigs into linkage groups
 * omnic_read_extracter-------# given the groups of contigs, separate Hi-C read pairs into the groups
 * long_read_separator_v2-----# given the groups of contigs, separate long reads into the groups
 * hic_binning------------------# this is the core of this pipeline, which separates contigs into haplotypes.
 * z_suppl_sample_10cultivars_checking_after_hic_binining.R # check size of groupped contigs after running Hi-C based contig binning

#### Besides, basic tools cat, grep, awk and sed should be installed in the system.

#### Step.0 Prepare data

##### step 0.1. prepare sequencing data of the focal genome. All data are from a tetraploid (potato cultivar) of interest, including PacBio HiFi reads, Hi-C reads from somatic tissues. In this example pipeline, suppose all raw sequencing data are collected in the path below (Illumina reads for coverage analysis of contig window markers are not given here, please refer to [Sun_and_Jiao_et_al_2021](https://nature.com/articles/s41588-022-01015-0) ).

    read_path=/your/work/directory/reads/
    cd ${read_path}

##### PacBio HiFi [available here](https://www.ncbi.nlm.nih.gov/sra/SRX11512735[accn])

    ls -l 4396_A_CCS.fastq

##### Hi-C reads, with 3 raw subsets for 'Otava' (sample 'O'):

    wget https://websafe.mpipz.mpg.de/d/0r5CGxtrn3/V300052863_L01_read_1.fq.gz
    wget https://websafe.mpipz.mpg.de/d/7BPV2CMwfg/V300052863_L01_read_2.fq.gz
    wget https://websafe.mpipz.mpg.de/d/acCNzIVNxS/V300052863_L02_read_1.fq.gz
    wget https://websafe.mpipz.mpg.de/d/iRwZYPqxRC/V300052863_L02_read_2.fq.gz
    wget https://websafe.mpipz.mpg.de/d/AwAq1VeuXi/V300052863_L03_read_1.fq.gz
    wget https://websafe.mpipz.mpg.de/d/Ze1rcfPZQF/V300052863_L03_read_2.fq.gz
    wget https://websafe.mpipz.mpg.de/d/GcadeWRrvD/md5.txt # check md5 by yourself

    mv V300052863_L01_read_1.fq.gz O_L1_1.fq.gz
    mv V300052863_L01_read_2.fq.gz O_L1_2.fq.gz
    mv V300052863_L02_read_1.fq.gz O_L2_1.fq.gz
    mv V300052863_L02_read_2.fq.gz O_L2_2.fq.gz
    mv V300052863_L03_read_1.fq.gz O_L3_1.fq.gz
    mv V300052863_L03_read_2.fq.gz O_L3_2.fq.gz

    ls -l O_L1_1.fq.gz
    ls -l O_L1_2.fq.gz
    ls -l O_L2_1.fq.gz
    ls -l O_L2_2.fq.gz
    ls -l O_L3_1.fq.gz
    ls -l O_L3_2.fq.gz

##### step 0.2. prepare DM v6p1 reference data, with gff and cdna of gene models.

    ref_seq_path=/your/work/directory/ref_dm6p1/
    cd ${ref_seq_path}
    wget http://spuddb.uga.edu/data/DM_1-3_516_R44_potato_genome_assembly.v6.1.fa.gz
    gunzip DM_1-3_516_R44_potato_genome_assembly.v6.1.fa.gz
    dm_ref="DM_1-3_516_R44_potato_genome_assembly.v6.1.fa"
    fasta_length ${dm_ref} | grep '>' | sed 's/>//g' > DM_1-3_516_R44_potato_genome_assembly.v6.1.chrsizes
    cat     DM_1-3_516_R44_potato_genome_assembly.v6.1.chrsizes | grep -v 'scaffold' | nl | sed 's/ //g' > DM_1-3_516_R44_potato_genome_assembly.v6.1_main12.chrsizes
    cut -f2 DM_1-3_516_R44_potato_genome_assembly.v6.1_main12.chrsizes                                   > DM_1-3_516_R44_potato_genome_assembly.v6.1_main12.chrids

##### save gff and cdna data of gene model for each lg [available here: DM-lg-wise-cda-gff](https://mega.nz/folder/GscjQawR#QgwKEbbVjIdghp7r0xeWIg)

#### Step 1. create and index preliminary assembly, and align HiFi reads to the assembly

##### step 1.1. for the generation of the initial assembly.

    asm_seq_path=/your/work/directory/assembly/
    cd ${asm_seq_path}
    sample="O"

    hifi=/path/to/s0_reads/4396_A_CCS.fastq
    hifiasm -t 10 -o otava ${hifi} >hifiasm.log

##### #Note, redundant contigs representing the same genomic regions were purged - please check supplementary information: section "Initial tetraploid genome assembly, polishing and purging" for details, in our previous work [Sun_and_Jiao_et_al_2021](https://nature.com/articles/s41588-022-01015-0). Here we provide the purged version of the initial assembly for the test [available here](https://mega.nz/folder/GktXEYCR#F3I8uTKvKO0Fu8VY2yc2WA).
    
    mv HiFiasm_ref_6366long_ctgs_selected.fasta clipped4_${sample}_hifiasm.p_utg.gfa.fa

    bwa index -a bwtsw clipped4_${sample}_hifiasm.p_utg.gfa.fa
    samtools faidx clipped4_${sample}_hifiasm.p_utg.gfa.fa
    bowtie2-build --threads 8 clipped4_${sample}_hifiasm.p_utg.gfa.fa clipped4_${sample}_hifiasm.p_utg.gfa.fa

##### step 1.2. align HiFi reads to the contigs

    sample=O
    genome=clipped4_${sample}_hifiasm.p_utg.gfa.fa
    otavahifi=/your/work/directory/reads/4396_A_CCS.fastq

    minimap2 -ax map-pb -t 24 -N 1 --secondary=no ${genome} ${otavahifi} | samtools view -@ 24 -bS - | samtools sort -@ 24 -o ${sample}_hifiasm.bam -
    samtools flagstat -@ 24 ${sample}_hifiasm.bam
    samtools index -@ 24 ${sample}_hifiasm.bam

#### step 2: prepare new window marker ([available here](https://github.com/HeQSun/tetraDecoder/tree/main/aux_intermediate_data)) : for generating this file please check [Step 4-6 here](https://github.com/schneebergerlab/GameteBinning_tetraploid/blob/main/Pipeline.md)

    marker_path=/your/work/directory/win_marker/
    cd ${marker_path}
    ls -l ${sample}_cnv_winsize10000_step10000_hq_merged_vs_hifi_markers_20221102_wsize500kb_final.txt

#### step 3. group contigs with a reference genome: here we use DM v6.1: ([available here](http://solanaceae.plantbiology.msu.edu/dm_v6_1_download.shtml)) ([ref](https://academic.oup.com/gigascience/article/9/9/giaa100/5910251?searchresult=1#207670451))

    wd=/your/work/directory/
    cd ${wd}
    mkdir a4_alignment_based_linkage_grouping

##### step 3.1. align assembled contigs to DM genonme

    for sample in O; do
        group_lg_path=/your/work/directory/a4_alignment_based_linkage_grouping/
        cd ${group_lg_path}
        mkdir zhifiasm_v0p7_dm_${sample}_asm
        cd zhifiasm_v0p7_dm_${sample}_asm
        dm_ref=/your/work/directory/ref_dm6p1/DM_1-3_516_R44_potato_genome_assembly.v6.1.fa
        assembly=/your/work/directory/assembly/clipped4_${sample}_hifiasm.p_utg.gfa.fa
        thread=4
        # paf format: for grouping
        minimap2 -cx asm20 -t ${thread} ${dm_ref} ${assembly} > ${sample}_against_dm.paf
        # bam format: for later allelic check
        minimap2 -ax asm20 -t ${thread} ${dm_ref} ${assembly} | samtools view -@ ${thread} -bS - | samtools sort -@ ${thread} -o ${sample}_against_dm.bam - 
        cd ${wd}
    done

##### step 3.2. separate alignment to linkage group => DM_based.log

    for sample in O; do
        group_lg_path=/your/work/directory/a4_alignment_based_linkage_grouping/
        cd ${group_lg_path}
        cd zhifiasm_v0p7_dm_${sample}_asm
        ref_linkage_grouper ${sample}_against_dm.paf /your/work/directory/ref_dm6p1/DM_1-3_516_R44_potato_genome_assembly.v6.1_main12.chrids > DM_based.log
    done

##### step 3.3. get grouping details and summary

    for sample in O; do
        group_lg_path=/your/work/directory/a4_alignment_based_linkage_grouping/
        cd ${group_lg_path}
        cd zhifiasm_v0p7_dm_${sample}_asm
        grep 'res' DM_based.log   | sed 's/   //g' > dm_res_grouping_details_${sample}.txt
        grep 'final' DM_based.log | sed 's/   //g' > dm_res_total_group_size_${sample}.txt
    done

##### step 3.4. get chr-wise group of contigs for omni-c read extraction - this is for omnic_read_extracter => dm_res_grouping_details_${sample}_${lg}.txt

    for sample in O; do
       group_lg_path=/your/work/directory/a4_alignment_based_linkage_grouping/
       cd ${group_lg_path}
       cd zhifiasm_v0p7_dm_${sample}_asm
       for lg in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12; do
           awk -v chr="$lg" '$7==chr ' dm_res_grouping_details_${sample}.txt | cut -d' ' -f2,4,7,11 | sed 's/ /\t/g' > dm_res_grouping_details_${sample}_${lg}.txt
       done
       cat dm_res_grouping_details_${sample}_*.txt | cut -d' ' -f1 | sort | uniq -c | wc -l
       cat dm_res_grouping_details_${sample}_*.txt | wc -l
    done

#### step 4. Hi-C read alignment to ungrouped contigs, and extract to lg-groups

##### step 4.1. align hi-c reads to the ungrouped contigs

    cd /your/work/directory/
    mkdir a3_hic_alignment

    for sample in O; do
       wd=/your/work/directory/a3_hic_alignment/
       cd ${wd}
       mkdir hic_${sample}
       cd hic_${sample}
       mkdir bwa_align
       cd bwa_align
       genome=/your/work/directory/assembly/clipped4_${sample}_hifiasm.p_utg.gfa.fa
       readpath=/your/work/directory/reads/
       Rs=`ls ${readpath}*_*.fq.gz`
       for readset in ${sample}_L1 ${sample}_L2 ${sample}_L3; do  
           R1="${readpath}/${readset}_1.fq.gz"
           R2="${readpath}/${readset}_2.fq.gz"
           nthreads=30
           bowtie2 -x ${genome} -1 ${R1} -2 ${R2} -p ${nthreads} --local | samtools view -@ ${nthreads} -bS - | samtools sort -n -@ ${nthreads} -o ${readset}_PE_RN_sorted_local.bam -
           samtools view -@ ${nthreads} -C -T ${genome} -o ${readset}_PE_RN_sorted_local.cram ${readset}_PE_RN_sorted_local.bam
           # rm *.bam # you can do this if space needs to be saved.
       done
    done

##### step 4.2. extract reads to ref-lg-grouped contigs

    cd /your/work/directory/a3_hic_alignment

    for sample in O; do
       wd=/your/work/directory/a3_hic_alignment/
       cd ${wd}
       cd hic_${sample}
       cd bwa_align
       for readset in ${sample}_L1 ${sample}_L2 ${sample}_L3; do
           for lg in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12; do
               group_lg_path=/your/work/directory/a4_alignment_based_linkage_grouping/
               group=${group_lg_path}/dm_res_grouping_details_${sample}_${lg}.txt
               samtools view -F 3840 ${readset}_PE_RN_sorted_local.cram | omnic_read_extracter - ${group} ${readset}_${lg}
               samtools sort -n -@ 1 -o ${readset}_${lg}_extract.bam ${readset}_${lg}_extract.sam
               rm ${readset}_${lg}_extract.sam
           done
       done
    done

##### step 4.3. extract paired-end hic reads to fastqs to assigned to each linkage group.

    cd /your/work/directory/a3_hic_alignment

    for sample in O; do
       wd=/your/work/directory/a3_hic_alignment/
       cd ${wd}
       cd hic_${sample}
       cd bwa_align
       for readset in ${sample}_L1 ${sample}_L2 ${sample}_L3; do
           for lg in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12; do        
               samtools fastq -1 ${readset}_${lg}_extract_R1.fastq.gz -2 ${readset}_${lg}_extract_R2.fastq.gz -0 ${readset}_${lg}_OTHER1.fastq.gz -s ${readset}_${lg}_OTHER2.fastq.gz -n -F 0x900 ${readset}_${lg}_extract.bam
           done
       done
    done

#### step 5. Hi-C read alignment to each lg

##### step 5.1. extract lg-wise assigned contigs (based on DM reference)

    cd /your/work/directory/a3_hic_alignment

    for sample in O; do
        wd=/your/work/directory/a3_hic_alignment/
        cd ${wd}
        cd hic_${sample}
        cd bwa_align
        mkdir lg_wise_contigs
        cd lg_wise_contigs
        genome=/your/work/directory/assembly/clipped4_${sample}_hifiasm.p_utg.gfa.fa
        for lg in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12; do
            group_lg_path=/your/work/directory/a4_alignment_based_linkage_grouping/
            cut -f 1 ${group_lg_path}/zhifiasm_v0p7_dm_${sample}_asm/dm_res_grouping_details_${sample}_${lg}.txt > dm_res_grouping_details_${sample}_${lg}_ids.txt
            fasta_name_selecter ${genome} dm_res_grouping_details_${sample}_${lg}_ids.txt
            mv clipped4_${sample}_hifiasm.p_utg.gfa_selected.fa clipped4_${sample}_${lg}.fa
        done
       #
    done

##### step 5.2. index lg-wise chrs for bwa

    wd=/your/work/directory/a3_hic_alignment/
    cd ${wd}

    for sample in O; do
        wd=/your/work/directory/a3_hic_alignment/
        cd ${wd}
        cd hic_${sample}
        cd bwa_align
        cd lg_wise_contigs
        for lg in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12; do
            bwa index -a bwtsw clipped4_${sample}_${lg}.fa
            samtools faidx clipped4_${sample}_${lg}.fa
            bowtie2-build --threads 8 clipped4_${sample}_${lg}.fa clipped4_${sample}_${lg}.fa
        done
    done

##### step 5.3. bwa align to initial assembly of 48 chrs

    wd=/your/work/directory/a3_hic_alignment/
    cd ${wd}

    for sample in O; do
        #
        wd=/your/work/directory/a3_hic_alignment/
        cd ${wd}
        cd hic_${sample}
        cd bwa_align
        mkdir lg_wise_contigs_read_align
        cd lg_wise_contigs_read_align
        for readset in ${sample}_L1 ${sample}_L2 ${sample}_L3; do
            for lg in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12; do
                mkdir align_${lg}
                cd align_${lg}
                genome=../../lg_wise_contigs/clipped4_${sample}_${lg}.fa
                R1="../../${readset}_${lg}_extract_R1.fastq.gz"
                R2="../../${readset}_${lg}_extract_R2.fastq.gz"
                bwa aln -t 8 ${genome} ${R1} > ${readset}_${lg}_R1.sai
                bwa aln -t 8 ${genome} ${R2} > ${readset}_${lg}_R2.sai
                cd ..
            done
        done
    done

##### step 5.4. merge read pairs

    wd=/your/work/directory/a3_hic_alignment/
    cd ${wd}

    for sample in O; do
        for lg in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12; do
            for readset in ${sample}_L1 ${sample}_L2 ${sample}_L3; do
                wd=/your/work/directory/a3_hic_alignment/
                cd ${wd}
                cd hic_${sample}
                cd bwa_align
                cd lg_wise_contigs_read_align
                cd align_${lg}
                genome=../../lg_wise_contigs/clipped4_${sample}_${lg}.fa
                R1="../../${readset}_${lg}_extract_R1.fastq.gz"
                R2="../../${readset}_${lg}_extract_R2.fastq.gz"
                bwa sampe ${genome} ${readset}_${lg}_R1.sai ${readset}_${lg}_R2.sai ${R1} ${R2} | samtools view -b -F12 -o ${readset}_${lg}_bwa_aln.bam
                samtools flagstat ${readset}_${lg}_bwa_aln.bam > ${readset}_${lg}_bwa_aln_bam_flagstat.txt
            done
        done
    done

##### step 5.5. filter reads in bam

    wd=/your/work/directory/a3_hic_alignment/
    cd ${wd}

    for sample in O; do
        for lg in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12; do
            wd=/your/work/directory/a3_hic_alignment/
            cd ${wd}
            cd hic_${sample}
            cd bwa_align
            cd lg_wise_contigs_read_align
            cd align_${lg}
            #
            for readset in ${sample}_L1 ${sample}_L2 ${sample}_L3; do
                filterBAM_forHiC.pl ${readset}_${lg}_bwa_aln.bam ${readset}_${lg}_bwa_aln_clean.sam
                samtools view -bt ../../lg_wise_contigs/clipped4_${sample}_${lg}.fa.fai ${readset}_${lg}_bwa_aln_clean.sam -o ${readset}_${lg}_bwa_aln_clean.bam
                samtools flagstat ${readset}_${lg}_bwa_aln_clean.bam > ${readset}_${lg}_bwa_aln_clean_flagstat.txt
                rm ${readset}_${lg}_bwa_aln_clean.sam ${readset}_${lg}_bwa_aln.bam
            done
        done
    done

##### step 5.6. merge clean bams (if multiple subsets of Hi-C reads are aligned, here we just use the one alignment as example)

    wd=/your/work/directory/a3_hic_alignment/
    cd ${wd}

    for sample in O; do
        #
        for lg in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12; do
            wd=/your/work/directory/a3_hic_alignment/
            cd ${wd}
            cd hic_${sample}
            cd bwa_align
            cd lg_wise_contigs_read_align
            cd align_${lg}
            s1=${sample}_L1_${lg}
            s2=${sample}_L2_${lg}
            s3=${sample}_L3_${lg}
            samtools view -H ${s1}_bwa_aln_clean.bam > ${s1}_header.sam
            samtools merge -nrf -h ${s1}_header.sam ${sample}_${lg}_L123_bwa_aln_clean.bam ${s1}_bwa_aln_clean.bam ${s2}_bwa_aln_clean.bam ${s3}_bwa_aln_clean.bam
            samtools flagstat ${sample}_${lg}_L123_bwa_aln_clean.bam > ${sample}_${lg}_L123_bwa_aln_clean_flagstat.txt
        done
    done

##### step 5.7. identify allelic contigs, see https://github.com/tangerzhang/ALLHiC/wiki/ALLHiC:-identify-allelic-contigs

    wd=/your/work/directory/a3_hic_alignment/
    cd ${wd}
    ln -s /path/to/src_bash_scripts/gmap.sh ln_gmap.sh

    for sample in O; do
        for lg in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12; do
            wd=/your/work/directory/a3_hic_alignment/
            cd ${wd}
            cd hic_${sample}
            cd bwa_align
            cd lg_wise_contigs_read_align
            mkdir z_allelic_map_${lg}
            cd z_allelic_map_${lg}
            genome=../../lg_wise_contigs/clipped4_${sample}_${lg}.fa
            ref_cds=/your/work/directory/ref_dm6p1/DM_gene_models_cdna_${lg}.fa
            chr_gff=/your/work/directory/ref_dm6p1/DM_gene_models_${lg}.gff3
            nthread=8
            # this will generate Allele.ctg.table in three steps wrapped in ln_gmap.sh
            ../../../../ln_gmap.sh ${genome} ${ref_cds} ${chr_gff} ${nthread}
        done
    done

##### step 5.8. use Allele.ctg.table to prune 1) signals that link alleles and 2) weak signals from BAM files: caution: ALLHiC_prune in the original version generatd Tb-level intermeidate files. Here we used a developing version of prune: https://github.com/sc-zhang/ALLHiC_components

    wd=/your/work/directory/a3_hic_alignment/
    cd ${wd}

    for sample in O; do
        for lg in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12; do
            wd=/your/work/directory/a3_hic_alignment/
            cd ${wd}
            cd hic_${sample}
            cd bwa_align
            cd lg_wise_contigs_read_align
            mkdir z_bam_prunning_${lg}
            cd z_bam_prunning_${lg}
            allele_table=../z_allelic_map_${lg}/Allele.ctg.table
            bam=../align_${lg}/${sample}_${lg}_L123_bwa_aln_clean.bam
            /path/to/ALLHiC_dev/ALLHiC_components/Prune/ALLHiC_prune -i ${allele_table} -b ${bam}
        done
    done

#### step 6. hic_binning to separate contigs into haplotype-specific groups - CORE FUNCTION

    sample="O"
    wd=/your/work/directory/a3_hic_alignment/
    cd ${wd}
    for sample in O; do
        for lg in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12; do
            wd=/your/work/directory/a3_hic_alignment/
            cd ${wd}
            cd hic_${sample}
            cd bwa_align
            cd lg_wise_contigs_read_align
            #
            mkdir hic_binning_${lg}
            cd hic_binning_${lg}
            # contigs cannot be linked according to allelic table
            allele_table=../z_allelic_map_${lg}/Allele.ctg.table
            # allhic-pruned hic bam according to allelic table
            bam=../z_bam_prunning_${lg}/prunning.bam
            # lg-wise contigs
            group_lg_path=/your/work/directory/a4_alignment_based_linkage_grouping/
            dm_res_grouping_details=${group_lg_path}/zhifiasm_v0p7_dm_${sample}_asm/dm_res_grouping_details_${sample}_${lg}.txt
            # alignment of contigs to DM reference => alignment-based allelic contigs
            ctg_alignment_to_dm=${group_lg_path}/zhifiasm_v0p7_dm_${sample}_asm/${sample}_against_dm.paf
            # window-level tig markers
            marker_path=/your/work/directory/win_marker/
            tig_win_marker=${marker_path}/${sample}_cnv_winsize10000_step10000_hq_merged_vs_hifi_markers_20221102_wsize500kb_final.txt
            # minimum hic contact to build up initial clusters of haplotype-wise contigs
            mkdir running_logs
            for min_hic_contact in 5 10 15 20 25 35 45 55 65 70 80; do
                # minimum size of haplotigs to build up initial clusters of haplotype-wise contigs
                for min_hapctg_size in 10000 50000 100000 150000 200000 250000 300000; do
                    for min_hic_contact_i in 0.55 1.1 2.2 3.3; do
                        for max_allelic_ratio in 0.01 0.02 0.03 0.04 0.05 0.10 0.15 0.20; do
                            samtools view -F 3840 ${bam} | hic_binning_vt --min-hic ${min_hic_contact} --min-hic-i ${min_hic_contact_i} --hap-win-r 1 --min-hap-size ${min_hapctg_size} --max-allelic-r ${max_allelic_ratio} --ctg ${dm_res_grouping_details} --win-marker ${tig_win_marker} --gmap-allelic ${allele_table} --align-paf ${ctg_alignment_to_dm} -o ${lg}_${min_hic_contact}_${min_hic_contact_i}_${min_hapctg_size}_${max_allelic_ratio} > ./running_logs/run_${lg}_${min_hic_contact}_${min_hic_contact_i}_${min_hapctg_size}_${max_allelic_ratio}.log
                        done
                    done
                done
            done
        done
    done

#### step 7. check the binning result

##### step 7.1. get group-contig sizes of chromosomes => s8_grouping_window_markers_refined_1st.txt

    sample="O"
    wd=/your/work/directory/a3_hic_alignment/
    cd ${wd}
    for sample in O; do
        wd=/your/work/directory/a3_hic_alignment/
        cd ${wd}
        cd hic_${sample}
        cd bwa_align
        cd lg_wise_contigs_read_align
        >interdiate_res_sample_${sample}_hic_binining_marker_size.txt

        for lg in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12; do
            wd=/your/work/directory/a3_hic_alignment/
            cd ${wd}
            cd hic_${sample}
            cd bwa_align
            cd lg_wise_contigs_read_align
            cd hic_binning_${lg}
            for min_hic_contact_i in 0.55 1.1 2.2 3.3; do
                for min_hic_contact in 5 10 15 20 25 35 45 55 65 70 80; do
                    for min_hapctg_size in 10000 50000 100000 150000 200000 250000 300000; do
                        for max_allelic_ratio in 0.01 0.02 0.03 0.04 0.05 0.10 0.15 0.20; do
                            echo ${sample}":"${lg}_${min_hic_contact}_${min_hic_contact_i}_${min_hapctg_size}_${max_allelic_ratio}
                            echo ${lg}_${min_hic_contact}_${min_hic_contact_i}_${min_hapctg_size}_${max_allelic_ratio} >> ../interdiate_res_sample_${sample}_hic_binining_marker_size.txt
                            for i in 1 2 3 4; do
                                awk -v hapi=$i '$5==hapi' ./s6_${lg}_${min_hic_contact}_${min_hic_contact_i}_${min_hapctg_size}_${max_allelic_ratio}_raw_tig_marker_cross_link_count/s8_grouping_window_markers_refined_1st.txt | awk '{s+=$3-$2+1} END {print s}' >> ../interdiate_res_sample_${sample}_hic_binining_marker_size.txt;
                            done
                        done
                    done
                done
            done
        done

        # reformatting: some binning with linkage number < 4 will be removed from below process
        wd=/your/work/directory/a3_hic_alignment/
        cd ${wd}
        cd hic_${sample}
        cd bwa_align
        cd lg_wise_contigs_read_align
        sed ':a;N;$!ba;s/\n/\t/g' interdiate_res_sample_${sample}_hic_binining_marker_size.txt | sed 's/chr/\nchr/g' | sed 's/_/\t/g' | grep 'chr' | awk '$9!=""' > interdiate_res_sample_${sample}_hic_binining_marker_size_reformat.txt
    done

##### step 7.2. select the best representative phasing for each chr.

		# sort by allelic - smaller better (= lower chance to mis-join contigs into linkage groups)
    #		    tmp_x <- data_chr_sorted_cleaned[order(data_chr_sorted_cleaned$V5, decreasing=FALSE), ]
    #		    tmp_x <- tmp_x[1:min(150, length(tmp_x$V1)), ]
		# sort by hic-strength - larger better - backbone construction
    #		    tmp_x <- tmp_x[order(tmp_x$V2, decreasing=TRUE), ]
    #		    tmp_x <- tmp_x[1:min(75, length(tmp_x$V1)), ]
		# sort by hic-strength - larger better - integrate
    #		    tmp_x <- tmp_x[order(tmp_x$V3, decreasing=TRUE), ]
    #		    tmp_x <- tmp_x[1:min(30, length(tmp_x$V1)), ]
		# sort by intra-group hic / inter-group hic - larger better
    #		    tmp_x <- tmp_x[order(tmp_x$V15, decreasing=TRUE), ]

    # you need to update path in the R script
    Rscript /path/to/tetraDecoder/src_hb_stage1_hic_binning/aux/z_suppl_sample_10cultivars_checking_after_hic_binining.R

##### step 7.3. merge phased markers within different chrs as a genome-level marker set    
    # caution: check max allelic ratio in best cases file must be two digits like 0.xx (, even for 0.10 or 0.20)    
    # example: lg=chr01; min_hic_contact=45; min_hic_contact_i=3.3; min_hapctg_size=200000; max_allelic_ratio=0.03

    wd=/your/work/directory/a3_hic_alignment/
    cd ${wd}
    for sample in O; do
        echo ${sample}
        wd=/your/work/directory/a3_hic_alignment/
        cd ${wd}
        cd hic_${sample}
        cd bwa_align
        cd lg_wise_contigs_read_align
        chri=0
        >final_res_${sample}_window_markers_12chrs.txt
        while read lg min_hic_contact min_hic_contact_i min_hapctg_size max_allelic_ratio hs1 hs2 hs3 hs4 hic_r size_sd; do
            phased_win_marker=./hic_binning_${lg}/s6_${lg}_${min_hic_contact}_${min_hic_contact_i}_${min_hapctg_size}_${max_allelic_ratio}_raw_tig_marker_cross_link_count/s8_grouping_window_markers_refined_1st.txt
            chri=$((chri+1))
            awk '$5==-1' ${phased_win_marker} > tmp_non_grouped_markers.txt
            awk '$5!=-1' ${phased_win_marker} > tmp_all_grouped_markers.txt
            #
            awk -v i="$chri" '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,$4,$5+(i-1)*4,$6,$7,$8,$9)}' tmp_all_grouped_markers.txt | sed "s/${lg}_${min_hic_contact}_${min_hic_contact_i}_${min_hapctg_size}_${max_allelic_ratio}/${lg}\t${min_hic_contact}_${min_hic_contact_i}_${min_hapctg_size}_${max_allelic_ratio}/g" >> final_res_${sample}_window_markers_12chrs.txt
            sed "s/${lg}_${min_hic_contact}_${min_hic_contact_i}_${min_hapctg_size}_${max_allelic_ratio}/-1\t${lg}_${min_hic_contact}_${min_hic_contact_i}_${min_hapctg_size}_${max_allelic_ratio}/g" tmp_non_grouped_markers.txt | sort -k1,1 -k2,2n >> final_res_${sample}_window_markers_12chrs.txt
            rm tmp_non_grouped_markers.txt tmp_all_grouped_markers.txt
            #
        done < ./final_res_${sample}_binning_best_cases.txt
        # calculate total marker size
        awk '{s+=$3-$2+1} END {print s}' final_res_${sample}_window_markers_12chrs.txt
        awk '$5!=-1' final_res_${sample}_window_markers_12chrs.txt | awk '{s+=$3-$2+1} END {print s}'
    done

#### step 8. extract HiFi reads according to marker grouping - group HiFi reads to 1-48 linkage groups according to marker phasing/grouping.

    sample="O"
    wd=/your/work/directory/a3_hic_alignment/
    cd ${wd}
    for sample in O; do
        wd=/your/work/directory/a3_hic_alignment/
        cd ${wd}
        cd hic_${sample}
        cd bwa_align
        cd lg_wise_contigs_read_align
        mkdir s_asm_phased_reads
        cd s_asm_phased_reads
        bam=/your/work/directory/assembly/${sample}_hifiasm.bam
        phased_win_marker=../final_res_${sample}_window_markers_12chrs.txt
        samtools view ${bam} | long_read_separator_v2 - ${phased_win_marker} hifi_lgs > hifi_separation.log
        tail -n 64 hifi_separation.log > hifi_separation_simple.log
        cd ./hifi_lgs_window_marker_separated_reads
        gzip *.fa
    done

#### step 9. LG-wise assembly

    wd=/your/work/directory/a3_hic_alignment/
    cd ${wd}

    for sample in O; do
        chri=0
        for lg in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12; do
            wd=/your/work/directory/a3_hic_alignment/
            cd ${wd}
            cd hic_${sample}
            cd bwa_align
            cd lg_wise_contigs_read_align
            cd s_asm_phased_reads
            chri=$((chri+1))
            for hapi in 1 2 3 4; do
                real_hapi=$(((chri-1)*4+hapi))
                mkdir zasm_homLG_${lg}_LG_${real_hapi}
                cd zasm_homLG_${lg}_LG_${real_hapi}
                ln -s /path/to/v0.7/hifiasm hifiasm
                hifi=../hifi_lgs_window_marker_separated_reads/homLG_${lg}_LG_${real_hapi}_reads.fa.gz
                ./hifiasm -t 4 -o homLG_${lg}_LG_${real_hapi}_asm ${hifi}; rm *.bin
                cd ..
            done
        done
    done

#### step 10. contig scaffolding can be done with RagTag/Ragoo, ALLHiC, Juicer/Juicebox, which is not covered here.
