# Here we develop the new pipeline for hic-based phasing and linkage groupping in tetaploids => haplotype-resolved tetraploid genome assembly

cd /your/path//a3_hic_alignment/bgi_deep_seq/

# Note, for hic alignment to clipped4_${sample}_hifiasm.p_utg.gfa.fa and separation of reads into lgs.fastq.gz: ./s2_a3_hic_alignment_allhic.sh

# step 1. extract lg-wise contigs (groupped based on DM reference)
#
wd=/your/path//a3_hic_alignment/bgi_deep_seq/
cd ${wd}
for sample in A B C D E F G H I J O; do
    #
    wd=/your/path//a3_hic_alignment/bgi_deep_seq/
    cd ${wd}
    cd hic_${sample}
    cd bwa_align
    mkdir lg_wise_contigs
    cd lg_wise_contigs
    genome=../../../../../a2_initial_assembly/sample_${sample}/hifiasm_asm_v0p7/subset4_illu_re_align_purge_ovl/DNB_align/clipped4_${sample}_hifiasm.p_utg.gfa.fa
    ls -l ${genome}
    bsub -q short -o fasta_name_selecter.log -e fasta_name_selecter.err "for lg in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12; do cut -f 1 ../../../../../a4_alignment_based_linkage_grouping/zhifiasm_v0p7_dm_${sample}_asm/dm_res_grouping_details_${sample}_\${lg}.txt > dm_res_grouping_details_${sample}_\${lg}_ids.txt; fasta_name_selecter ${genome} dm_res_grouping_details_${sample}_\${lg}_ids.txt; mv clipped4_${sample}_hifiasm.p_utg.gfa_selected.fa clipped4_${sample}_\${lg}.fa; done"
   #
done

# step 2. index lg-wise ctgs for BWA,bowtie2

wd=/your/path//a3_hic_alignment/bgi_deep_seq/
cd ${wd}
for sample in A B C D E F G H I J O; do
    wd=/your/path//a3_hic_alignment/bgi_deep_seq/
    cd ${wd}
    cd hic_${sample}
    cd bwa_align
    cd lg_wise_contigs
    for lg in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12; do
        ll  clipped4_${sample}_${lg}.fa*
        bsub -o bwa_index.log -q short -R "rusage[mem=10000]" -R "span[hosts=1]" -M 10000 "bwa index -a bwtsw clipped4_${sample}_${lg}.fa"
        bsub -o samtools_index.log -q short -R "rusage[mem=10000]" -R "span[hosts=1]" -M 10000 "samtools faidx clipped4_${sample}_${lg}.fa"
        bsub -o bowtie2_index.log -q ioheavy -n 8 -R "rusage[mem=10000]" -R "span[hosts=1]" -M 10000 "bowtie2-build --threads 8 clipped4_${sample}_${lg}.fa clipped4_${sample}_${lg}.fa"
    done
    cd ${wd}
done

# step 3. bwa align

wd=/your/path//a3_hic_alignment/bgi_deep_seq/
cd ${wd}
for sample in A B C D E F G H I J O; do
    #
    wd=/your/path//a3_hic_alignment/bgi_deep_seq/
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
            #
            ll ${genome} ${R1} ${R2}
            #
            bsub -o bwa_aln.log -e bwa_aln.err -q multicore20 -n 8 -R "rusage[mem=8000]" -R "span[hosts=1]" -M 8000 "bwa aln -t 8 ${genome} ${R1} > ${readset}_${lg}_R1.sai"
            bsub -o bwa_aln.log -e bwa_aln.err -q multicore20 -n 8 -R "rusage[mem=8000]" -R "span[hosts=1]" -M 8000 "bwa aln -t 8 ${genome} ${R2} > ${readset}_${lg}_R2.sai"
            #
            cd ..
        done
    done
   #
done

# step 4. merge read pairs

wd=/your/path//a3_hic_alignment/bgi_deep_seq/
cd ${wd}
for sample in A B C D E F G H I J O; do
    for lg in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12; do
        for readset in ${sample}_L1 ${sample}_L2 ${sample}_L3; do
            wd=/your/path//a3_hic_alignment/bgi_deep_seq/
            cd ${wd}
            cd hic_${sample}
            cd bwa_align
            cd lg_wise_contigs_read_align
            cd align_${lg}
            genome=../../lg_wise_contigs/clipped4_${sample}_${lg}.fa
            R1="../../${readset}_${lg}_extract_R1.fastq.gz"
            R2="../../${readset}_${lg}_extract_R2.fastq.gz"
            #
            ll ${genome} ${R1} ${R2} *.sai
            #
            bsub -o bwa_sampe.log -e bwa_sampe.err -q normal -R "rusage[mem=5000]" -R "span[hosts=1]" -M 5000 "bwa sampe ${genome} ${readset}_${lg}_R1.sai ${readset}_${lg}_R2.sai ${R1} ${R2} | samtools view -b -F12 -o ${readset}_${lg}_bwa_aln.bam; samtools flagstat ${readset}_${lg}_bwa_aln.bam > ${readset}_${lg}_bwa_aln_bam_flagstat.txt"
        done
    done
   #
done

# step 5. filter reads in bam

wd=/your/path//a3_hic_alignment/bgi_deep_seq/
cd ${wd}
for sample in A B C D E F G H I J O; do
    for lg in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12; do
        wd=/your/path//a3_hic_alignment/bgi_deep_seq/
        cd ${wd}
        cd hic_${sample}
        cd bwa_align
        cd lg_wise_contigs_read_align
        cd align_${lg}
        #
        for readset in ${sample}_L1 ${sample}_L2 ${sample}_L3; do
            bsub -o Filtering.log -e Filtering.err -q normal -R "rusage[mem=1000]" -R "span[hosts=1]" -M 1000 "filterBAM_forHiC.pl ${readset}_${lg}_bwa_aln.bam ${readset}_${lg}_bwa_aln_clean.sam; samtools view -bt ../../lg_wise_contigs/clipped4_${sample}_${lg}.fa.fai ${readset}_${lg}_bwa_aln_clean.sam -o ${readset}_${lg}_bwa_aln_clean.bam; samtools flagstat ${readset}_${lg}_bwa_aln_clean.bam > ${readset}_${lg}_bwa_aln_clean_flagstat.txt; rm ${readset}_${lg}_bwa_aln_clean.sam ${readset}_${lg}_bwa_aln.bam"
        done
    done
done
#
# step 6. merge clean bams:
#
#         note after merging and getting file, if everything works well, three types of intermediate can be removed (to save storage):
#              ${wd}/hic_*/bwa_align/*_L*_chr*_extract_R*.fastq.gz
#              ${wd}/hic_*/bwa_align/lg_wise_contigs_read_align/align_chr*/*.sai
#              ${wd}/hic_*/bwa_align/lg_wise_contigs_read_align/align_chr*/*_L*_chr*_bwa_aln_clean.bam
#         result:
#              ${wd}/hic_*/bwa_align/lg_wise_contigs_read_align/align_chr*/*_chr*_L123_bwa_aln_clean.bam
#
wd=/your/path//a3_hic_alignment/bgi_deep_seq/
cd ${wd}
for sample in A B C D E F G H I J O; do
    #
    for lg in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12; do
        wd=/your/path//a3_hic_alignment/bgi_deep_seq/
        cd ${wd}
        cd hic_${sample}
        cd bwa_align
        cd lg_wise_contigs_read_align
        cd align_${lg}
        #
        s1=${sample}_L1_${lg}
        s2=${sample}_L2_${lg}
        s3=${sample}_L3_${lg}
        #
        bsub -o merge_bam.log -e merge_bam.err -q normal -R "rusage[mem=1000]" -R "span[hosts=1]" -M 1000 "samtools view -H ${s1}_bwa_aln_clean.bam > ${s1}_header.sam; samtools merge -nrf -h ${s1}_header.sam ${sample}_${lg}_L123_bwa_aln_clean.bam ${s1}_bwa_aln_clean.bam ${s2}_bwa_aln_clean.bam ${s3}_bwa_aln_clean.bam; samtools flagstat ${sample}_${lg}_L123_bwa_aln_clean.bam > ${sample}_${lg}_L123_bwa_aln_clean_flagstat.txt"
    done
    #
done

# step 7. identify allelic contigs, see https://github.com/tangerzhang/ALLHiC/wiki/ALLHiC:-identify-allelic-contigs
#         see here for pre-process of gff and cds files: /Potato_single_cells/reference_DH_line_2020/reference_sequence/protein_sequences/separation_to_lgs_for_multicultivar_analysis.sh
          # step 1: ~5 mins
          # step 2: note: ~30 mins
          #         1) target.genome is the contig level of polyploid genome assembly
          #         2) $N could be the ploidy of your target genome, for example $N=4 if it is a tetraploid
          #         3) ref_cds is the coding sequences of diploid genome, which could be reference to get allelic table
          # step 3: note: 1 second
          #         1) ref_gff3 is the gff3 annotation of the diploid genome
          #         2) gmap2Alleletable.pl could be get in the following link: https://github.com/tangerzhang/ALLHiC/blob/master/scripts/gmap2AlleleTable.pl
          #         3) format of Allele.ctg.table: chromosomeID position contig1 [contig2 contig 3 contig4].
          #                                        Prune step will remove the Hi-C linked reads between allelic contigs.
          # caution: I updated gmap2AlleleTableBED2 to get the expected info: chr pos ctg1[ ctg2 ctg3 ctg4]
          # gmap.gff3 generated above is a hardcoded input of gmap2AlleleTableBED2.pl
#
wd=/your/path//a3_hic_alignment/bgi_deep_seq/
cd ${wd}
ln -s /biodata/dep_mercier/grp_schneeberger/projects/methods/src_shq/z_GameteBinning_potato_multivar/src_bash_scripts/gmap.sh ln_gmap.sh
for sample in A B C D E F G H I J O; do
    for lg in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12; do
        wd=/your/path//a3_hic_alignment/bgi_deep_seq/
        cd ${wd}
        cd hic_${sample}
        cd bwa_align
        cd lg_wise_contigs_read_align
        mkdir z_allelic_map_${lg}
        cd z_allelic_map_${lg}
        #
        genome=../../lg_wise_contigs/clipped4_${sample}_${lg}.fa
        ref_cds=../../../../../../../../Potato_single_cells/reference_DH_line_2020/reference_sequence/protein_sequences/DM_gene_models_cdna_${lg}.fa
        chr_gff=../../../../../../../../Potato_single_cells/reference_DH_line_2020/reference_sequence/protein_sequences/DM_gene_models_${lg}.gff3
        nthread=8
        # this will generate Allele.ctg.table in three steps wrapped in ln_gmap.sh
        bsub -o ln_gmap.log -e ln_gmap.err -q normal -R "rusage[mem=10000]" -R "span[hosts=1]" -M 10000 "../../../../ln_gmap.sh ${genome} ${ref_cds} ${chr_gff} ${nthread}"
        #
        cd ..
    done
    #
done

# step 8. use Allele.ctg.table to prune 1) signals that link alleles and 2) weak signals from BAM files: caution: ALLHiC_prune in the original version generatd Tb-level intermeidate files. Here is a developing version of prune used: https://github.com/sc-zhang/ALLHiC_components
#
wd=/your/path//a3_hic_alignment/bgi_deep_seq/
cd ${wd}
for sample in A B C D E F G H I J O; do
    for lg in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12; do
        wd=/your/path//a3_hic_alignment/bgi_deep_seq/
        cd ${wd}
        cd hic_${sample}
        cd bwa_align
        cd lg_wise_contigs_read_align
        mkdir z_bam_prunning_${lg}
        cd z_bam_prunning_${lg}
        #
        allele_table=../z_allelic_map_${lg}/Allele.ctg.table
        bam=../align_${lg}/${sample}_${lg}_L123_bwa_aln_clean.bam
        bsub -o ALLHiC_prune.log -e ALLHiC_prune.err -q normal -R "rusage[mem=1000]" -R "span[hosts=1]" -M 1000 "/netscratch/dep_mercier/grp_schneeberger/bin/ALLHiC_dev/ALLHiC_components/Prune/ALLHiC_prune -i ${allele_table} -b ${bam}"
        #/netscratch/dep_mercier/grp_schneeberger/bin/ALLHiC_dev/ALLHiC_components/Prune/ALLHiC_prune -i ${allele_table} -b ${bam} &>ALLHiC_prune.log&
        #
        cd ..
    done
   #
done

########################################################################################################################
# step 9. now I use hic_binning to separate contigs into haplotype-specific groups

wd=/your/path//a3_hic_alignment/bgi_deep_seq/
cd ${wd}
for sample in A B C D E F G H I J O; do
    for lg in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12; do
        wd=/your/path//a3_hic_alignment/bgi_deep_seq/
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
        dm_res_grouping_details=../../../../../../a4_alignment_based_linkage_grouping/zhifiasm_v0p7_dm_${sample}_asm/dm_res_grouping_details_${sample}_${lg}.txt
        # alignment of contigs to DM reference => alignment-based allelic contigs
        ctg_alignment_to_dm=../../../../../../a4_alignment_based_linkage_grouping/zhifiasm_v0p7_dm_${sample}_asm/${sample}_against_dm.paf
        # window-level tig markers
        tig_win_marker=../../../../../../../s2_10_Cultivars_PacBio_HiFi/a2_initial_assembly/sample_${sample}/hifiasm_asm_v0p7/subset4_illu_re_align_purge_ovl/${sample}_cnv_winsize10000_step10000_hq_merged_vs_hifi_markers_20221102_wsize500kb_final.txt
        # minimum hic contact to build up initial clusters of haplotype-wise contigs - TODO we need a statistical analysis to determine value like '70' here!
        mkdir running_logs
        for min_hic_contact in 5 10 15 20 25 35 45 55 65 70 80; do
        # minimum size of haplotigs to build up initial clusters of haplotype-wise contigs
        for min_hapctg_size in 10000 50000 100000 150000 200000 250000 300000; do
        bsub -m "hpc001 hpc002 hpc003 hpc004 hpc005 hpc006" -o hic_binning.log -e hic_binning.err -q normal -R "rusage[mem=1000]" -R "span[hosts=1]" -M 2000 "for min_hic_contact_i in 0.55 1.1 2.2 3.3; do for max_allelic_ratio in 0.01 0.02 0.03 0.04 0.05 0.10 0.15 0.20; do samtools view -F 3840 ${bam} | hic_binning_vt --min-hic ${min_hic_contact} --min-hic-i \${min_hic_contact_i} --hap-win-r 1 --min-hap-size ${min_hapctg_size} --max-allelic-r \${max_allelic_ratio} --ctg ${dm_res_grouping_details} --win-marker ${tig_win_marker} --gmap-allelic ${allele_table} --align-paf ${ctg_alignment_to_dm} -o ${lg}_${min_hic_contact}_\${min_hic_contact_i}_${min_hapctg_size}_\${max_allelic_ratio} > ./running_logs/run_${lg}_${min_hic_contact}_\${min_hic_contact_i}_${min_hapctg_size}_\${max_allelic_ratio}.log; done; done"
        done
        done
        #
        cd ..
    done
   #
done

# step 10. select the binning result: some binning clearly failed, some successful. Which one to use?
#          According to Otava,
#          if I select in the following way,
#          I got an overal phasing accuracy of 98.8% - though not every chr was good. - Can I further improve phasing?

wd=/your/path//a3_hic_alignment/bgi_deep_seq/
cd ${wd}
for sample in A B C D E F G H I J O; do
    wd=/your/path//a3_hic_alignment/bgi_deep_seq/
    cd ${wd}
    cd hic_${sample}
    cd bwa_align
    cd lg_wise_contigs_read_align
    >interdiate_res_sample_${sample}_hic_binining_marker_size.txt
    for lg in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12; do
        #
        wd=/your/path//a3_hic_alignment/bgi_deep_seq/
        cd ${wd}
        cd hic_${sample}
        cd bwa_align
        cd lg_wise_contigs_read_align
        cd hic_binning_${lg}
        #
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
        #
    done
    #
    # reformatting: some binning with linkage number < 4 will be removed from below process
    wd=/your/path//a3_hic_alignment/bgi_deep_seq/
    cd ${wd}
    cd hic_${sample}
    cd bwa_align
    cd lg_wise_contigs_read_align
    sed ':a;N;$!ba;s/\n/\t/g' interdiate_res_sample_${sample}_hic_binining_marker_size.txt | sed 's/chr/\nchr/g' | sed 's/_/\t/g' | grep 'chr' | awk '$9!=""' > interdiate_res_sample_${sample}_hic_binining_marker_size_reformat.txt
   #
done

# step 11. select the best representative phasing for each chr for each cultivar.
#          note 20221229, selection criteria:
#             11.1 C:chr10,E:chr06,G:chr03: if(this_diff > 30000000 | data_chr_sorted[ci,16] > 18000000) #
#             11.2 all other chrs.........: if(this_diff > 20000000 | data_chr_sorted[ci,16] > 15000000)
#             and
#		# sort by allelic - smaller better (= lower chance to mis-join contigs into linkage groups)
#		    tmp_x <- data_chr_sorted_cleaned[order(data_chr_sorted_cleaned$V5, decreasing=FALSE), ]
#		    tmp_x <- tmp_x[1:min(150, length(tmp_x$V1)), ]
#		# sort by hic-strength - larger better - backbone construction
#		    tmp_x <- tmp_x[order(tmp_x$V2, decreasing=TRUE), ]
#		    tmp_x <- tmp_x[1:min(75, length(tmp_x$V1)), ]
#		# sort by hic-strength - larger better - integrate
#		    tmp_x <- tmp_x[order(tmp_x$V3, decreasing=TRUE), ]
#		    tmp_x <- tmp_x[1:min(30, length(tmp_x$V1)), ]
#		# sort by intra-group hic / inter-group hic - larger better
#		    tmp_x <- tmp_x[order(tmp_x$V15, decreasing=TRUE), ]
#		# select top 1
#		    min_case  <- find_min( tmp_x[1:min(1, length(tmp_x$V1)), c(1:9, 15, 16)] )
# Corrected cases: sort tmp_x with hap_size_sd etc - see z_suppl_sample_10cultivars_checking_after_hic_binining_manual_correction.R
# note, contigs require further correction according to Hi-C signal between hap-chrs after chr-level scaffolding: a ctg in hap-a but showing higher hic to hap-b, should be re-assigned.

# update path in the script before running, and let's discuss how to do the selection once you have the phasing result.

Rscript ${scripath}/z_suppl_sample_10cultivars_checking_after_hic_binining.R                   # for general selection for A,B,D,F,H,I,J

# step 12. merge phased markers within different chrs as a genome-level marker set
#          caution: check max allelic ratio in best cases file must be two digits like 0.xx (, even for 0.10 or 0.20)
#          example: lg=chr01; min_hic_contact=45; min_hic_contact_i=3.3; min_hapctg_size=200000; max_allelic_ratio=0.03

wd=/your/path//a3_hic_alignment/bgi_deep_seq/
cd ${wd}
for sample in A B C D E F G H I J O; do
    echo ${sample}
    wd=/your/path//a3_hic_alignment/bgi_deep_seq/
    cd ${wd}
    cd hic_${sample}
    cd bwa_align
    cd lg_wise_contigs_read_align
    #
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
        #
        rm tmp_non_grouped_markers.txt tmp_all_grouped_markers.txt
        #
    done < ./final_res_${sample}_binning_best_cases.txt
    # calculate total marker size
    awk '{s+=$3-$2+1} END {print s}' final_res_${sample}_window_markers_12chrs.txt
    awk '$5!=-1' final_res_${sample}_window_markers_12chrs.txt | awk '{s+=$3-$2+1} END {print s}'
   #
done

# step 13. extract HiFi reads according to marker grouping
#          code-not-used: bsub -o hifi_read_extraction.log -e hifi_read_extraction.err -q normal -R "rusage[mem=8000]" -R "span[hosts=1]" -M 8000 "samtools view ${cram} | hic_long_read_separator - ${phased_win_marker} ${lg} > hifi_read_extraction_detail.log; gzip *.fa" # this only works for individual chrs
#          C: has a lot reads 7.43 Gb / 87.3 Gb missed, what are they? Need to assemble them.

wd=/your/path//a3_hic_alignment/bgi_deep_seq/
cd ${wd}
for sample in A B C D E F G H I J O; do
    wd=/your/path//a3_hic_alignment/bgi_deep_seq/
    cd ${wd}
    cd hic_${sample}
    cd bwa_align
    cd lg_wise_contigs_read_align
    mkdir s_asm_phased_reads
    cd s_asm_phased_reads
    cram=../../../../../../a2_initial_assembly/sample_${sample}/hifiasm_asm_v0p7/subset4_illu_re_align_purge_ovl/HIFI_align/ptt${sample}_hiasm_ref_pilon_subset4.cram
    phased_win_marker=../final_res_${sample}_window_markers_12chrs.txt
    ls -l ${cram} ${phased_win_marker}
    bsub -m "hpc001 hpc002 hpc003 hpc004 hpc005 hpc006" -o hifi_read_extraction.log -e hifi_read_extraction.err -q normal -R "rusage[mem=8000]" -R "span[hosts=1]" -M 8000 "samtools view ${cram} | long_read_separator_v2 - ${phased_win_marker} hifi_lgs > hifi_separation.log; tail -n 64 hifi_separation.log > hifi_separation_simple.log; cd ./hifi_lgs_window_marker_separated_reads; gzip *.fa"
   #
done

# step 14. lg-wise re-assembly

wd=/your/path//a3_hic_alignment/bgi_deep_seq/
cd ${wd}
for sample in A B C D E F G H I J O; do
    chri=0
    for lg in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12; do
        wd=/your/path//a3_hic_alignment/bgi_deep_seq/
        cd ${wd}
        cd hic_${sample}
        cd bwa_align
        cd lg_wise_contigs_read_align
        cd s_asm_phased_reads
        #
        chri=$((chri+1))
        #
        for hapi in 1 2 3 4; do
            #
            real_hapi=$(((chri-1)*4+hapi))
            mkdir zasm_homLG_${lg}_LG_${real_hapi}
            cd zasm_homLG_${lg}_LG_${real_hapi}
            #
            ln -s /netscratch/dep_mercier/grp_schneeberger/projects/Method/gamete_binning_book_chapter_v9_ms/software/bin/hifiasm hifiasm
            hifi=../hifi_lgs_window_marker_separated_reads/homLG_${lg}_LG_${real_hapi}_reads.fa.gz
            ll ${hifi}
            # hifiasm v0p7, with fewer redundant unitigs and higher N50 (than version 0.16) - finally used 20220921
            bsub -m "hpc001 hpc002 hpc003 hpc004 hpc005 hpc006" -q ioheavy -R "rusage[mem=24000]" -M 24000 -n 4 -o homLG_${lg}_LG_${real_hapi}.log -e homLG_${lg}_LG_${real_hapi}.err "./hifiasm -t 4 -o homLG_${lg}_LG_${real_hapi}_asm ${hifi}; rm *.bin"
            #
            cd ..
        done
        #
    done
   #
done

# step 15. get assembly statistics
           # cat homLG_${lg}_LG_${real_hapi}_asm.p_ctg.gfa | grep '^S' | cut -f2,3 | awk '{print ">"$1"\n"$2}' > homLG_${lg}_LG_${real_hapi}_asm.p_ctg.fasta
           # fasta_length homLG_${lg}_LG_${real_hapi}_asm.p_ctg.fasta | grep '>' | sed 's/>//g' > homLG_${lg}_LG_${real_hapi}_asm.p_ctg.ctgsizes
           # /srv/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/bin/GitHub-schneeberger-group-jac/toolbox/Analysis/Assembly/calc_CN50.pl homLG_${lg}_LG_${real_hapi}_asm.p_ctg.fasta 110000000 1 &> homLG_${lg}_LG_${real_hapi}_asm.p_ctg.N50.calc.result.txt&

wd=/your/path//a3_hic_alignment/bgi_deep_seq/
cd ${wd}
for sample in A B C D E F G H I J O; do
    chri=0
    for lg in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12; do
        wd=/your/path//a3_hic_alignment/bgi_deep_seq/
        cd ${wd}
        cd hic_${sample}
        cd bwa_align
        cd lg_wise_contigs_read_align
        cd s_asm_phased_reads
        #
        chri=$((chri+1))
        #
        for hapi in 1 2 3 4; do
            #
            real_hapi=$(((chri-1)*4+hapi))
            cd zasm_homLG_${lg}_LG_${real_hapi}
            # Get Assembly stats
            bsub -m "hpc001 hpc002 hpc003 hpc004 hpc005 hpc006" -q short -R "rusage[mem=2000]" -M 4000 -o get_stats.log -e get_stats.err "cat homLG_${lg}_LG_${real_hapi}_asm.p_ctg.gfa | grep '^S' | cut -f2,3 | awk '{print \">\"\$1\"\\n\"\$2}' > homLG_${lg}_LG_${real_hapi}_asm.p_ctg.fasta; fasta_length homLG_${lg}_LG_${real_hapi}_asm.p_ctg.fasta | grep '>' | sed 's/>//g' > homLG_${lg}_LG_${real_hapi}_asm.p_ctg.ctgsizes; /srv/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/bin/GitHub-schneeberger-group-jac/toolbox/Analysis/Assembly/calc_CN50.pl homLG_${lg}_LG_${real_hapi}_asm.p_ctg.fasta 110000000 1 > homLG_${lg}_LG_${real_hapi}_asm.p_ctg.N50.calc.result.txt"
            #
            cd ..
        done
        #
    done
   #
done


# step 16. hic scaffolding: see ./a6_lg_wise_scaffolding/
#                               a6p1_filter_short_organelle_contigs.sh
#                               a6p2_lg_wise_scaffolding_omnic.sh
#                               a6p3_lg_wise_scaffolding_correction_with_alignment_to_DM_ref.sh
#                               a6p4_lg_wise_2nd_scaffolding_omnic.sh
#                               a6p5_lg_wise_scaffolding_alignment_to_DM_ref_and_visualization.sh


































# get number of hic reads

cd /your/path//a3_hic_alignment/bgi_deep_seq

>hic_read_number.txt
for sample in A B C D E F G H I J O; do
   grep 'reads; of these:' hic_${sample}/bwa_align/bowtie2.err | sed "s/reads; of these:/${sample}/g" >> hic_read_number.txt
done

cd /Users/sun/Desktop/z_10potato_project/xjtu_a3_hic_alignment
scp sun@dell-node-1.mpipz.mpg.de:/your/path//a3_hic_alignment/bgi_deep_seq/hic_read_number.txt ./
