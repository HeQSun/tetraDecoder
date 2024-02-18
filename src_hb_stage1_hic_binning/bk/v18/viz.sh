cd /biodata/dep_mercier/grp_schneeberger/projects/methods/src_shq/z_GameteBinning_potato_multivar/src_hb_stage1_hic_binning/

cd s6_pruned_bam_raw_tig_marker_cross_link_count/

fdp -Tpdf s6_haplotig_hic_contact_matrix_subclusters_raw.dot > s6_haplotig_hic_contact_matrix_subclusters_raw.pdf
fdp -Tpdf s6_haplotig_hic_contact_matrix_subclusters_raw_merged.dot > s6_haplotig_hic_contact_matrix_subclusters_raw_merged.pdf
fdp -Tpdf s7_haplotig_hic_contact_matrix_subclusters_hap_extended.dot > s7_haplotig_hic_contact_matrix_subclusters_hap_extended.pdf

dot_comparator ../gamete_binning_result_s3_refLG1_haplotig_GT_similarity_matrix_subclusters_final.dot s6_haplotig_hic_contact_matrix_subclusters_raw_merged.dot
dot_comparator ../gamete_binning_result_s3_refLG1_haplotig_GT_similarity_matrix_subclusters_final.dot s7_haplotig_hic_contact_matrix_subclusters_hap_extended.dot

dot_evaluator ../ln_O_s4p6_refine_grouping_final_window_markers_sorted.txt s6_haplotig_hic_contact_matrix_subclusters_raw.dot
dot_evaluator ../ln_O_s4p6_refine_grouping_final_window_markers_sorted.txt s7_haplotig_hic_contact_matrix_subclusters_hap_extended.dot


# working case with all default in coding 20221101: chr08 and chr10 are extreme cases with one haplotype with nearly no haplotigs.

samtools view -F 3840 ./data/ln_O_L123_prunning_chr01.bam | ./hic_binning - 70 ./data/dm_res_grouping_details_O_chr01.txt ./data/O_cnv_winsize10000_wsize500kb_final.txt 300000 ./data/ln_Allele.ctg.chr01.table ./data/ln_O_against_dm.paf 1 CHR1 > torm

samtools view -F 3840 ./data/ln_O_L123_prunning_chr02.bam | ./hic_binning - 70 ./data/dm_res_grouping_details_O_chr02.txt ./data/O_cnv_winsize10000_wsize500kb_final.txt 10000 ./data/ln_Allele.ctg.chr02.table ./data/ln_O_against_dm.paf 1 CHR2 > torm_CHR2

samtools view -F 3840 ./data/ln_O_L123_prunning_chr03.bam | ./hic_binning - 70 ./data/dm_res_grouping_details_O_chr03.txt ./data/O_cnv_winsize10000_wsize500kb_final.txt 10000 ./data/ln_Allele.ctg.chr03.table ./data/ln_O_against_dm.paf 1 CHR3 > torm_CHR3

samtools view -F 3840 ./data/ln_O_L123_prunning_chr04.bam | ./hic_binning - 70 ./data/dm_res_grouping_details_O_chr04.txt ./data/O_cnv_winsize10000_wsize500kb_final.txt 10000 ./data/ln_Allele.ctg.chr04.table ./data/ln_O_against_dm.paf 1 CHR4 > torm_CHR4

samtools view -F 3840 ./data/ln_O_L123_prunning_chr07.bam | ./hic_binning - 80 ./data/dm_res_grouping_details_O_chr07.txt ./data/O_cnv_winsize10000_wsize500kb_final.txt 10000 ./data/ln_Allele.ctg.chr07.table ./data/ln_O_against_dm.paf 1 CHR7 > torm_CHR7

samtools view -F 3840 ./data/ln_O_L123_prunning_chr08.bam | ./hic_binning - 25 ./data/dm_res_grouping_details_O_chr08.txt ./data/O_cnv_winsize10000_wsize500kb_final.txt 300000 ./data/ln_Allele.ctg.chr08.table ./data/ln_O_against_dm.paf 1 CHR8 > torm_CHR8

samtools view -F 3840 ./data/ln_O_L123_prunning_chr09.bam | ./hic_binning - 70 ./data/dm_res_grouping_details_O_chr09.txt ./data/O_cnv_winsize10000_wsize500kb_final.txt 300000 ./data/ln_Allele.ctg.chr09.table ./data/ln_O_against_dm.paf 1 CHR9 > torm_CHR9

# chr10 can be improved - not good
samtools view -F 3840 ./data/ln_O_L123_prunning_chr10.bam | ./hic_binning - 70 ./data/dm_res_grouping_details_O_chr10.txt ./data/O_cnv_winsize10000_wsize500kb_final.txt 300000 ./data/ln_Allele.ctg.chr10.table ./data/ln_O_against_dm.paf 1 CHR10 > torm_CHR10

samtools view -F 3840 ./data/ln_O_L123_prunning_chr11.bam | ./hic_binning - 70 ./data/dm_res_grouping_details_O_chr11.txt ./data/O_cnv_winsize10000_wsize500kb_final.txt 300000 ./data/ln_Allele.ctg.chr11.table ./data/ln_O_against_dm.paf 1 CHR11 > torm_CHR11

# chr12 can be improved - not good
samtools view -F 3840 ./data/ln_O_L123_prunning_chr12.bam | ./hic_binning - 70 ./data/dm_res_grouping_details_O_chr12.txt ./data/O_cnv_winsize10000_wsize500kb_final.txt 300000 ./data/ln_Allele.ctg.chr12.table ./data/ln_O_against_dm.paf 1 CHR12 > torm_CHR12

# size checking
# lg1 96140525
# lg2 100233849
# lg3 97703350
# lg4 97107539
# non-grouped 2481853

for i in 1 2 3 4; do awk -v hapi=$i '$5==hapi' s8_grouping_window_markers_refined_1st.txt | awk '{s+=$3-$2+1} END {print s}'; done
for i in 1 2 3 4; do awk -v hapi=$i '$5==hapi' s8_grouping_window_markers_refined_2nd.txt| awk '{s+=$3-$2+1} END {print s}'; done
for i in 1 2 3 4; do awk -v hapi=$i '$5==hapi' s8_grouping_window_markers.txt | awk '{s+=$3-$2+1} END {print s}'; done


# resolved.....: 393667116 = 96140525+100233849+97703350+97107539+2481853
# unique marker: 309579213


awk '$5==1' s8_grouping_window_markers_refined.txt | awk '{s+=$3-$2+1} END {print s}'
awk '$5==2' s8_grouping_window_markers_refined.txt | awk '{s+=$3-$2+1} END {print s}'
awk '$5==3' s8_grouping_window_markers_refined.txt | awk '{s+=$3-$2+1} END {print s}'
awk '$5==4' s8_grouping_window_markers_refined.txt | awk '{s+=$3-$2+1} END {print s}'
awk '$5==-1' s8_grouping_window_markers_refined.txt | awk '{s+=$3-$2+1} END {print s}'

cut -f1-4 s8_grouping_window_markers_refined.txt | awk '{s+=$3-$2+1} END {print s}'
cut -f1-4 s8_grouping_window_markers_refined.txt | uniq | awk '{s+=$3-$2+1} END {print s}'

# phasing acc (test): 0.997 - Potato_multipleCultivars/sc_add_check_hb_grouped_kmers/run_20221105/ 
