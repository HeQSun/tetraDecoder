cd /path/to/src_hb_stage1_hic_binning/

cd s6_pruned_bam_raw_tig_marker_cross_link_count/

fdp -Tpdf s6_haplotig_hic_contact_matrix_subclusters_raw.dot > s6_haplotig_hic_contact_matrix_subclusters_raw.pdf
fdp -Tpdf s6_haplotig_hic_contact_matrix_subclusters_raw_merged.dot > s6_haplotig_hic_contact_matrix_subclusters_raw_merged.pdf
fdp -Tpdf s7_haplotig_hic_contact_matrix_subclusters_hap_extended.dot > s7_haplotig_hic_contact_matrix_subclusters_hap_extended.pdf

dot_comparator ../gamete_binning_result_s3_refLG1_haplotig_GT_similarity_matrix_subclusters_final.dot s6_haplotig_hic_contact_matrix_subclusters_raw_merged.dot
dot_comparator ../gamete_binning_result_s3_refLG1_haplotig_GT_similarity_matrix_subclusters_final.dot s7_haplotig_hic_contact_matrix_subclusters_hap_extended.dot

dot_evaluator ../ln_O_s4p6_refine_grouping_final_window_markers_sorted.txt s6_haplotig_hic_contact_matrix_subclusters_raw.dot
dot_evaluator ../ln_O_s4p6_refine_grouping_final_window_markers_sorted.txt s7_haplotig_hic_contact_matrix_subclusters_hap_extended.dot


# working case with all default in coding 20221101: chr08 and chr10 are extreme cases with one haplotype with nearly no haplotigs.

sample="O"

samtools view -F 3840 ./data_${sample}/ln_${sample}_L123_prunning_chr01.bam | ./hic_binning --min-hic 70 --hap-win-r 1 --min-hap-size 300000 --max-allelic-r 0.02 --ctg ./data_${sample}/dm_res_grouping_details_${sample}_chr01.txt --win-marker ./data_${sample}/${sample}_cnv_winsize10000_wsize500kb_final.txt --gmap-allelic ./data_${sample}/ln_Allele.ctg.chr01.table --align-paf ./data_${sample}/ln_${sample}_against_dm.paf -o CHR1_${sample} > torm_CHR1_${sample}

sample="A"

samtools view -F 3840 ./data_${sample}/ln_${sample}_L123_prunning_chr01.bam | ./hic_binning --min-hic 45 --hap-win-r 1 --min-hap-size 100000 --max-allelic-r 0.02 --ctg ./data_${sample}/dm_res_grouping_details_${sample}_chr01.txt --win-marker ./data_${sample}/${sample}_cnv_winsize10000_wsize500kb_final.txt --gmap-allelic ./data_${sample}/ln_Allele.ctg.chr01.table --align-paf ./data_${sample}/ln_${sample}_against_dm.paf -o CHR1_${sample} > torm_CHR1_${sample}

sample="B"

samtools view -F 3840 ./data_${sample}/ln_${sample}_L123_prunning_chr02.bam | ./hic_binning --min-hic 10.0 --min-hic-i 0.55 --hap-win-r 1 --min-hap-size 10000 --max-allelic-r 0.05 --ctg ./data_${sample}/dm_res_grouping_details_${sample}_chr02.txt --win-marker ./data_${sample}/${sample}_cnv_winsize10000_wsize500kb_final.txt --gmap-allelic ./data_${sample}/ln_Allele.ctg.chr02.table --align-paf ./data_${sample}/ln_${sample}_against_dm.paf -o CHR2_${sample} > torm_CHR2_${sample}

sample="C"
lg="chr03"

samtools view -F 3840 ./data_${sample}/ln_${sample}_L123_prunning_${lg}.bam | ./hic_binning --min-hic 20.0 --min-hic-i 1.1 --hap-win-r 1 --min-hap-size 10000 --max-allelic-r 0.1 --ctg ./data_${sample}/dm_res_grouping_details_${sample}_${lg}.txt --win-marker ./data_${sample}/${sample}_cnv_winsize10000_wsize500kb_final.txt --gmap-allelic ./data_${sample}/ln_Allele.ctg.${lg}.table --align-paf ./data_${sample}/ln_${sample}_against_dm.paf -o ${lg}_${sample} > torm_${lg}_${sample}

samtools view -F 3840 ./data_${sample}/ln_${sample}_L123_prunning_${lg}.bam | ./hic_binning --min-hic 20.0 --min-hic-i 1.1 --hap-win-r 1 --min-hap-size 200000 --max-allelic-r 0.05 --ctg ./data_${sample}/dm_res_grouping_details_${sample}_${lg}.txt --win-marker ./data_${sample}/${sample}_cnv_winsize10000_wsize500kb_final.txt --gmap-allelic ./data_${sample}/ln_Allele.ctg.${lg}.table --align-paf ./data_${sample}/ln_${sample}_against_dm.paf -o ${lg}_${sample} > torm_${lg}_${sample}

sample="C"
lg="chr05"

samtools view -F 3840 ./data_${sample}/ln_${sample}_L123_prunning_${lg}.bam | ./hic_binning --min-hic 20.0 --min-hic-i 1.1 --hap-win-r 1 --min-hap-size 200000 --max-allelic-r 0.05 --ctg ./data_${sample}/dm_res_grouping_details_${sample}_${lg}.txt --win-marker ./data_${sample}/${sample}_cnv_winsize10000_wsize500kb_final.txt --gmap-allelic ./data_${sample}/ln_Allele.ctg.${lg}.table --align-paf ./data_${sample}/ln_${sample}_against_dm.paf -o ${lg}_${sample} > torm_${lg}_${sample}

sample="C"
lg="chr06"

samtools view -F 3840 ./data_${sample}/ln_${sample}_L123_prunning_${lg}.bam | ./hic_binning --min-hic 30.0 --min-hic-i 1.1 --hap-win-r 1 --min-hap-size 100000 --max-allelic-r 0.1 --ctg ./data_${sample}/dm_res_grouping_details_${sample}_${lg}.txt --win-marker ./data_${sample}/${sample}_cnv_winsize10000_wsize500kb_final.txt --gmap-allelic ./data_${sample}/ln_Allele.ctg.${lg}.table --align-paf ./data_${sample}/ln_${sample}_against_dm.paf -o ${lg}_${sample} > torm_${lg}_${sample}

sample="C"
lg="chr10"

samtools view -F 3840 ./data_${sample}/ln_${sample}_L123_prunning_${lg}.bam | ./hic_binning --min-hic 10.0 --min-hic-i 0.55 --hap-win-r 1 --min-hap-size 100000 --max-allelic-r 0.05 --ctg ./data_${sample}/dm_res_grouping_details_${sample}_${lg}.txt --win-marker ./data_${sample}/${sample}_cnv_winsize10000_wsize500kb_final.txt --gmap-allelic ./data_${sample}/ln_Allele.ctg.${lg}.table --align-paf ./data_${sample}/ln_${sample}_against_dm.paf -o ${lg}_${sample} > torm_${lg}_${sample}


# size checking
# lg1 96140525
# lg2 100233849
# lg3 97703350
# lg4 97107539
# non-grouped 2481853

for i in 1 2 3 4; do awk -v hapi=$i '$5==hapi' s8_grouping_window_markers_refined_1st.txt | awk '{s+=$3-$2+1} END {print s}'; done
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
