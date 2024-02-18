# here we link the data for tests

#
sample="A"
#
cd /biodata/dep_mercier/grp_schneeberger/projects/methods/src_shq/z_GameteBinning_potato_multivar/src_hb_stage1_hic_binning/data_${sample}/

cp -a /netscratch/dep_mercier/grp_schneeberger/projects/Potato_multipleCultivars/s2_10_Cultivars_PacBio_HiFi/a2_initial_assembly/sample_${sample}/hifiasm_asm_v0p7/subset4_illu_re_align_purge_ovl/${sample}_cnv_winsize10000_step10000_hq_merged_vs_hifi_markers_20221102_wsize500kb_final.txt ${sample}_cnv_winsize10000_wsize500kb_final.txt


lg=chr01

ln -s /netscratch/dep_mercier/grp_schneeberger/projects/Potato_multipleCultivars/s2_10_Cultivars_PacBio_HiFi/a4_alignment_based_linkage_grouping/zhifiasm_v0p7_dm_${sample}_asm/${sample}_against_dm.paf ln_${sample}_against_dm.paf

ln -s /netscratch/dep_mercier/grp_schneeberger/projects/Potato_multipleCultivars/s2_10_Cultivars_PacBio_HiFi/a3_hic_alignment/bgi_deep_seq/hic_${sample}/bwa_align/lg_wise_contigs_read_align//z_allelic_map_${lg}/Allele.ctg.table ln_Allele.ctg.${lg}.table

ln -s /netscratch/dep_mercier/grp_schneeberger/projects/Potato_multipleCultivars/s2_10_Cultivars_PacBio_HiFi/a3_hic_alignment/bgi_deep_seq/hic_${sample}/bwa_align/lg_wise_contigs_read_align/z_bam_prunning_${lg}/prunning.bam ln_${sample}_L123_prunning_${lg}.bam

ln -s /netscratch/dep_mercier/grp_schneeberger/projects/Potato_multipleCultivars/s2_10_Cultivars_PacBio_HiFi/a4_alignment_based_linkage_grouping/zhifiasm_v0p7_dm_${sample}_asm/dm_res_grouping_details_${sample}_${lg}.txt dm_res_grouping_details_${sample}_${lg}.txt
