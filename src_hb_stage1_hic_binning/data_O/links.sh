# here we link the data for tests

#
cd /biodata/dep_mercier/grp_schneeberger/projects/methods/src_shq/z_GameteBinning_potato_multivar/src_hb_stage1_hic_binning/
ln -s /netscratch/dep_mercier/grp_schneeberger/projects/Potato_multipleCultivars/s2_10_Cultivars_PacBio_HiFi/a4_alignment_based_linkage_grouping/zhifiasm_v0p7_dm_O_asm/O_against_dm.paf ln_O_against_dm.paf


#
cd /biodata/dep_mercier/grp_schneeberger/projects/methods/src_shq/z_GameteBinning_potato_multivar/src_hb_stage1_hic_binning/data/
lg=chr12
ln -s /netscratch/dep_mercier/grp_schneeberger/projects/Potato_multipleCultivars/s2_10_Cultivars_PacBio_HiFi/a3_hic_alignment/bgi_deep_seq/hic_O/bwa_align/lg_wise_contigs_read_align//z_allelic_map_${lg}/Allele.ctg.table ln_Allele.ctg.${lg}.table

ln -s /netscratch/dep_mercier/grp_schneeberger/projects/Potato_multipleCultivars/s2_10_Cultivars_PacBio_HiFi/a3_hic_alignment/bgi_deep_seq/hic_O/bwa_align/lg_wise_contigs_read_align/z_bam_prunning_${lg}/prunning.bam ln_O_L123_prunning_${lg}.bam

ln -s /netscratch/dep_mercier/grp_schneeberger/projects/Potato_multipleCultivars/s2_10_Cultivars_PacBio_HiFi/a4_alignment_based_linkage_grouping/zhifiasm_v0p7_dm_O_asm/dm_res_grouping_details_O_${lg}.txt dm_res_grouping_details_O_${lg}.txt
