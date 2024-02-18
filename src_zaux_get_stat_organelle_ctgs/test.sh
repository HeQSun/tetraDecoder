ctgsize=/netscratch/dep_mercier/grp_schneeberger/projects/Potato_single_cells/reference_manish_assembled_haplotye_aware/hifiasm_assembly_polish_and_select_v2/s0_select_long_contigs_v2/exclude_subset_5/HiFiasm_ref_pilon_6366long_ctgs_selected.ctgsizes

get_organelle_ctg_stat utg010854l_pilon_againt_ncbi_nt_nt_all.oblast ${ctgsize} > utg010854l_pilon_againt_ncbi_nt_nt_all.report

grep -A 5 'Report: ' utg010854l_pilon_againt_ncbi_nt_nt_all.report > ctg_blast_top3_stat.txt
grep -A 1 'Report' ctg_blast_top3_stat.txt | grep -v 'Report' | grep -v '\-\-' | sed 's/           //g' > ctg_blast_top1_stat.txt.txt
grep -v 'plastid' ctg_blast_top1_stat.txt.txt | grep -v 'chloropla' | grep -v 'mitochon' > non_organelle.txt
