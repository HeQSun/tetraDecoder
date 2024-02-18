tig_marker_finder cnv_winsize10000_step10000_hq_merged_vs_hifi.txt 113 30 50000 test

echo "   Info: checking intermediate output file: cnv_winsize10000_step10000_hq_merged_vs_hifi_markers_test.txt"
grep -v "#" cnv_winsize10000_step10000_hq_merged_vs_hifi_markers_test.txt | awk '{s+=$3-$2+1} END {print "   Info: total marker covering size: "s" bp\n"}' 

echo "   Info: checking final output file: cnv_winsize10000_step10000_hq_merged_vs_hifi_markers_test_wsize50kb_final.txt"
grep -v "#" cnv_winsize10000_step10000_hq_merged_vs_hifi_markers_test_wsize50kb_final.txt | awk '{s+=$3-$2+1} END {print "   Info: total marker covering size: "s" bp\n"}' 

# check consistency in ctg sizes and window sizes
# cut -f1,6 cnv_winsize10000_step10000_hq_markers_test.txt | grep 'utg' | uniq | cut -f1 > ctgid_torm
# >tmp;while read r; do grep $r cnv_winsize10000_step10000_hq_merged_vs_hifi_markers_test.txt | awk -v a="$r" '{s+=$3-$2+1} END {print a"\t"s}' >> tmp; done < ctgid_torm
# awk '{s+=$2} END {print s} ' tmp 

#=> /netscratch/dep_mercier/grp_schneeberger/projects/Potato_single_cells/reference_manish_assembled_haplotye_aware/hifiasm_assembly_polish/polished_asmsize_2437914205bp
