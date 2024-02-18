links <- read.table("/biodata/dep_mercier/grp_schneeberger/projects/methods/src_shq/z_GameteBinning_potato_multivar/src_hb_stage1_hic_binning/s5_test_raw_tig_marker_cross_link_count/all.txt")

hist(links$V11, col="deepskyblue", breaks = 1000, border=NA, xlim=c(0, 200), ylim=c(0, 100), main="", xlab="read link per 100 kb")

links_hap <- read.table("/biodata/dep_mercier/grp_schneeberger/projects/methods/src_shq/z_GameteBinning_potato_multivar/src_hb_stage1_hic_binning/s5_test_raw_tig_marker_cross_link_count/hap.txt")
links_hap_select <- links_hap[links_hap$V11>5, ]

hist(links_hap$V11, col="red", breaks = 1000, border=NA, xlim=c(0, 50), ylim=c(0, 100),main="", xlab="read link per 100 kb")
hist(links_hap_select$V11, col="orangered", breaks = 1000, border=NA, xlim=c(0, 50), ylim=c(0, 100),main="", xlab="read link per 100 kb")
unique(c(as.character(links_hap_select$V1), as.character(links_hap_select$V3)))

link_spe <- links_hap[links_hap$V1=="utg000035lc" & links_hap$V11>5, ]