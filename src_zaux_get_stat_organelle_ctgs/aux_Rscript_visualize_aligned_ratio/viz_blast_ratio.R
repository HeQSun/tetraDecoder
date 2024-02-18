# source code is here: /biodata/dep_mercier/grp_schneeberger/projects/methods/src_shq/z_GameteBinning_potato_version/src_zaux_get_stat_organelle_ctgs/aux_Rscript_visualize_aligned_ratio/

# working path
path <- "/netscratch/dep_mercier/grp_schneeberger/projects/Potato_single_cells/reference_manish_assembled_haplotye_aware/hifiasm_assembly_polish_and_select_v2/z_addi_s2_blast_limited_pollen_support_contigs/individual_ctgs_blast/ctgseqs_blast_results/"
#
pdf(paste(path, "/751ctgs_top1_hit_organelles_ratio.pdf", sep=""), height=8, width=9.26)
par(mfrow=c(2,1))
par(mai = c(1.0, 1.0, 0.5, 0.5)); # margin: bottom, left, top, right

# Note: 1.contig-id 2.size-bp 3.subject 4.coverage-bp 5.identity-bp; tab-separated. 
blast <- read.table("/netscratch/dep_mercier/grp_schneeberger/projects/Potato_single_cells/reference_manish_assembled_haplotye_aware/hifiasm_assembly_polish_and_select_v2/z_addi_s2_blast_limited_pollen_support_contigs/individual_ctgs_blast/ctgseqs_blast_results/all_751_ctg_blast_top1_stat.txt", sep="\t")

thecolors <-c()
ratio     <-c()
total_ctg_size      <- 0
total_coverage_size <- 0
for (bi in 1:length(blast$V1))
{
  # get ratio
  ratio[bi] = blast[bi, 4]/ blast[bi, 2]
  total_ctg_size = total_ctg_size + blast[bi, 2]
  if(ratio[bi] > 1)
  {
    ratio[bi] = 1
    total_coverage_size = total_coverage_size + blast[bi, 2]
  }else
  {
    total_coverage_size = total_coverage_size + blast[bi, 4]
  }
  # prepare color
  if(grepl("mitochon", as.character(blast[bi, 3]), fixed=T) == TRUE)
  {
    thecolors[bi] <- "gold"
  }else
  if(grepl("chloro", as.character(blast[bi, 3]), fixed=T) == TRUE)
  {
    thecolors[bi] <- "purple"
  }else    
  if(grepl("plastid", as.character(blast[bi, 3]), fixed=T) == TRUE)
  {
    thecolors[bi] <- "purple"
  }else    
  {
    thecolors[bi] <- "gray"
  }
}

#
plot(log10(blast$V2), 
     ratio, 
     col=thecolors, 
     frame.plot = F, 
     ylim=c(0, 1),
     xlab="Log10(ctg size: bp)", 
     ylab="Ratio",
     pch=4,
     main="Coverage ratio of 751 contigs by organelles - top hit only", cex.main=0.95)
##
##
legend("topright",
       pch    = c(15,15,15),
       col    = c("gold", "purple","gray"),
       legend = c("Mitochondrion", "Plastid (chloroplast)", "Other"),
       horiz  = FALSE,
       border = "NA",
       bty    = "n",
       cex    = 0.8)
#
legend(4.67, 0.57,
       pch    = c(15,15,15),
       col    = c("black", "black", "black"),
       legend = c(paste("Total limited-pollen-covered contig size: ", total_ctg_size, " bp", sep=""),
                  paste("Total top-1-hit of organelle covered length: ", total_coverage_size, " bp", sep=""),
                  paste("Total top-1-hit of organelle covered ratio: ", round(total_coverage_size/total_ctg_size, digits = 4), sep="") ),
       horiz  = FALSE,
       border = "red",
       #bty    = "n",
       cex    = 0.8)
#
hist(ratio, breaks=50, border = "NA", col="red", main="Coverage ratio of 751 contigs by organelles", cex.main=0.95, xlim=c(0, 1), xlab="Ratio")

#
dev.off()







