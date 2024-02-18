# visualize the grouping details: 1. how much bp grouped. 2. how each contig is grouped
#
# shq path if you would like to check data format: path="/path/to/a4_alignment_based_linkage_grouping/"
path="/your/working/path"
setwd(path)
#
pdf(paste(path, "/viz_s2_a4_ref_based_grouping_res_hifiasm_v0p7_c.pdf", sep=""), family="Helvetica", height= 8, width=8)
#par(mfrow=c(2,2), mai = c(0.7, 0.7, 0.2, 0.3)) # margin: bottom, left, top, right
m <- rbind(c(rep(1, 8), rep(2, 8)),
           c(rep(3, 8), rep(3, 8))
)
layout(m, heights = rep(1, 4) )
#
library("scales")
#
dm_12_ids    <- c("chr01", "chr02", "chr03", "chr04", "chr05", "chr06", "chr07", "chr08", "chr09", "chr10", "chr11", "chr12")
dm_12_lg_size<- read.table("/path/to/reference_sequence/DM_1-3_516_R44_potato_genome_assembly.v6.1_main12.chrids")
#
stie1_12_ids <- c("Stie_Chr_1_1_mapped_homLG_6_LG_20","Stie_Chr_2_1_mapped_homLG_11_LG_30","Stie_Chr_3_1_mapped_homLG_4_LG_17","Stie_Chr_4_1_mapped_homLG_9_LG_28","Stie_Chr_5_1_mapped_homLG_12_LG_1","Stie_Chr_6_1_mapped_homLG_3_LG_13","Stie_Chr_7_1_mapped_homLG_8_LG_37","Stie_Chr_8_1_mapped_homLG_10_LG_3","Stie_Chr_9_1_mapped_homLG_5_LG_29","Stie_Chr_10_1_mapped_homLG_1_LG_42","Stie_Chr_11_1_mapped_homLG_2_LG_12","Stie_Chr_12_1_mapped_homLG_7_LG_24")
stie1_12_lg_size <- read.table("/path/to/reference_genome/Stie1_genome.chr")
#
lg_12_ids <- c(dm_12_ids, stie1_12_ids)
#
dm_LG_group_size              <- data.frame(matrix(data = 0, nrow = 12, ncol = 11)) # row: LGs, column: samples
rownames(dm_LG_group_size)    <- c(paste("LG", 1:12, sep=""))
colnames(dm_LG_group_size)    <- c(LETTERS[1:10], "O")
#
stie1_LG_group_size           <- data.frame(matrix(data = 0, nrow = 12, ncol = 11)) # row: LGs, column: samples
rownames(stie1_LG_group_size) <- c(paste("LG", 1:12, sep=""))
colnames(stie1_LG_group_size) <- c(LETTERS[1:10], "O")
#
group_ctgcols  <- c(rep(rainbow(12), rep(2,12))) # grouped
lgctgcols      <- matrix("NA", nrow = 12, ncol = 4)
lgctgcols_all  <- matrix("NA", nrow = 12, ncol = 24)
for(lgx in 1:12)
{
  lgctgcols[lgx, 1]   <- "gray"     #  dm ref size col
  lgctgcols[lgx, 2]   <- group_ctgcols[2*lgx-1]
  lgctgcols[lgx, 3]   <- "gray50" #  stie1 ref size col
  lgctgcols[lgx, 4]   <- group_ctgcols[2*lgx]
  #
  lgctgcols_all[lgx, 1]     <- "gray"
  lgctgcols_all[lgx, 2:12]  <- group_ctgcols[2*lgx-1]
  lgctgcols_all[lgx, 13]    <- "gray50" #  stie1 ref size col
  lgctgcols_all[lgx, 14:24] <- group_ctgcols[2*lgx]
}
#
#
si <- 0
#for (sample in c("C"))
for(sample in c(LETTERS[1:10], "O"))  
{
  si <- si + 1
  # plot 1 & 2
  for (ref in c("dm", "stie1"))
  {
    grp_detail <- read.table(paste(path, "/zhifiasm_v0p7_", ref, "_", sample, "_asm/", ref, "_res_grouping_details_", sample, ".txt", sep=""))
    plot(grp_detail$V4, 
         grp_detail$V11, 
         xlim=c(0, max(grp_detail$V4)),
         ylim=c(0, max(grp_detail$V4)),
         col="red", 
         frame.plot = F, 
         xlab="ctg size", 
         ylab="ctg-cover-ref size", 
         main = paste(sample, " vs ", ref, sep=""))
    #
    grp_detail_main12 <- grp_detail[grp_detail$V7 %in% lg_12_ids, ]
    # legend
    legend("topleft",
           pch    = c(15, 15),
           col    = c("gray", "gray"),
           legend = c(paste("Grouped size: ", round(sum(grp_detail$V4)/1000000000, digits = 2), " Gb", sep=""),
                      paste("Grouped size: ", round(sum(grp_detail_main12$V4)/1000000000, digits = 2), " Gb in 12 LGs", sep="")  ),
           horiz  = FALSE,
           border = "NA",
           bty    = "n",
           cex    = 0.7)    
    # get size of each linkage group
    if(ref=="dm")
    {
      for(lgi in 1:12)
      {
        this_lg = dm_12_ids[lgi]
        dm_LG_group_size[lgi, si] <- sum(grp_detail_main12[grp_detail_main12$V7==this_lg, 4])
      }
    }else
    {
      for(lgi in 1:12)
      {
        this_lg = stie1_12_ids[lgi]
        stie1_LG_group_size[lgi, si] <- sum(grp_detail_main12[grp_detail_main12$V7==this_lg, 4])
      }      
    }
    #
  }
  # plot 3: lg-wise size comparison with two references
  #
  lgctgsizes <- round(cbind(dm_12_lg_size[, 2]*4, dm_LG_group_size[, si], stie1_12_lg_size[, 2]*4, stie1_LG_group_size[, si])/1000000, digits = 4)
  mp <- barplot(t(lgctgsizes), 
                beside=TRUE, 
                col=t(lgctgcols), 
                xlab="", 
                ylab="", 
                cex.lab=1, 
                axes = FALSE, 
                axisnames  = FALSE, 
                border=NA, 
                ylim=c(0,400), 
                xpd=FALSE)
  #
  axis(2, at=seq(0, 400,100), cex.axis=1)
  axis(1, at=seq(0.5, 60.5, 5), labels = FALSE, cex.axis=.8)
  # x-label on coverage
  text(mp[1,]+0.5, par("usr")[1]-4, labels = paste("LG-", 1:12, sep=""), srt = 60, adj = c(1.1,1.1), xpd = TRUE, cex=1)
  # axis labels
  mtext(text = "Size (Mb)", side = 2, line = 3.0, at = 200, cex=0.7, srt=90, col="black")  
  #
  legend(22,400,
         pch = c(15, 15, 15),
         col = c("gray", "gray50"),
         bty = "n",
         legend = c("DM size x4", "Stie1 size x4"),
         horiz = T,
         cex=0.9,
         border=NA)
  # "Otava (x4)"
  rect(xleft = seq(38,40.2,0.2), ybottom =  378.0, xright = seq(38.4,40.6, 0.2), ytop = 388.5,col = rainbow(12), border = "NA")
  text(x = 47, y = 382, labels = paste("Sample-", sample, " grouped size", sep=""), cex=0.9)
  #
}
# overall sizes
if(1)
{
  plot(NA,
       xlim=c(0, max(grp_detail$V4)),
       ylim=c(0, max(grp_detail$V4)),
       axes = FALSE, 
       xlab="",
       ylab="",
       frame.plot = F)
  plot(NA,
       xlim=c(0, max(grp_detail$V4)),
       ylim=c(0, max(grp_detail$V4)),
       axes = FALSE, 
       xlab="",
       ylab="",
       frame.plot = F)  
  #
  all_sizes <- cbind( dm_12_lg_size[, 2]*4, dm_LG_group_size, stie1_12_lg_size[, 2]*4, stie1_LG_group_size)/1000000
  colnames(all_sizes) <- c("dm", paste("dm_", c(LETTERS[1:10], "O"), sep=""), "stie1", paste("st1_", c(LETTERS[1:10], "O"), sep=""))
  mp <- barplot(t(all_sizes), 
                beside=TRUE, 
                col=t(lgctgcols_all), 
                xlab="", 
                ylab="", 
                cex.lab=1, 
                axes = FALSE, 
                axisnames  = FALSE, 
                border=NA, 
                ylim=c(0,400), 
                xpd=FALSE)
  #
  axis(2, at=seq(0, 400,100), cex.axis=1)
  axis(1, at=seq(0.5, 300.5, 25), labels = FALSE, cex.axis=.8)
  # x-label on coverage
  text(mp[1,]+10, par("usr")[1]-4, labels = paste("LG-", 1:12, sep=""), srt = 60, adj = c(1.1,1.1), xpd = TRUE, cex=1)
  # axis labels
  mtext(text = "Size (Mb)", side = 2, line = 3.0, at = 200, cex=0.7, srt=90, col="black")  
  #
  legend(22,400,
         pch = c(15, 15),
         col = c("gray", "gray50"),
         bty = "n",
         legend = c("DM size x4", "Stie1 size x4"),
         horiz = T,
         cex=0.9,
         border=NA)
  # "Otava (x4)"
  rect(xleft = seq(38,40.2,0.2)*4, ybottom =  378.0, xright = seq(38.4,40.6, 0.2)*4, ytop = 388.5,col = rainbow(12), border = "NA")
  text(x = 47*4, y = 382, labels = paste("All sample grouped size", sep=""), cex=0.9)
}
#
dev.off()
