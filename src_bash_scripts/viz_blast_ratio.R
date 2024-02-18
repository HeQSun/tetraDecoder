#!/bin/env Rscript

## get options from cmd
###############################################################################
args<-commandArgs(TRUE)
if(length(args) < 3) {
  cat(length(args))
  cat("\n   Info: visualize blast result of contigs (to NCBI nucleotide database).\n")
  cat("   USAGE: viz_blast_ratio.R ctg_blast_top1_stat.txt working_path out_prefix\n\n")
}else
{
  ## check args
  top1_stat_file = args[[1]]
  working_path   = args[[2]]
  out_prefix     = args[[3]]
  #
  selected_as_organ_ratio = 0.95
  # working path
  path <- working_path
  #
  pdf(paste(path, "/", out_prefix, "_ctgs_top1_hit_organelles_ratio.pdf", sep=""), height=8, width=9.26)
  par(mfrow=c(2,1))
  par(mai = c(1.0, 1.0, 0.5, 0.5)); # margin: bottom, left, top, right  
  # Note: 1.contig-id 2.size-bp 3.subject 4.coverage-bp 5.identity-bp; tab-separated. 
  if(file.info(top1_stat_file)$size > 0)
  {
    blast <- read.table(top1_stat_file, sep="\t")
    cov_sub <- blast[blast$V4/blast$V2>selected_as_organ_ratio, ]
    write.table(x=unique(cov_sub[, 1:2]), file=paste(path, "/", out_prefix, "_ctgs_top1_hit_subset_organelle_unique.ctgsizes", sep=""), quote=F,row.names=F, col.names=F, append=F)  
    write.table(x=unique(cov_sub$V1), file=paste(path, "/", out_prefix, "_ctgs_top1_hit_subset_organelle_unique.ctgids", sep=""), quote=F,row.names=F, col.names=F, append=F)
    #
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
         xlim=c(4, log10(350000)),
         xlab="Log10(ctg size: bp)", 
         ylab="Ratio",
         pch=4,
         main=paste("Coverage ratio of ", length(blast$V1)," contigs by organelles - top hit only", sep=""), cex.main=0.95)
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
           legend = c(paste("Total ......................... contig size: ", total_ctg_size, " bp", sep=""),
                      paste("Total top-1-hit of organelle covered length: ", total_coverage_size, " bp", sep=""),
                      paste("Total top-1-hit of organelle covered ratio.: ", round(total_coverage_size/total_ctg_size, digits = 4), sep="") ),
           horiz  = FALSE,
           border = "red",
           #bty    = "n",
           cex    = 0.8)
    #
    hist(ratio, breaks=min(10, length(blast$V1)), border = "NA", col="red", main=paste("Coverage ratio of ", length(blast$V1)," contigs by organelles - top hit only", sep=""), cex.main=0.95, xlim=c(0, 1), xlab="Ratio")  
    abline(v=selected_as_organ_ratio, col="gray")  
    text(x=0.48, y=length(blast$V1)/3, label=paste("Note, ", length(cov_sub[,1]), " contigs overlapping >", selected_as_organ_ratio, " organelle sequences: ", sum(cov_sub$V2), " bp; to be excluded!", sep=""), col="red")
  }else
  {
    #
    plot(0, 
         0, 
         frame.plot = F, 
         ylim=c(0, 1),
         xlim=c(4, log10(350000)),
         xlab="Log10(ctg size: bp)", 
         ylab="Ratio",
         pch=4,
         main=paste("No contigs hit by organelles", sep=""), cex.main=0.95)    
    #
    write.table(x=c("no-ctg-hit-by-organelle 0"), file=paste(path, "/", out_prefix, "_ctgs_top1_hit_subset_organelle_unique.ctgsizes", sep=""), quote=F,row.names=F, col.names=F, append=F)  
    write.table(x="no-ctg-hit-by-organelle", file=paste(path, "/", out_prefix, "_ctgs_top1_hit_subset_organelle_unique.ctgids", sep=""), quote=F,row.names=F, col.names=F, append=F)
  }
  #
  dev.off()
}






