#
library(scales)
##
my_plot<-function(df, chr, chrsize, speciescol, this_chr_cov)
{
  plot(0,0, col="white", cex=0.1, xlim=c(0, chrsize[chr]), ylim=c(-1.25, +1.25), axes = FALSE, xpd  = FALSE,  xlab="", ylab="")
  # 0.9-1.2 plot coverage
  ylg  = 0.8
  yref = -1.2
  #
  # LG
  segments(x0 = 1, y0 = ylg,    x1 = max(df$starta, df$enda),     y1 = ylg,  col = alpha("cornflowerblue", 0.8), lwd=3, lend=0)
  # Ref
  segments(x0 = 1, y0 = yref, x1 = max(df$startb, df$endb),       y1 = yref, col = alpha("orangered", 0.8),      lwd=3, lend=0)
  axis(side = 1, at=seq(0, chrsize[chr], 10e+06), labels=rep("", length(seq(0, chrsize[chr], 10e+06)) ), lty=1, col="gray", tck=-0.05)
  axis(side = 1, at=seq(0, chrsize[chr], 10e+05), labels=rep("", length(seq(0, chrsize[chr], 10e+05)) ), lty=1, col="gray", tck=-0.02)
  text(c(seq(0, chrsize[chr], 10e+06), chrsize[chr]+5e+06)*1.02, -1.8, 
       labels = c(seq(0, chrsize[chr], 10e+06)/1e+06, "(Mb)"),
       srt = 0, adj = c(1.1,1.1), xpd = TRUE, cex=0.8, srt = 00, col="gray")    
  #
  axis(side = 3, at=seq(0, chrsize[chr], 10e+06), labels=rep("", length(seq(0, chrsize[chr], 10e+06)) ), lty=1, col="gray", tck=-0.05)
  axis(side = 3, at=seq(0, chrsize[chr], 10e+05), labels=rep("", length(seq(0, chrsize[chr], 10e+05)) ), lty=1, col="gray", tck=-0.02)
  # tune position of variations plots
  positioning1 <- c(-1, -1, +1, +1)*0.10
  ######################################### plot lg-wise assembly and ref ####################################################
  plot_lg_ref = 1
  if(plot_lg_ref)
  {
    # syn over 5 kb
    syn <- df[df$enda-df$starta > 0 & df$type=="SYN", ]
    num = length(syn$chra)
    if(num>0)
      for(i in c(1:num))
      {
        x<-c(syn[i, 2], syn[i, 3], syn[i, 6], syn[i, 5])
        y<-c(ylg,  ylg,  yref,  yref)  + positioning1
        polygon(x, y,
                col = alpha("gray", 1),
                border = NA)
      }
  }
  ########## add coverage info
  meancov <- mean(this_chr_cov$V6)+0.000001 # avoid 0
  meancov <- mean(this_chr_cov[this_chr_cov$V6>=0.2*meancov & this_chr_cov$V6<=1.8*meancov, 6])+0.000001
  points(x = (this_chr_cov[, 2]+this_chr_cov[, 3])/2, 
         y = this_chr_cov[, 6]/meancov*0.2 + rep(ylg, length(this_chr_cov[, 6]) ), 
         col="purple", type="p", cex=0.02)
  ##########
  return (0)
}
alignmraw_process<-function(alignmraw)
{
  len <- length(alignmraw$V1)  
  for(r in (1:len))
  {
    if(alignmraw[r, 4] == "-")
    {
      tmp             = alignmraw[r, 7]
      alignmraw[r, 7] = alignmraw[r, 6]
      alignmraw[r, 6] = tmp
    }
  }
  return(alignmraw)
}
#
## get options from cmd
###############################################################################
args<-commandArgs(TRUE)
if(length(args) < 5) {
  cat(length(args))
  cat("\n   Info: visualize alignment of new chr assembly to DM reference genome.\n")
  cat("   USAGE: Rscript visualize_alignment_four_hapChr_and_refDM_interhap_phase_correction.R 1.cultivar 2.strand 3.align_file 4.out_path 5.out_prefix\n\n")
}else
{
  ## check args: e.g., cultivar="A";strand="_phase_corr";aln_file="fw_data";path="./";out_prefix="A_phase_corr_refDM_aln"
  cultivar   = args[[1]]
  strand     = args[[2]]
  aln_file   = args[[3]]
  path       = args[[4]]
  out_prefix = args[[5]]
  #
  # cultivar="A"
  # strand="_phase_corr"
  # aln_file="fw_data"
  # path="../"
  # out_prefix="A_phase_corr_refDM_aln2_with_cov"
  #
  ################ output pdf file ###########################################
  pdf(paste(path,"/",out_prefix,".pdf", sep=""), family="Helvetica", height=5, width=7.08*2)
  #par(mfrow=c(4,2), mai = c(0.4, 0.5, 0.4, 0.5)); # margin: bottom, left, top, right
  par(mfrow=c(2,2), mai = c(0.5, 0.1, 0.5, 0.3)); # margin: bottom, left, top, right
  #  
  #
  chri <- 0
  for (lg in c(paste("chr0", 1:9, sep=""), "chr10", "chr11", "chr12"))
  ####for (lg in c(paste("chr0", 1, sep="")))
  {
    chri=chri+1    
    real_hap1=(chri-1)*4+1
    real_hap2=(chri-1)*4+2
    real_hap3=(chri-1)*4+3
    real_hap4=(chri-1)*4+4
    hapstr   =paste(real_hap1, real_hap2, real_hap3, real_hap4, sep=",")
    haps <- strsplit(hapstr, ',')
    haps <- as.vector(haps[[1]])
    #
    #### cov <- read.table(paste("/Users/sun/Desktop/z_10potato_project/xjtu_a9_find_regions_not_assembled/src_find_gaps/data/s_read_realign_chr01/cnv_winsize50000_step50000_hq.txt",sep=""))
    cov <- read.table(paste("/netscratch/dep_mercier/grp_schneeberger/projects/Potato_multipleCultivars/s2_10_Cultivars_PacBio_HiFi/a8_make_full_genome_tmp/tmp_full_", cultivar,"/shallow_DNB_aln/cnv_winsize50000_step50000_hq.txt", sep=""))
  # example input
  # cultivar="A"
  # lg="chr01"
  # hapstr="1,2,3,4"
  # strand="_fw" # or "_rev"
  # aln_file="fw_data"
  # path<-"/netscratch/dep_mercier/grp_schneeberger/projects/Potato_multipleCultivars/s2_10_Cultivars_PacBio_HiFi/a6_lg_wise_scaffolding/scaf_A/s_read_realign_chr01"
  # out_prefix="A_chr01"

    ####################################################################################################################
    ################################# double monoploid ref: totalN_12047830_of_860168943 #################################
    for(hapi in haps)
    {
      # cov
      this_chr_cov <- cov[cov$V1==paste(lg, "_hap", hapi, sep=""), ]
      # aln
      if(strand=="_zfinal")
      {
        alignmraw2       <- read.table(paste(path, "/s_read_realign_",lg,"/juicebox_correction_",hapi,"/zcorrection_",hapi,"_refDM",strand,"/", aln_file, sep=""))
      }    
      if(strand=="_zfinal.oneseq")
      {
        alignmraw2       <- read.table(paste(path, "/s_read_realign_",lg,"/juicebox_correction_",hapi,"/zcorrection_",hapi,"_refDM",strand,"/", aln_file, sep=""))
      }
      if(strand=="_phase_corr")
      {
        alignmraw2       <- read.table(paste(path, "/s_read_realign_",lg,"/zcorrection_four_hap_juicebox_DM_strand_interhap_correction/zcorrection_",hapi,strand,"_4th_refDM/", aln_file, sep=""))
      }                
      alignmraw        <- alignmraw_process(alignmraw2)
      len              <- length(alignmraw$V1)
      alignm2t         <- cbind(alignmraw[, c(1,2,3,5,6,7)], paste(rep("SYN", len), 1:len, sep=""), rep("SYN", len))
      alignm           <- alignm2t[alignm2t$V3 - alignm2t$V2 > 1000, ]
      #alignm          <- alignm2t
      colnames(alignm) <- c("chra", "starta", "enda", "chrb", "startb", "endb", "id", "type")
      #
      alignm6         <- as.data.frame(alignm[alignm$chra==paste(lg, "_", "hap", hapi, sep=""), ])
      df              <- alignm6
      colnames(df)    <- c("chra", "starta", "enda", "chrb", "startb", "endb", "id", "type")
      #df             <- alignm6
      chr             <- 1
      chrsize         <- max(alignm6$enda, alignm6$endb, 100000000)
      speciescol      <- "mediumorchid3"
      my_plot(df, chr, chrsize, speciescol, this_chr_cov)
      title(main = paste("DM-Ref v.s. ", cultivar, ": ", lg, "-", hapi, sep=""), cex.main=0.8, font.main=1)
      # legend
      legend(90000000, -0.2,
             pch    = c(15,15),
             col    = c("cornflowerblue", "orangered"),
             legend = c(paste("Cultivar-", cultivar, sep=""),
                        paste("DM-Ref", sep="")
             ),
             horiz  = FALSE,
             border = "NA",
             bty    = "n",
             cex    = 0.6)  
      #
    }
  }
  #
  ##########################################################################################################################
  #
  dev.off()
}









