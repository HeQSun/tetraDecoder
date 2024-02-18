#
library(scales)
##
my_plot<-function(df, chr, chrsize, speciescol)
{
  plot(0,0, col="white", cex=0.1, xlim=c(0, chrsize[chr]), ylim=c(-1.25, +1.25), axes = FALSE, xpd  = FALSE,  xlab="", ylab="")
  #
  ylg  = 1.2
  yref = -1.2
  #
  # LG
  segments(x0 = 1, y0 = ylg,    x1 = max(df$starta, df$enda),     y1 = ylg,  col = alpha("cornflowerblue", 0.8), lwd=3, lend=0)
  # Ref
  segments(x0 = 1, y0 = yref, x1 = max(df$startb, df$endb),       y1 = yref, col = alpha("orangered", 0.8),      lwd=3, lend=0)
  axis(side = 1, at=seq(0, chrsize[chr], 10e+06), labels=rep("", length(seq(0, chrsize[chr], 10e+06)) ), lty=1, col="gray", tck=-0.05)
  axis(side = 1, at=seq(0, chrsize[chr], 10e+05), labels=rep("", length(seq(0, chrsize[chr], 10e+05)) ), lty=1, col="gray", tck=-0.02)
  text(c(seq(0, chrsize[chr], 10e+06), 11e+07)*1.04, -1.8, 
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
path<-"/netscratch/dep_mercier/grp_schneeberger/projects/Potato_multipleCultivars/s2_10_Cultivars_PacBio_HiFi/a6_lg_wise_scaffolding/scaf_A/s_read_realign_chr01/juicebox_correction_1/zcorrection_1_refDM_rev/"

aln_file="/netscratch/dep_mercier/grp_schneeberger/projects/Potato_multipleCultivars/s2_10_Cultivars_PacBio_HiFi/a6_lg_wise_scaffolding/scaf_A/s_read_realign_chr01/juicebox_correction_1/zcorrection_1_refDM/fw_data"
reflg=1
cultivar="A"

# one page including 4 plots of the same homologous linkage groups
pdf(paste(path, "/refDM_outprefix.pdf", sep=""), family="Helvetica", height=6, width=7.08)
par(mfrow=c(4,2), mai = c(0.4, 0.5, 0.4, 0.5)); # margin: bottom, left, top, right
#
#
if(reflg==1)
{
  newlgs <- c("chr01_hap1", "chr01_hap2","chr01_hap3","chr01_hap4")
}else
  if(reflg==2)
  {
    newlgs <- c("chr02_hap5", "chr02_hap6","chr02_hap7","chr02_hap8")
  }else
    if(reflg==3)
    {
      newlgs <- c("chr03_hap9", "chr03_hap10","chr03_hap11","chr03_hap12")
    }else    
      if(reflg==4)
      {
        newlgs <- c("chr04_hap13", "chr04_hap14","chr04_hap15","chr04_hap16")
      }else    
        if(reflg==5)
        {
          newlgs <- c("chr05_hap17", "chr05_hap18","chr05_hap19","chr05_hap20")
        }else    
          if(reflg==6)
          {
            newlgs <- c("chr06_hap21", "chr06_hap22","chr06_hap23","chr06_hap24")
          }else    
            if(reflg==7)
            {
              newlgs <- c("chr07_hap25", "chr07_hap26","chr07_hap27","chr07_hap28")
            }else    
              if(reflg==8)
              {
                newlgs <- c("chr08_hap29", "chr08_hap30","chr08_hap31","chr08_hap32")
              }else    
                if(reflg==9)
                {
                  newlgs <- c("chr09_hap33", "chr09_hap34","chr09_hap35","chr09_hap36")
                }else    
                  if(reflg==10)
                  {
                    newlgs <- c("chr10_hap37", "chr10_hap38","chr10_hap39","chr10_hap40")
                  }else    
                    if(reflg==11)
                    {
                      newlgs <- c("chr11_hap41", "chr11_hap42","chr11_hap43","chr11_hap44")
                    }else    
                      if(reflg==12)
                      {
                        newlgs <- c("chr12_hap45", "chr12_hap46","chr12_hap47","chr12_hap48")
                      }
#
lgorder <- 0
#
for(lg in newlgs )
{
  lgorder <- lgorder + 1
  ####################################################################################################################
  ################################# double monoploid ref: totalN_12047830_of_860168943 #################################
  if(1)
  {
    alignmraw2       <- read.table(aln_file)
    alignmraw        <- alignmraw_process(alignmraw2)
    len              <- length(alignmraw$V1)
    alignm2t         <- cbind(alignmraw[, c(1,2,3,5,6,7)], paste(rep("SYN", len), 1:len, sep=""), rep("SYN", len))
    alignm           <- alignm2t[alignm2t$V3 - alignm2t$V2 > 10000, ]
    #alignm          <- alignm2t
    colnames(alignm) <- c("chra", "starta", "enda", "chrb", "startb", "endb", "id", "type")
    #
    alignm6         <- as.data.frame(alignm[alignm$chra==lg, ])
    df              <- alignm6
    colnames(df)    <- c("chra", "starta", "enda", "chrb", "startb", "endb", "id", "type")
    #df             <- alignm6
    chr             <- 1
    chrsize         <- max(alignm6$enda, alignm6$endb, 100000000)
    speciescol      <- "mediumorchid3"
    my_plot(df, chr, chrsize, speciescol)
    title(main = paste("DM LG ", reflg, " v.s. ", cultivar, lg, sep=""), cex.main=0.8, font.main=1)
    # legend
    legend(90000000, -0.2,
             pch    = c(15,15),
             col    = c("cornflowerblue", "orangered"),
             legend = c(paste("cultivar-", cultivar, sep=""),
                        paste("DM", sep="")
             ),
             horiz  = FALSE,
             border = "NA",
             bty    = "n",
             cex    = 0.6)  
  }
  ##########################################################################################################################
  #
}
dev.off()











