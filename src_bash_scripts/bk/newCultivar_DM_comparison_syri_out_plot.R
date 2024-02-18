##
library(scales)
##
my_plot<-function(df, chr, speciescol)
{
  plot(0,0, col="white", cex=0.1, xlim=c(0, 100e+06), ylim=c(-1.2, +1.2), axes = FALSE, xpd  = FALSE,  xlab="", ylab="")
  ##
  # Rojo Pasion apricot: currot is at bottom
  segments(x0 = 1, y0 = -1.14, x1 = max(df$starta, df$enda),   y1 = -1.14, col = alpha("cornflowerblue", 0.8), lwd=3, lend=1)
  # other related species
  segments(x0 = 1, y0 = +1.14, x1 = max(df$startb, df$endb),   y1 = +1.14, col = alpha(speciescol, 0.8), lwd=3, lend=1)
  #
  # axis
  axis(side = 1, at=seq(0, chrsize[chr], 10e+06), labels=rep("", length(seq(0, chrsize[chr], 10e+06)) ), lty=1, col="gray", tck=-0.05)
  axis(side = 1, at=seq(0, chrsize[chr], 10e+05), labels=rep("", length(seq(0, chrsize[chr], 10e+05)) ), lty=1, col="gray", tck=-0.02)
  text(c(seq(0, chrsize[chr], 10e+06), 11e+07)*1.02, -1.5, 
       labels = c(seq(0, chrsize[chr], 10e+06)/1e+06, "(Mb)"),
       srt = 0, adj = c(1.1,1.1), xpd = TRUE, cex=0.4, srt = 00, col="gray")     
  #
  axis(side = 3, at=seq(0, chrsize[chr], 10e+06), labels=rep("", length(seq(0, chrsize[chr], 10e+06)) ), lty=1, col="gray", tck=-0.05)
  axis(side = 3, at=seq(0, chrsize[chr], 10e+05), labels=rep("", length(seq(0, chrsize[chr], 10e+05)) ), lty=1, col="gray", tck=-0.02)   
  #
  denolg_id <- as.character(unique(df$chrb))
  reflg_id  <- as.character(unique(df$chra))
  mtext(text = reflg_id,  side = 1, line = +0.8, at = 100e+06, cex=0.7, srt=90, font=1, col="gray")
  mtext(text = denolg_id, side = 1, line = -10.5, at = 100e+06, cex=0.7, srt=90, font=1, col="gray")
  # syn over 5 kb
  syn <- df[df$endb-df$startb > 0 & df$type=="SYN", ]
  num = length(syn$chra)
  num = length(syn$chra)
  if(num>0)
    for(i in c(1:num))
    {
      x<-c(syn[i, 5], syn[i, 6], syn[i, 3], syn[i, 2])
      y<-c(+1,          +1,          -1,          -1)  + c(0.05, 0.05, -0.05, -0.05)
      polygon(x, y,
              col = alpha("gray", 1),
              border = NA)
    }
  min_sv_size=100
  # trans over 5 kb
  trans <- df[df$endb-df$startb > min_sv_size & (df$type=="TRANS" | df$type=="INVTR"), ]
  num = length(trans$chra)
  if(num>0)
    for(i in c(1:num))
    {
      x<-c(trans[i, 5], trans[i, 6], trans[i, 3], trans[i, 2])
      y<-c(+1,          +1,          -1,          -1)  + c(0.05, 0.05, -0.05, -0.05)
      polygon(x, y,
              col = alpha("yellowgreen", 0.8),
              border = NA)
    }
  # dup over 5 kb
  dup <- df[df$endb-df$startb > min_sv_size & (df$type=="DUP" | df$type=="INVDP" ), ]
  num = length(dup$chra)
  if(num>0)
    for(i in c(1:num))
    {
      x<-c(dup[i, 5], dup[i, 6], dup[i, 3], dup[i, 2])
      y<-c(+1,          +1,          -1,          -1)  + c(0.05, 0.05, -0.05, -0.05)
      polygon(x, y,
              col = alpha("deepskyblue", 0.5),
              border = NA)
    }
  # inv over 5 kb
  inv <- df[df$endb-df$startb > min_sv_size & (df$type=="INV"), ]
  num = length(inv$chra)
  if(num>0)
    for(i in c(1:num))
    {
      x<-c(inv[i, 5], inv[i, 6], inv[i, 2], inv[i, 3])
      y<-c(+1,          +1,          -1,          -1)  + c(0.05, 0.05, -0.05, -0.05)
      polygon(x, y,
              col = alpha("orange", 0.9),
              border = NA)
    }
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
path<-"/netscratch/dep_mercier/grp_schneeberger/projects/Potato_single_cells/reference_manish_assembled_haplotye_aware/hifiasm_assembly_polish_and_select_v2/s12_LG_wise_assembly_scaffolding_omnic_correction/zStie_1_manual_correction/"
#
# one page
pdf(paste(path, "/s12_St1_DM_comparison_syri_out_plot.pdf", sep=""), family="Helvetica", height=2.5, width=7.08)
par(mfrow=c(1,1), mai = c(0.4, 0.1, 0.3, 0.3)); # margin: bottom, left, top, right
# otava LGs below are at the order of reference LGs 1-12/
#for(reflg in c(6, 11, 4, 9, 12, 3, 8, 10, 5, 1, 2, 7))
for(reflg in c(1:12))
{
  if(reflg==1)
  {
    denovolg=6
    otavalgs <- c("homLG_6_LG_20")
    otavalgpairs <- c('homLG_6_LG_20_vs_chr01') # ac
    parentalinfo <- c()
  }else
    if(reflg==2)
    {
      denovolg=11
      otavalgs <- c("homLG_11_LG_30")
      otavalgpairs <- c('homLG_11_LG_30_vs_chr02') # ac
    }else
      if(reflg==3)
      {
        denovolg=4
        otavalgs <- c("homLG_4_LG_17")
        otavalgpairs <- c('homLG_4_LG_17_vs_chr03') # ac
      }else    
        if(reflg==4)
        {
          denovolg=9
          otavalgs <- c("homLG_9_LG_28")
          otavalgpairs <- c('homLG_9_LG_28_vs_chr04') # ac
        }else    
          if(reflg==5)
          {
            denovolg=12
            otavalgs <- c("homLG_12_LG_1")
            otavalgpairs <- c('homLG_12_LG_1_vs_chr05')  # ac
          }else    
            if(reflg==6)
            {
              denovolg=3
              otavalgs <- c("homLG_3_LG_13")
              otavalgpairs <- c('homLG_3_LG_13_vs_chr06') # ac
            }else    
              if(reflg==7)
              {
                denovolg=8
                otavalgs <- c("homLG_8_LG_37")
                otavalgpairs <- c('homLG_8_LG_37_vs_chr07') # ac
              }else    
                if(reflg==8)
                {
                  denovolg=10
                  otavalgs <- c("homLG_10_LG_3")
                  otavalgpairs <- c('homLG_10_LG_3_vs_chr08') #ac
                }else    
                  if(reflg==9)
                  {
                    denovolg=5
                    otavalgs <- c("homLG_5_LG_29")
                    otavalgpairs <- c('homLG_5_LG_29_vs_chr09') # ac
                  }else    
                    if(reflg==10)
                    {
                      denovolg=1
                      otavalgs <- c("homLG_1_LG_42")
                      otavalgpairs <- c('homLG_1_LG_42_vs_chr10')# ac
                    }else    
                      if(reflg==11)
                      {
                        denovolg=2
                        otavalgs <- c("homLG_2_LG_12")
                        otavalgpairs <- c('homLG_2_LG_12_vs_chr11') # ac
                      }else    
                        if(reflg==12)
                        {
                          denovolg=7
                          otavalgs <- c("homLG_7_LG_24")
                          otavalgpairs <- c('homLG_7_LG_24_vs_chr12') # ac
                        }
  #############################################################################################################################    
  #############################################################################################################################  
  #############################################################################################################################  
  ############# pair 1 ########################################################################################################
  ipair = 1
  alignmraw              <- read.table(paste(path, "/2022_corrected_ALLHiC_build_",otavalgs,"/synteny_check_",otavalgs,"_DH_corrected/", "/syri_selected_minimap2_sorted.txt", sep=""))
  len                    <- length(alignmraw$V1)
  alignm2t               <- cbind(alignmraw[, c(1,2,3,  6,7,8, 9,11)]) # careful on ordering ab
  colnames(alignm2t)     <- c("chra", "starta", "enda", "chrb", "startb", "endb", "id", "type")  
  alignm                 <- alignm2t
  ##
  chrsize                <- max(100000000)
  chr        <- 1
  #############################################################################################################################        
  #############################################################################################################################     
  speciescol <- "mediumorchid3"
  df         <- as.data.frame(alignm)
  pp7        <- my_plot(df, chr, speciescol)
  mtext(text = paste("Reference LG ", reflg, sep=""), side = 3, line = 0.5, at = -chrsize/6, cex = 0.7, col = "gray")  
  #
}
dev.off()











