# here we check depth distribution along initially assembled contigs: hifiasm version 0.7
#
library("scales")
#
for (sample in LETTERS[c(1,2,3,4,5,6,7,8,9,10,15)])
#for (sample in LETTERS[c(2,3,4,6,7,8,9)])
#for (sample in LETTERS[c(10)])
{
        subset=4
        subfolder=paste("subset",subset,"_illu_re_align_purge_ovl/DNB_align", sep="")
        #subfolder=paste("subset",subset,"_hifi_re_align_purge_dup_ctg/DNB_align", sep="")
        #
        path <- paste("/netscratch/dep_mercier/grp_schneeberger/projects/Potato_multipleCultivars/s2_10_Cultivars_PacBio_HiFi/a2_initial_assembly/sample_",sample,"/hifiasm_asm_v0p7/",subfolder,"/", sep="")
        setwd(path)
        #
        pdf(paste(path, "sample_",sample,"_cnv_initial_contigs_winsize10kb_step10kb_depth.pdf", sep=""), family="Helvetica", height= 3, width=7.08661)
        #
        par(mai = c(0.7, 0.8, 0.25, 0.1)) # margin: bottom, left, top, right
        winsize  <- 10000
        winstep  <- winsize
        # raw before any purge
        path_raw <- paste("/netscratch/dep_mercier/grp_schneeberger/projects/Potato_multipleCultivars/s2_10_Cultivars_PacBio_HiFi/a2_initial_assembly/sample_",sample,"/hifiasm_asm_v0p7/z_DNB_no_purge_depth_illu/", sep="")        
        win_raw  <- read.table(paste(path_raw, "/cnv_winsize", winsize, "_step", winstep,"_hq.txt", sep=""))
        if(sample=="O")
        {
                win_raw$V6 <- win_raw$V6*1.08 # one with hifi, raw not, so I am scaling raw a bit here!
        }        
        #
        path_win <- path
        win3     <- read.table(paste(path_win, "/cnv_winsize", winsize, "_step", winstep,"_hq.txt", sep=""))
        # test 
        if(0)
        {
           # test only
           low_cov_list <- read.table("/netscratch/dep_mercier/grp_schneeberger/projects/Potato_multipleCultivars/s2_10_Cultivars_PacBio_HiFi/a2_initial_assembly/sample_C/hifiasm_asm/subset1_illu_re_align_purge_ovl/DNB_align/low_cov_ctgs_to_excluded.list")
           win2 <- win3[!(win3$V1 %in% low_cov_list$V1), ]
        }else
        {
                winx <- win3[win3$V3-win3$V2>5000, ]
                winl <- win3[win3$V6<30, ]
                win2 <- win3
        }
        #
        for(exlen in c(0))
        {        
                # select non-end windows
                exlen_end <- exlen
                if(exlen > 0 )
                {
                        win     <- win2[(win2$V2>1 & win2$V7-win2$V3>exlen_end), ]
                }else
                {
                        win         <- win2
                }        
                v_scale  <- 113 # genome-wide average
                win$V5   <- round(win$V6 / v_scale * v_scale)
                win_keep <- win
                #
                # set up break number
                breaknum <- round(max(win$V6))
                hwin_raw <- hist(win$V6, breaks=breaknum, plot=F) # raw depth
                hwin_nor <- hist(win$V5, breaks=breaknum, plot=F) # normalized with genome-wide average
                ylim_both<- max(hwin_nor$counts, hwin_raw$counts)*1.3
                if(winsize==10000) ylim_both = 20000
                #
                yle40<-win[win$V6<40, ]
                yge40<-win[win$V6>=40, ]
                yle40_len <- sum(yle40$V3-yle40$V2+1)
                yge40_len <- sum(yge40$V3-yge40$V2+1)
                yle40_len/(yle40_len+yge40_len)
                # plot
                hwin_raw <- hist(win$V6,
                                 breaks = breaknum,
                                 col    = alpha("orangered", 0.8),
                                 xlim   = c(0, 1000),
                                 ylim   = c(0, ylim_both),
                                 border = NA,
                                 plot   = T,
                                 xlab   = "",
                                 ylab   = "",
                                 axes = FALSE, 
                                 main   = paste("Depth at win ", winsize, 
                                                " bp (step: ", winstep," bp)", sep=""),
                                 font.main=1, # not bold
                                 cex.main=1
                )
                cnt_max_index <- which.max(hwin_raw$counts)
                abline(v=hwin_raw$mids[cnt_max_index], col="gray", lty=2)
                text(x = hwin_raw$mids[cnt_max_index]-35, y=ylim_both*0.8, 
                     labels = paste(hwin_raw$mids[cnt_max_index], "x", sep=""), col="gray", cex=0.7)                
                if(winsize==10000) write(x = cnt_max_index, file = paste(path, "sample_",sample, "_avg_cov", sep=""), append = F)
                cat("sample_",sample, "_avg_cov", cnt_max_index, "\n")
                if(exlen==0)
                {
                   # this_gse <- round(sum( (win$V3-win$V2+1)*win$V6/cnt_max_index/1e+09 ), digits = 3)
                   # cat("gse=",this_gse, "\n")
                }
                # hap size
                win_hap      <- win[win$V6<=1.7*cnt_max_index, ]
                win_hap_size <- round(sum( (win_hap$V3-win_hap$V2+1)*win_hap$V6/cnt_max_index/1e+09), digits = 3)                
                # non-hap size
                win_non_hap      <- win[win$V6>1.7*cnt_max_index, ]
                win_non_hap_size <- round(sum( (win_non_hap$V3-win_non_hap$V2+1)*win_non_hap$V6/cnt_max_index/1e+09), digits = 2)
                # [0, 204.8]
                hap_upper        <- 1.53*cnt_max_index
                win_hap          <- win[win$V6<=hap_upper, ]
                win_hap_size     <- sum(win_hap$V3-win_hap$V2+1)
                # (204.8, 332.8]
                dip_upper        <- 2.45*cnt_max_index
                win_dip          <- win[win$V6>hap_upper & win$V6<=dip_upper, ]
                win_dip_size     <- sum(win_dip$V3-win_dip$V2+1)*2
                # (332.8, 460.8]
                trip_upper       <- 3.45*cnt_max_index
                win_trip         <- win[win$V6>dip_upper & win$V6<=trip_upper, ]
                win_trip_size    <- sum(win_trip$V3-win_trip$V2+1)*3
                # (460.8, 588.8]
                tetrap_upper     <- 4.45*cnt_max_index
                win_tetrap       <- win[win$V6>trip_upper & win$V6<=tetrap_upper, ]
                win_tetrap_size  <- sum(win_tetrap$V3-win_tetrap$V2+1)*4
                # (588.8, --)
                win_rep          <- win[win$V6>tetrap_upper, ]
                win_rep_size     <- round(sum( (win_rep$V3-win_rep$V2+1)*win_rep$V6/cnt_max_index ))
                #
                asm_size <- sum(win$V3-win$V2+1)
                this_gse <- round((win_hap_size+win_dip_size+win_trip_size+win_tetrap_size+win_rep_size)/1e+09, digits = 3)
                #
                # ratio of type-contigs in the assembly
                hap_ratio    <- win_hap_size/asm_size
                dip_ratio    <- win_dip_size/asm_size/2
                trip_ratio   <- win_trip_size/asm_size/3
                tetrap_ratio <- win_tetrap_size/asm_size/4
                rep_ratio    <- sum((win_rep$V3-win_rep$V2+1))/asm_size                
                hap_ratio+dip_ratio+trip_ratio+tetrap_ratio+rep_ratio
                #
                #abline(v=180, col="gray", lty=2)
                abline(v=hap_upper, col="gray", lty=2)
                abline(v=dip_upper, col="gray", lty=2)
                abline(v=trip_upper, col="gray", lty=2)
                abline(v=tetrap_upper, col="gray", lty=2)  
                #
                text(x = cnt_max_index*1+34, y=10000, labels = paste(round(hap_ratio*100, digits = 1), "%", sep=""),    col="black", cex=0.7) 
                text(x = cnt_max_index*2+00, y=2000,  labels = paste(round(dip_ratio*100, digits = 1), "%", sep=""),    col="black", cex=0.7) 
                text(x = cnt_max_index*3+00, y=2000,  labels = paste(round(trip_ratio*100, digits = 1), "%", sep=""),   col="black", cex=0.7)
                text(x = cnt_max_index*4+00, y=2000,  labels = paste(round(tetrap_ratio*100, digits = 1), "%", sep=""), col="black", cex=0.7)
                text(x = cnt_max_index*5+00, y=2000,  labels = paste(round(rep_ratio*100, digits = 1), "%", sep=""),    col="black", cex=0.7)
                tig_ratio <- c(hap_ratio, 
                               dip_ratio,
                               trip_ratio,
                               tetrap_ratio,
                               rep_ratio,
                               round(asm_size/1e+09, digits = 3),
                               this_gse)
                #write(x = round(tig_ratio, digits = 6), file = paste(path, "sample_",sample, "_tig_ratio", sep=""), append = F)
                options(scipen=999)
                cat( c(sample, round(tig_ratio, digits = 4)), file = paste(path, "sample_",sample, "_tig_ratio", sep=""))
                #
                # axis
                axis(1, at=seq(0, 1000, 125), label=seq(0, 1000, 125), cex.axis=1)
                axis(2, at=seq(0, 20000, 10000),  label=c("0", "10k", "20k"), cex.axis=1, line = -0.5)
                # axis labels
                mtext(text = "10 kb-window count", side = 2, line = 2.0, at = 8200, cex=1.0, srt=90, col="black")
                mtext(text = "Sequencing depth (x)",        side = 1, line = 2.2, at =  500, cex=1, srt=90, col="black")
                #
                mtext("", side = 3, cex = 0.7, at= -120, line=0.8, adj = 0, family = "Helvetica", font=2)
                #
                # plot
                breaknum <- round(max(win_raw$V6))
                hwin_b4_purge <- hist(win_raw$V6,
                                 breaks = breaknum,
                                 plot   = F)
                lines(hwin_b4_purge$mids, hwin_b4_purge$counts, col="black", type="l", lwd=0.2)
                # raw gs
                cnt_max_index_raw <- which.max(hwin_b4_purge$counts)
                if(exlen==0)
                {
                    this_gse_raw      <- round(sum( (win_raw$V3-win_raw$V2+1)*win_raw$V6/cnt_max_index_raw/1e+09 ), digits = 3)
                    cat("gse_raw=", this_gse_raw, "\n")    
                }
                # legend
                legend("topright",
                       pch    = c(15, 15, 15 ,15,15,15,15,15),
                       col    = c("black", "orangered", "gray", "gray", "gray", "gray", "gray", "gray"),
                       legend = c(paste("Before purging: ", round(sum(win_raw$V3-win_raw$V2+1)/1e+09, digits=2), 
                                        " Gb", sep=""),
                                  paste("After purging (DNB): ", round(sum(win$V3-win$V2+1)/1e+09, digits=2), 
                                        " Gb", sep=""),
                                  paste("Genome size est.: ", this_gse, " Gb", sep=""),
                                  paste("hap: ", win_hap_size, " bp", sep=""),
                                  paste("dip: ", win_dip_size, " bp", sep=""),
                                  paste("trip: ", win_trip_size, " bp", sep=""),
                                  paste("tetrap: ", win_tetrap_size, " bp", sep=""),
                                  paste("rep: ", win_rep_size, " bp", sep="")  ),
                       horiz  = FALSE,
                       border = "NA",
                       bty    = "n",
                       cex    = 0.7)                      
        }
        dev.off()
        rm(list=ls()) 
}

