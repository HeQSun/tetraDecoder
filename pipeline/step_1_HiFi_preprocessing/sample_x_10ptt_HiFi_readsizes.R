path <-"/your/working/path/"

for (sample in  LETTERS[8])
{
y<-as.data.frame(read.table(paste(path, "/sample_",sample, "/sample_", sample,"_zcat_HiFi_length_4cells.txt", sep="")))
x<-y$V1
pdf(paste(path, "/sample_",sample, "/sample_", sample,"_zcat_HiFi_length_4cells.pdf", sep=""), height=4, width=8.26)
par(mai = c(1.0, 1.0, 0.5, 0.5)); # margin: bottom, left, top, right
h<-hist(x,
     breaks = 500,
     xlim   = c(0, 5e+04),
     xlab   = "Read length",
     ylab   = "Frequency",
#    main   = "raw pb read length hist\n(1 kb bin)",
     main   = " ",
     col    = "red",
     border = NA)
abline(v = 15000, lty=2, col="gray")
maxindex <- which.max(h$counts)
peakvalue <- h$mids[maxindex]
legend("topright",
       pch    = c(15,15,15,15,15,15,15),
       legend = c(paste("All ", sum(x>0e+04), " reads: ", round(sum(x[x>0]  /  1e+09), digits=2), " Gb", sep=""),
                  paste(">10 kb ", sum(x>1e+04), " reads: ", round(sum(x[x>1e+04]/1e+09), digits=2), " Gb", sep=""),
                  paste(">20 kb ", sum(x>2e+04), " reads: ", round(sum(x[x>2e+04]/1e+09), digits=2), " Gb", sep=""),
                  paste("Peak (over 10kb): ", round(peakvalue/1000), " kb", sep=""),
                  paste("Median (overall): ", round(median(x)/1000), " kb", sep=""),
                  paste("Mean (overall): ", round(mean(x)/1000),   " kb", sep="")
                 ),
       horiz  = FALSE,
       border = "NA",
       bty    = "n",
       cex    = 0.8)
dev.off()


# non trimmed

# y<-as.data.frame(read.table(paste(path, "/sample_",sample, "/sample_", sample,"_zcat_HiFi_length_4cells_nontrimmed.txt", sep="")))
# x<-y$V1
# pdf(paste(path, "/sample_",sample, "/sample_", sample,"_zcat_HiFi_length_4cells_nontrimmed.pdf", sep=""), height=4, width=8.26)
# par(mai = c(1.0, 1.0, 0.5, 0.5)); # margin: bottom, left, top, right
# h<-hist(x,
#         breaks = 500,
#         xlim   = c(0, 5e+04),
#         xlab   = "Read length",
#         ylab   = "Frequency",
#         #    main   = "raw pb read length hist\n(1 kb bin)",
#         main   = " ",
#         col    = "red",
#         border = NA)
# abline(v = 15000, lty=2, col="gray")
# maxindex <- which.max(h$counts)
# peakvalue <- h$mids[maxindex]
# legend("topright",
#        pch    = c(15,15,15,15,15,15,15),
#        legend = c(paste("All ", sum(x>0e+04), " reads: ", round(sum(x[x>0]  /  1e+09), digits=2), " Gb", sep=""),
#                   paste(">10 kb ", sum(x>1e+04), " reads: ", round(sum(x[x>1e+04]/1e+09), digits=2), " Gb", sep=""),
#                   paste(">20 kb ", sum(x>2e+04), " reads: ", round(sum(x[x>2e+04]/1e+09), digits=2), " Gb", sep=""),
#                   paste("Peak (over 10kb): ", round(peakvalue/1000), " kb", sep=""),
#                   paste("Median (overall): ", round(median(x)/1000), " kb", sep=""),
#                   paste("Mean (overall): ", round(mean(x)/1000),   " kb", sep="")
#        ),
#        horiz  = FALSE,
#        border = "NA",
#        bty    = "n",
#        cex    = 0.8)
# dev.off()
}
