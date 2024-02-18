# Step 2:
# Here we manually select the best grouping of markers: for C, E, G, for some chrs! 20221229
# find the case where the difference among the haplotypes is minimal
this_date=20221229
#
find_min <- function(data_chr)
{
  this_dim <- dim(data_chr)
  min_diff <- 1e+20
  min_case <- data_chr[1, ]
  for(i in 1:this_dim[1])
  {
    this_line <- data_chr[i, ]
    this_diff <- 0
    skip <- 0
    for(j in 4:6)
    {
      if(this_line[, j] == 0)
      {
        skip = 1
        break
      }
    }
    if(skip == 1)
    {
      next
    }
    for(j in 4:6)
    {
      for(k in (j+1):7)
      {
        this_diff <- this_diff + abs(this_line[, j] - this_line[, k] )
      }
    }
    if(min_diff > this_diff)
    {
      min_diff <- this_diff
      min_case <- this_line
    }
  }
  return(min_case)
}
#
binning_path <- "/netscratch/dep_mercier/grp_schneeberger/projects/Potato_multipleCultivars/s2_10_Cultivars_PacBio_HiFi/a3_hic_alignment/bgi_deep_seq/"
setwd(binning_path)
#
for(sample in c("G"))
{
  data <- read.table(paste(binning_path, "/hic_",sample,"/bwa_align/lg_wise_contigs_read_align/interdiate_res_sample_",sample,"_hic_binining_marker_size_reformat.txt", sep=""))
  #
  chrs <- as.character(unique(data$V1))  
  #
  min_cases <- c()
  for(chr in "chr12")
  {
    # sample "O"
    size_chr <- data[data$V1==chr & data$V2!=65 & data$V3>=2.2 & data$V5<=0.20, ]
    # general
    # size_chr <- data[data$V1==chr & data$V2!=65 & data$V3>=2.2 & data$V5<=0.20, ]
    # read intra-group and inter-group hic stats
    hic_chr <- as.data.frame(matrix(data = NA, nrow = length(size_chr$V1), ncol = 7))
    for (ri in 1:length(size_chr$V1))
    {
        min_hic_contact   <- size_chr[ri, 2]
        min_hic_contact_i <- size_chr[ri, 3]
        min_hapctg_size   <- size_chr[ri, 4]
        max_allelic_ratio <- format(size_chr[ri, 5],nsmall = 2)
        this_hic_stat_file <- paste(binning_path, "/hic_", sample, "/bwa_align/lg_wise_contigs_read_align/hic_binning_", chr, "/s6_",chr,"_", min_hic_contact,"_",min_hic_contact_i,"_", min_hapctg_size,"_", max_allelic_ratio, "_raw_tig_marker_cross_link_count/s8_grouping_window_markers_refined_1st_group_hic_stats.txt", sep="")
        this_hic_stat      <- read.table(this_hic_stat_file)
        if(length(this_hic_stat$V1) != 10) # failure case - no need to keep in further checking
        {
          this_hic_stat    <- this_hic_stat[this_hic_stat$V2<=4 & this_hic_stat$V3<=4, ]
          this_hic_stat$V4 <- 0
          this_hic_stat$V6 <- 0
          cat(paste(chr,"_", min_hic_contact,"_", min_hapctg_size,"_", max_allelic_ratio, " excluded.\n", sep=""))
        }
        # Later, use the sum of intra-group hic divided by the sum of inter-group hic as the indicator for selecting the best phasing, 
        # with some control on size diff of hap chrs.
        avg_intra_grp      <- sum(this_hic_stat[c(1,5,8,10), 4])/4
        avg_inter_grp      <- sum(this_hic_stat[c(2,3,4,6,7,9), 4])/6
        this_hic_stat_46   <- c(avg_intra_grp,
                                avg_inter_grp,
                                round( avg_intra_grp / (avg_inter_grp+0.1), digits = 2) )
        #
        hic_chr[ri, 1]    <- chr
        hic_chr[ri, 2]    <- min_hic_contact
        hic_chr[ri, 3]    <- min_hapctg_size
        hic_chr[ri, 4:6]  <- this_hic_stat_46
        hic_chr[ri, 7]    <- round( sd(as.numeric(size_chr[ri, 6:9]) ) )
    }  
    #
    data_chr <- cbind(size_chr,  hic_chr)
    colnames(data_chr) <- paste("V", 1:16, sep="")
    #
    data_chr_sorted <- data_chr[order(data_chr$V15, decreasing = TRUE), ]
    # rm 
    # clean large variations among chrs
    data_chr_sorted_cleaned <- data_chr_sorted
    to_be_removed <- c()
    for (ci in 1:length(data_chr_sorted$V1) )
    {
      this_line <- data_chr_sorted[ci, ]
      rm_line   <- 0
      for(j in 6:8)
      {
        for(k in (j+1):9)
        {
          this_diff <- abs(this_line[, j] - this_line[, k] )
          # if the size difference between any paired chrs is larger than 25 Mb, or the sd of all chr-sizes is larger than 15 Mb, skip it
          # if(this_diff > 30000000 | data_chr_sorted[ci,16] > 18000000) # C:chr10,E:chr06, G:ch03 only
          if(this_diff > 20000000 | data_chr_sorted[ci,16] > 15000000)
          {
            rm_line <- 1   
            break
          }
        }
        if(rm_line == 1)
        {
          break
        }
      }
      if(rm_line == 1)
      {
        #print(paste("cleaning ", ci, " row", sep=""))
        to_be_removed <- c(to_be_removed, ci)
      }
    }
    data_chr_sorted_cleaned <- data_chr_sorted_cleaned[-1*to_be_removed, ]   
    data_chr_sorted_cleaned <- data_chr_sorted_cleaned[data_chr_sorted_cleaned$V15>0, ]
    # sort by allelic - smaller better (lower chance to mis-join contigs into linkage groups)
    tmp_x <- data_chr_sorted_cleaned[order(data_chr_sorted_cleaned$V5, decreasing=FALSE), c(1:9, 15, 16)]
    # tmp_x <- tmp_x[1:min(150, length(tmp_x$V1)), ]         
    # # sort by hic-strength - larger better - backbone construction
    # # tmp_x <- data_chr_sorted_cleaned[order(data_chr_sorted_cleaned$V2, decreasing=TRUE), ]
    # tmp_x <- tmp_x[order(tmp_x$V2, decreasing=TRUE), ]
    # tmp_x <- tmp_x[1:min(150, length(tmp_x$V1)), ]     
    # # sort by hic-strength - larger better - integrate
    # tmp_x <- tmp_x[order(tmp_x$V3, decreasing=TRUE), ]
    # tmp_x <- tmp_x[1:min(150, length(tmp_x$V1)), ] 
    # # sort by intra-group hic / inter-group hic - larger better
    # tmp_x <- tmp_x[order(tmp_x$V15, decreasing=TRUE), ]
    # rm done
    colnames(tmp_x) <- c("chr", "min_hic_contact", "min_hic_contact_i", "min_hapctg_size", "max_allelic_ratio", "h1_size", "h2_size", "h3_size", "h4_size", "intra_inter_hic_ratio", "hap_size_sd")    
    # select one of the following sort after visually checking the tmp_x
    tmp_x <- tmp_x[order(tmp_x$hap_size_sd), ]
    tmp_x <- tmp_x[order(tmp_x$intra_inter_hic_ratio, decreasing=TRUE), ]
    tmp_x <- tmp_x[order(tmp_x$max_allelic_ratio, decreasing=FALSE), ]
    tmp_x <- tmp_x[order(tmp_x$min_hic_contact, decreasing=TRUE), ]
    # below is for your selection settings
    tmp_x <- tmp_x[tmp_x$hap_size_sd<5.1e+06 & tmp_x$intra_inter_hic_ratio>4.0 & tmp_x$max_allelic_ratio<0.03 & tmp_x$min_hic_contact>30, ]
    #
    min_case  <- find_min( tmp_x[1:min(1, length(tmp_x$chr)), ] ) 
    colnames(min_case)=NULL
    print(min_case)     
  }
}