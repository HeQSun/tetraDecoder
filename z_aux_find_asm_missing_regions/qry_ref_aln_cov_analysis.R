# Here we check those regions along a ref seq that are not covered by a query sequence.
# functions
merge_overlap_aln <- function(aln)
{
  # merge overlapping alignments, according to ref seq coordinates.
  changed <- 1
  while(changed != 0)
  {
    this_dim <- dim(aln)
    changed  <- 0
    for(i in 2:this_dim[1])
    {
      last_interval <- aln[i-1, ]
      this_interval <- aln[i, ]
      s1 <- last_interval$rsta
      e1 <- last_interval$rend
      s2 <- this_interval$rsta
      e2 <- this_interval$rend
      if(s2>=s1 & s2<=e1)
      {
        #cat("clean ", paste(aln[i, 4:7], sep=" "), " ovlp ", paste(aln[i-1, 4:7], sep=" "), "\n")
        if(e2>=e1)
        {
          aln[i-1, ]$rend <- e2
          aln[i-1, ]$qend <- this_interval$qend
        }
        aln <- aln[-1*i, ]
        changed <- 1
        break
      }
    }
  }
  return(aln)
}
#
find_mis_region<-function(aln, min_gap_size, min_fold)
{
  # find intervals with large gaps along ref sequence
  this_dim <- dim(aln)
  interval_involved <- c()
  for(i in 2:this_dim[1])
  {
    last_interval <- aln[i-1, ]
    this_interval <- aln[i, ]
    s1 <- last_interval$rsta
    e1 <- last_interval$rend
    s2 <- this_interval$rsta
    e2 <- this_interval$rend
    if(s2 - e1 > min_gap_size)
    {
      #cat(paste(aln[i-1, 4:7], sep=" "), " with ", paste(aln[i, 4:7], sep=" "),  " with gap ", s2 - e1, " bp", "\n")
      interval_involved <- c(interval_involved, i-1, i)
    }
  }
  interval_involved <- unique(interval_involved)
  return(aln[interval_involved, ])
}
#
sample="A"
lg="chr01"
real_hapi=4
#
for (sample in c(LETTERS[1:10], "O"))
#for (sample in "A")
{
  this_case <- read.table(paste("/netscratch/dep_mercier/grp_schneeberger/projects/Potato_multipleCultivars/s2_10_Cultivars_PacBio_HiFi/a6_lg_wise_scaffolding/scaf_",sample,"/",sample,"_juicebox_corr_manual_strand_selection.txt", sep=""))
  ####this_case <- read.table("/Users/sun/Desktop/z_10potato_project/xjtu_a9_find_regions_not_assembled/src_find_gaps/R_cov_check/A_juicebox_corr_manual_strand_selection.txt")
  #
  for(casei in 1:length(this_case$V1))
  {
    #
    lg        <- as.character(this_case[casei, 1])
    real_hapi <- this_case[casei, 2]
    #
    ####aln_raw <- as.data.frame(read.table(paste("/Users/sun/Desktop/z_10potato_project/xjtu_a9_find_regions_not_assembled/src_find_gaps/data/fw_data_",sample,"_",lg,"_", real_hapi, sep="")))
    aln_raw <- as.data.frame(read.table(paste("/netscratch/dep_mercier/grp_schneeberger/projects/Potato_multipleCultivars/s2_10_Cultivars_PacBio_HiFi/a6_lg_wise_scaffolding/scaf_",sample,"/s_read_realign_",lg,"/zcorrection_four_hap_juicebox_DM_strand_interhap_correction/zcorrection_",real_hapi,"_phase_corr_4th_refDM/fw_data", sep="")))
    colnames(aln_raw) <-c("qchr", "qsta", "qend", "strand", "rchr", "rsta", "rend")
    # merge overlapping intervals
    aln <- merge_overlap_aln(aln_raw)
    #
    min_gap_size = 500000 # size of gap in ref
    min_fold     = 5      # size of gap in ref is 5 times of gap in query
    # find large gaps
    gap_involved_interal <- find_mis_region(aln, min_gap_size, min_fold)
    # output involved intervals
    print(gap_involved_interal)
    write.table(x = gap_involved_interal, file = paste("/netscratch/dep_mercier/grp_schneeberger/projects/Potato_multipleCultivars/s2_10_Cultivars_PacBio_HiFi/a6_lg_wise_scaffolding/scaf_",sample,"/s_read_realign_",lg,"/zcorrection_four_hap_juicebox_DM_strand_interhap_correction/zcorrection_",real_hapi,"_phase_corr_4th_refDM/fw_data_potential_missing", sep=""), append = F, quote = F, row.names = F)
  }
}





