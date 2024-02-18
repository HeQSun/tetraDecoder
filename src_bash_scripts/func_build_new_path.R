# Here we figure out a new assembly with the help of ragoo, for chrs where all contigs are pieces.
# 
## get options from cmd
###############################################################################
args<-commandArgs(TRUE)
if(length(args) < 4) {
  cat(length(args))
  cat("\n   Info: rebuild the assembly tour in a format of juicer tool.\n")
  cat("   USAGE: Rscript func_build_new_path.R 1.sample 2.lg 3.real_hapi 4.path\n\n")
}else
{
  ## check args
  sample    = args[[1]]
  lg        = args[[2]]
  real_hapi = args[[3]]  
  path      = args[[4]]
  #setwd("/netscratch/dep_mercier/grp_schneeberger/projects/Potato_multipleCultivars/s2_10_Cultivars_PacBio_HiFi/zaux_ragoo/ragoo_D/hap_1/")
  setwd(path)
  #
  # given
  assembly<-read.table(paste("/netscratch/dep_mercier/grp_schneeberger/projects/Potato_multipleCultivars/s2_10_Cultivars_PacBio_HiFi/a6_lg_wise_scaffolding/scaf_",sample,"/s_read_realign_",lg,"/zcorrection_four_hap_juicebox_DM_strand_interhap_correction/updated/groups.asm.hap",real_hapi,"_strand_raw.review.oneseq.updated.assembly.ctg.info.only", sep=""))
  # ragoo
  ragoo_out<-read.table(paste(path, "/ragoo_output/orderings/",lg,"_orderings.txt", sep=""))
  ragoo_out_new <- as.data.frame(ragoo_out)
  # reorder the assembly according to ragoo assembly
  # reorder the assembly according to ragoo assembly
  ctg_num  <- length(ragoo_out$V1)
  assembly <- as.data.frame(assembly[1:ctg_num, ])
  for(ri in 1:ctg_num)
  {
    this_ctg             <- as.character(paste(">", as.character(ragoo_out[ri, 1]), sep=""))
    this_ctg_info        <- assembly[assembly[, 1]==this_ctg, ]
    if(length(this_ctg_info$V1)>0)
    {
      ragoo_out_new[ri, 3] <- this_ctg_info[, 2]
      ragoo_out_new[ri, 4] <- this_ctg_info[, 3]
    }else
    {
      ragoo_out_new[ri, 3] <- this_ctg
      ragoo_out_new[ri, 4] <- this_ctg 
    }
  }
  # format of new ragoo_out: >ctg_str_id strand ctg_int_id ctg_size
  tour <- c()
  for(ri in 1:ctg_num)
  {
    this_ctg_info <- ragoo_out_new[ri, ]
    if(this_ctg_info$V2=="+")
    {
      tour <- c(tour, paste("", this_ctg_info$V3, sep="") )
    }else
    {
      tour <- c(tour, paste("-", this_ctg_info$V3, sep="") )
    }
  }
  write.table(x = paste(tour, collapse=" "), file=paste(sample, "_", lg, "_hap", real_hapi, "_updated_from_ragoo.tour", sep=""), append = F, quote = F, row.names = F, col.names = F)
}
