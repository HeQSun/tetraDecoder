
# .libPaths("/dss/dsshome1/lxc03/di36guz2/R")
# mamba activate env_others
.libPaths()

setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/03_Nhaplotypes")

library(ggplot2)
library(tidyverse)

chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

DM_ref_chrSizes <- data.frame(chr = c("chr01", "chr02", "chr03", "chr04", "chr05", "chr06", "chr07", "chr08", "chr09", "chr10", "chr11", "chr12"), 
                              totalLen = c(88591686, 46102915, 60707570, 69236331, 55599697, 59091578, 57639317, 59226000, 67600300, 61044151, 46777387, 59670755))

chrN <- "01"

for (chrN in chr_list) {
  #lengthHaplBlocks_all = data.frame(chr="tem", haplotype1="tem", haplotype2="tem", startpos=0, endpos=0)
  table_hapG <- read_table(paste0("haplotypeGroups_min10SS_chr",chrN,".txt"), F)
  names(table_hapG) <- c("pos", "sample", "hapGroup")
  
  table_hapG_spread <- table_hapG %>%
    spread(key=sample, value=hapGroup) %>%
    arrange(pos)
  
  
  list_samples <- unique(table_hapG$sample)
  
  for (i in seq(1,length(list_samples))){
    for (j in seq(1,length(list_samples))){
      if (j > i){
        table_hapG_2Hapl <- table_hapG_spread %>%
          select(pos, list_samples[i], list_samples[j]) 
        
        sameHapl <- as.vector(((table_hapG_2Hapl[,2] == table_hapG_2Hapl[,3])*1))
        
        lengthHaplBlocks <- table_hapG_2Hapl %>%
          mutate(sameHapl) %>%
          #mutate(cont_haplo=(sameHapl==lag(sameHapl) & sameHapl==1)*1) 
          mutate(cont_haplo=(sameHapl + lag(sameHapl) > 0)*1) %>%
          mutate(sequence = data.table::rleid(cont_haplo == 1)) %>%
          filter(cont_haplo == 1) %>%
          group_by(sequence) %>%
          summarise(startpos = min(pos), 
                    endpos=max(pos)) %>%
          mutate(chr=paste0("chr", chrN), 
                 haplotype1=list_samples[i], 
                 haplotype2=list_samples[j]) %>%
          select(chr, haplotype1, haplotype2, startpos, endpos)
        #mutate(len_haploBlock=endpos-startpos) %>%
        #pull(len_haploBlock)
        #ggplot(aes(len_haploBlock)) +
        #geom_histogram()
        #head()
        write.table(lengthHaplBlocks, file=paste0("Table_SharedBlocks_chr", chrN, ".csv"), quote = F, row.names = F, col.names = F, sep = ",", append=T)
      }
    }
  }
}



chrN="10"
for (chrN in chr_list) {
lengthHaplBlocks_all <- read.csv(paste0("Table_SharedBlocks_chr",chrN,".csv"), F)
names(lengthHaplBlocks_all) <- c("chr", "haplotype1", "haplotype2", "startpos", "endpos")

lengthHaplBlocks_all %>%
  mutate(lenHaplBlocks=endpos-startpos) %>%
  group_by(chr, haplotype1, haplotype2) %>%
  summarise(NSharedBlocks=n(), 
            totalLenShared=sum(lenHaplBlocks)) %>%
  left_join(DM_ref_chrSizes, by="chr") %>%
  mutate(frac=totalLenShared/totalLen) %>%
  write.table(file=paste0("Table_PropSharedBlocks.csv"), quote = F, row.names = F, col.names = F, sep = ",", append=T)
}




lengthHaplBlocks_all <- read.csv(paste0("Table_PropSharedBlocks.csv"), F)
names(lengthHaplBlocks_all) <- c("chr", "haplotype1", "haplotype2", "NSharedBlocks", "totalLenShared", "totalLen", "frac")


lengthHaplBlocks_all_meanValues <- lengthHaplBlocks_all %>%
  group_by(chr) %>%
  summarise(mean_frac=mean(frac), 
            min_frac=min(frac), 
            max_frac=max(frac)) 


#   chr   mean_frac min_frac max_frac
#  Chr01     0.204   0.0672    0.850
#  Chr02     0.179   0.0243    0.773
#  Chr03     0.212   0.0693    0.621
#  Chr04     0.181   0.0391    0.794
#  Chr05     0.258   0.0487    0.855
#  Chr06     0.197   0.0425    0.743
#  Chr07     0.190   0.0231    0.938
#  Chr08     0.189   0.0474    0.830
#  Chr09     0.196   0.0251    0.901
#  Chr10     0.348   0.0429    0.800
#  Chr11     0.195   0.0624    0.880
#  Chr12     0.285   0.0295    0.875


lengthHaplBlocks_all %>%
  ggplot(aes(frac))+
  geom_histogram(fill="black")+
  geom_vline(data=lengthHaplBlocks_all_meanValues, aes(xintercept = mean_frac), 
             #linetype="dotted", 
             color = "red") +
  labs(x = "Fraction",
       y = "Count") +
  facet_grid(chr~., scale="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))

ggsave(paste0("09_MeanShareHaplBlocks_10kbWin_min10SS.png"), height=10, width=3)




lengthHaplBlocks_all %>%
  ggplot(aes(frac))+
  geom_histogram(fill="black")+
  geom_vline(data=lengthHaplBlocks_all_meanValues, aes(xintercept = mean_frac), 
             #linetype="dotted", 
             color = "red") +
  scale_x_continuous(breaks=seq(0,1,0.1))+
  labs(x = "Fraction",
       y = "Count") +
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor = element_line(color = "grey80"),
        axis.text = element_text(size=24), 
        axis.title = element_text(size=24))

ggsave(paste0("09_MeanShareHaplBlocks_ed_10kbWin_min10SS.png"), height=5, width=7)
