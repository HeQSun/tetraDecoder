# .libPaths("/dss/dsshome1/lxc03/di36guz2/R")
# mamba activate env_others
.libPaths()

setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/03_Nhaplotypes")

library(ggplot2)
library(tidyverse)

secondHighest <-  function(x) {
  u <- unique(x)
  sort(u, decreasing = TRUE)[2L]
}

centromere_table <- read.csv("centromeres.csv", header = T)

chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

chrID <- "01"

for (chrID in chr_list) {
  
  table_hapG <- read_table(paste0("haplotypeGroups_min10SS_chr",chrID,".txt"), F)
  names(table_hapG) <- c("pos", "sample", "hapGroup")
  
  
  table_hapG_sharedVal <- table_hapG %>%
    group_by(pos, hapGroup) %>%
    summarise(NShared=n()) 
  

  table_hapG %>%
    left_join(table_hapG_sharedVal, by=c("pos", "hapGroup")) %>%
    ggplot(aes(pos/1000000, sample)) +
    geom_tile(aes(fill=NShared))+
    #scico::scale_fill_scico(palette = "lajolla")+
    scale_fill_viridis_c(option = "magma", direction = -1, limits = c(1, 44))+
    scale_x_continuous(breaks=seq(0,100,5), name = "Position (Mb)", expand = c(0,0))+
    scale_y_discrete(name = "Haplotype", expand = c(0,0))+
    ggtitle(paste0("Chr:", chrID))+
    theme_classic()+
    theme(panel.grid.major = element_line(color = "grey80"),
          panel.grid.minor = element_line(color = "grey80"))


  # ggsave(paste0("06_haplotypeSharing_per10kbWin_chr", chr,".png"), height=7, width=7)

  WinSizeUsed = 100000
  table_hapG %>%
    left_join(table_hapG_sharedVal, by=c("pos", "hapGroup")) %>%
    mutate(win100=floor(pos/WinSizeUsed)) %>%
    group_by(sample, win100) %>%
    summarise(MedianNShared=mean(NShared)) %>%
    ggplot(aes(win100*WinSizeUsed/1000000, sample)) +
    geom_tile(aes(fill=MedianNShared))+
    #scico::scale_fill_scico(palette = "lajolla")+
    scale_fill_viridis_c(option = "magma", direction = -1, limits = c(1, 44))+
    scale_x_continuous(breaks=seq(0,100,5), name = "Position (Mb)", expand = c(0,0))+
    scale_y_discrete(name = "Haplotype", expand = c(0,0))+
    ggtitle(paste0("Chr:", chrID))+
    theme_classic()+
    theme(panel.grid.major = element_line(color = "grey80"),
          panel.grid.minor = element_line(color = "grey80"))


  # ggsave(paste0("06b_100kbmedian_haplotypeSharing_per10kbWin_chr", chr,".png"), height=7, width=7)


  WinSizeUsed = 1000000

  table_hapG %>%
    left_join(table_hapG_sharedVal, by=c("pos", "hapGroup")) %>%
    mutate(win100=floor(pos/WinSizeUsed)) %>%
    group_by(sample, win100) %>%
    summarise(MedianNShared=mean(NShared)) %>%
    ggplot(aes(win100*WinSizeUsed/1000000, sample)) +
    geom_tile(aes(fill=MedianNShared))+
    #scico::scale_fill_scico(palette = "lajolla")+
    scale_fill_viridis_c(option = "magma", direction = -1, limits = c(1, 44))+
    scale_x_continuous(breaks=seq(0,100,5), name = "Position (Mb)", expand = c(0,0))+
    scale_y_discrete(name = "Haplotype", expand = c(0,0))+
    ggtitle(paste0("Chr:", chrID))+
    theme_classic()+
    theme_classic()+
    theme(panel.grid.major = element_line(color = "grey80"),
          panel.grid.minor = element_line(color = "grey80"),
          axis.text = element_text(size=18),
          axis.title = element_text(size=18))


  # ggsave(paste0("06c_1Mbmedian_haplotypeSharing_per10kbWin_chr", chr,".png"), height=12, width=15*((max(table_hapG$pos)/1000000)/100))
  
  WinSizeUsed = 1000000
  
  table_hapG %>%
    left_join(table_hapG_sharedVal, by=c("pos", "hapGroup")) %>%
    mutate(win100=floor(pos/WinSizeUsed)) %>%
    group_by(sample, win100) %>%
    summarise(MedianNShared=mean(NShared)) %>%
    ggplot(aes(win100*WinSizeUsed/1000000, sample)) +
    geom_tile(aes(fill=MedianNShared))+
    #scico::scale_fill_scico(palette = "lajolla")+
    scale_fill_viridis_c(option = "magma", direction = -1, limits = c(1, 44))+
    scale_x_continuous(breaks=seq(0,100,5), name = "Position (Mb)", expand = c(0,0))+
    scale_y_discrete(name = "Haplotype", expand = c(0,0))+
    ggtitle(paste0("Chr:", chrID))+
    theme_classic()+
    theme(panel.grid.major = element_line(color = "grey80"), 
          panel.grid.minor = element_line(color = "grey80"),
          axis.text.x = element_text(size=26), 
          axis.text.y = element_blank(), 
          axis.title = element_text(size=26), 
          legend.position = "bottom")
  
  
  # ggsave(paste0("06d_1Mbmedian_haplotypeSharing_per10kbWin_chr", chr,".png"), height=12, width=15*((max(table_hapG$pos)/1000000)/100))
  

  table_hapG %>%
    group_by(pos, hapGroup) %>%
    summarise(NShared=n()) %>%
    group_by(pos) %>%
    summarise(MaxNShared=max(NShared)) %>%
    ggplot(aes(pos/1000000, MaxNShared))+
    geom_line(alpha=0.05)+
    geom_point(alpha=0.1)+
    scale_x_continuous(breaks=seq(0,100,5), name = "Position (Mb)", expand = c(0,0))+
    ggtitle(paste0("Chr:", chrID))+
    theme_classic()+
    theme(panel.grid.major = element_line(color = "grey80"),
          panel.grid.minor = element_line(color = "grey80"))

  # ggsave(paste0("07_MaxhaplotypeSharing_per10kbWin_chr", chr,".png"), height=5, width=12)

  WinSizeU=500000
  if (chrID == "01"){
    merged_table <- table_hapG %>%
      group_by(pos, hapGroup) %>%
      summarise(NShared=n()) %>%
      group_by(pos) %>%
      summarise(MaxNShared=max(NShared),
                SecondMaxNShared=secondHighest(NShared),
                NHaplotypes = length(unique(hapGroup))) %>%
      mutate(win_ed=floor(pos/WinSizeU)) %>%
      group_by(win_ed) %>%
      summarise(MeanMaxNShared=mean(MaxNShared),
                MeanSecondMaxNShared=mean(SecondMaxNShared),
                MeanNHaplotypes=mean(NHaplotypes)) %>%
      mutate(chr=paste0("chr", chrID), 
             pos=win_ed*WinSizeU) %>%
      select(chr, pos, MeanMaxNShared, MeanSecondMaxNShared, MeanNHaplotypes)
  } else{
    merged_table <- rbind(merged_table,
                          table_hapG %>%
                            group_by(pos, hapGroup) %>%
                            summarise(NShared=n()) %>%
                            group_by(pos) %>%
                            summarise(MaxNShared=max(NShared),
                                      SecondMaxNShared=secondHighest(NShared),
                                      NHaplotypes = length(unique(hapGroup))) %>%
                            mutate(win_ed=floor(pos/WinSizeU)) %>%
                            group_by(win_ed) %>%
                            summarise(MeanMaxNShared=mean(MaxNShared),
                                      MeanSecondMaxNShared=mean(SecondMaxNShared),
                                      MeanNHaplotypes=mean(NHaplotypes)) %>%
      mutate(chr=paste0("chr", chrID), pos=win_ed*WinSizeU) %>%
      select(chr, pos, MeanMaxNShared, MeanSecondMaxNShared, MeanNHaplotypes))
  }

  table_hapG %>%
    group_by(pos, hapGroup) %>%
    summarise(NShared=n()) %>%
    group_by(pos) %>%
    summarise(MaxNShared=max(NShared),
              SecondMaxNShared=secondHighest(NShared),
              NHaplotypes = length(unique(hapGroup))) %>%
    mutate(win100=floor(pos/100000)) %>%
    group_by(win100) %>%
    summarise(MeanMaxNShared=mean(MaxNShared),
              MeanSecondMaxNShared=mean(SecondMaxNShared),
              MeanNHaplotypes=mean(NHaplotypes)) %>%
    ggplot(aes(win100/10, MeanMaxNShared))+
    geom_line(alpha=0.5)+
    geom_point(alpha=0.1)+ 
    geom_line(aes(win100/10, MeanSecondMaxNShared), colour="red", alpha=0.5)+
    geom_point(aes(win100/10, MeanSecondMaxNShared), colour="red", alpha=0.1)+
    geom_line(aes(win100/10, MeanNHaplotypes), colour="darkgreen", alpha=0.5)+
    geom_point(aes(win100/10, MeanNHaplotypes), colour="darkgreen", alpha=0.1)+
    scale_x_continuous(breaks=seq(0,100,5), name = "Position (Mb)", expand = c(0,0))+
    ggtitle(paste0("Chr:", chrID))+
    theme_classic()+
    theme(panel.grid.major = element_line(color = "grey80"),
          panel.grid.minor = element_line(color = "grey80"))


  # ggsave(paste0("08_MeanMaxhaplotypeSharing_top2_per100kbWin_chr", chr,".png"), height=5, width=12)

}


merged_table %>%
  filter(MeanNHaplotypes<24) %>%
  ggplot(aes(pos/1000000, MeanNHaplotypes))+
  geom_line(alpha=0.8)+
  # geom_line(aes(win100/10, MeanSecondMaxNShared), colour="red", alpha=0.5)+
  # geom_point(aes(win100/10, MeanSecondMaxNShared), colour="red", alpha=0.1)+
  # geom_line(aes(win100/10, MeanNHaplotypes), colour="darkgreen", alpha=0.5)+
  # geom_point(aes(win100/10, MeanNHaplotypes), colour="darkgreen", alpha=0.1)+
  geom_line(data=centromere_table %>%
              gather(key="posB", value = "pos", -chr) %>%
              mutate(posy = 0), 
            aes(pos/1000000,posy), colour="darkgreen", size=1.5, alpha=0.5) +
  geom_hline(yintercept=mean(merged_table$MeanNHaplotypes), linetype="dashed", color = "red")+
  # scale_y_continuous(breaks=seq(0,1,0.01)) +
  scale_x_continuous(breaks=seq(0,100,15)) +
  # scale_x_continuous(breaks=seq(0,100,5), expand = c(0,0))+
  # ggtitle(paste0("Chr:", chrID))+
  labs(x = "Pos. (Mb)",
       y = "Mean N. Haplo.") +
  facet_grid(. ~ chr, scale="free", space="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12), 
        axis.title = element_text(size=14))


ggsave("10_ed_MeanNHaplot_per500kbWin.png", height=2, width=15)








merged_table %>%
  # filter(MeanNHaplotypes<24) %>%
  ggplot(aes(pos/1000000, MeanMaxNShared))+
  geom_line(alpha=0.8)+
  geom_line(data=centromere_table %>%
              gather(key="posB", value = "pos", -chr) %>%
              mutate(posy = 5), 
            aes(pos/1000000,posy), colour="darkgreen", size=1.5, alpha=0.5) +
  geom_hline(yintercept=mean(merged_table$MeanMaxNShared), linetype="dashed", color = "red")+
  scale_y_continuous(breaks=seq(0,45, 5)) +
  scale_x_continuous(breaks=seq(0,100,15)) +
  labs(x = "Pos. (Mb)",
       y = "Max N. Shared Hap.") +
  facet_grid(. ~ chr, scale="free", space="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12), 
        axis.title = element_text(size=14))


ggsave("11_ed_MeanMaxNShared_per500kbWin.png", height=2, width=15)


