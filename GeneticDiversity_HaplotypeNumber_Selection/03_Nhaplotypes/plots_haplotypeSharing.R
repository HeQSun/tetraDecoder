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

chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

chr <- "06"

for (chr in chr_list) {
table_hapG <- read_table(paste0("haplotypeGroups_min10SS_chr",chr,".txt"), F)
names(table_hapG) <- c("pos", "sample", "hapGroup")

head(table_hapG)

table_hapG_sharedVal <- table_hapG %>%
  group_by(pos, hapGroup) %>%
  summarise(NShared=n()) 


table_hapG %>%
  left_join(table_hapG_sharedVal, by=c("pos", "hapGroup")) %>%
  ggplot(aes(pos/1000000, sample)) +
  geom_tile(aes(fill=NShared))+
  #scico::scale_fill_scico(palette = "lajolla")+
  scale_fill_viridis_c(option = "magma", direction = -1)+
  scale_x_continuous(breaks=seq(0,100,5), name = "Position (Mb)", expand = c(0,0))+
  scale_y_discrete(name = "Haplotype", expand = c(0,0))+
  ggtitle(paste0("Chr:", chr))+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))


ggsave(paste0("06_haplotypeSharing_per10kbWin_chr", chr,".png"), height=7, width=7)


table_hapG %>%
  group_by(pos, hapGroup) %>%
  summarise(NShared=n()) %>%
  group_by(pos) %>%
  summarise(MaxNShared=max(NShared)) %>%
  ggplot(aes(pos/1000000, MaxNShared))+
  geom_line(alpha=0.05)+
  geom_point(alpha=0.1)+
  scale_x_continuous(breaks=seq(0,100,5), name = "Position (Mb)", expand = c(0,0))+
  ggtitle(paste0("Chr:", chr))+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))

ggsave(paste0("07_MaxhaplotypeSharing_per10kbWin_chr", chr,".png"), height=5, width=12)


table_hapG_SecondMax <- table_hapG %>%
  group_by(pos, hapGroup) %>%
  summarise(NShared=n()) %>%
  group_by(pos) %>%
  summarise(SecondMaxNShared=secondHighest(NShared)) 


table_hapG %>%
  group_by(pos, hapGroup) %>%
  summarise(NShared=n()) %>%
  group_by(pos) %>%
  summarise(MaxNShared=max(NShared)) %>%
  left_join(table_hapG_SecondMax, by="pos") %>%
  mutate(win100=floor(pos/100000)) %>%
  group_by(win100) %>%
  summarise(MeanMaxNShared=mean(MaxNShared), 
            MeanSecondMaxNShared=mean(SecondMaxNShared)) %>%
  ggplot(aes(win100/10, MeanMaxNShared))+
  geom_line(alpha=0.5)+
  geom_point(alpha=0.1)+
  geom_line(aes(win100/10, MeanSecondMaxNShared), colour="red", alpha=0.5)+
  geom_point(aes(win100/10, MeanSecondMaxNShared), colour="red", alpha=0.1)+
  scale_x_continuous(breaks=seq(0,100,5), name = "Position (Mb)", expand = c(0,0))+
  ggtitle(paste0("Chr:", chr))+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))
  

ggsave(paste0("08_MeanMaxhaplotypeSharing_top2_per100kbWin_chr", chr,".png"), height=5, width=12)

}

