# .libPaths("/dss/dsshome1/lxc03/di36guz2/R")
# mamba activate env_others
.libPaths()

setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/03_Nhaplotypes")

library(ggplot2)
library(tidyverse)

chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

chr <- "03"

for (chr in chr_list) {
  table_CountSS <- read_table(paste0("Total_NSS_chr",chr,".txt"), F)
  names(table_CountSS) <- c("NSS", "Count")
  
  table_CountSS %>%
    ggplot(aes(NSS, Count))+
    geom_point()+
    theme_classic()+
    theme(panel.grid.major = element_line(color = "grey80"),
          panel.grid.minor = element_line(color = "grey80"))
  
  ggsave(paste0("01_distribution_pairwiseSS_per10kbWin_perWindow_chr", chr,".png"), height=7, width=8)
  
  table_CountSS %>%
    filter(NSS>2) %>%
    filter(NSS<500) %>%
    ggplot(aes(NSS, Count))+
    geom_point()+
    geom_line() +
    theme_classic()+
    theme(panel.grid.major = element_line(color = "grey80"),
          panel.grid.minor = element_line(color = "grey80"))
  
  ggsave(paste0("02_distribution_pairwiseSS_per10kbWin_perWindow_chr", chr,"_2to500NSS.png"), height=7, width=8)
  
  table_CountSS %>%
    filter(NSS>2) %>%
    filter(NSS<100) %>%
    ggplot(aes(NSS, Count))+
    geom_point()+
    geom_line() +
    theme_classic()+
    theme(panel.grid.major = element_line(color = "grey80"),
          panel.grid.minor = element_line(color = "grey80"))
  
  ggsave(paste0("03_distribution_pairwiseSS_per10kbWin_perWindow_chr", chr,"_2to100NSS.png"), height=7, width=8)
}


table_CountSS <- read_table("all_Total_NSS.txt", F)
names(table_CountSS) <- c("chr", "NSS", "Count")

table_CountSS %>%
  filter(NSS>2) %>%
  filter(NSS<500) %>%
  ggplot(aes(NSS, Count/1000))+
  geom_point()+
  geom_line() +
  labs(x = "NSS",
       y = "Count (10^3)") +
  facet_grid(chr ~ ., scale="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))

ggsave("04_distribution_pairwiseSS_per10kbWin_perWindow_allchr_2to500NSS.png", height=12, width=6)





table_CountSS %>% 
  # filter(NSS>0) %>%
  # filter(NSS<500) %>%
  group_by(NSS) %>%
  summarise(CountTotal=sum(Count)) %>%
  ggplot(aes(NSS, CountTotal/1000))+
  geom_point()+
  geom_line() +
  labs(x = "NSS",
       y = "Count (10^3)") +
  # facet_grid(chr ~ ., scale="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor = element_line(color = "grey80"),
        axis.text = element_text(size=24), 
        axis.title = element_text(size=24))

ggsave("04_ed_distribution_pairwiseSS_per10kbWin_perWindow_all.png", height=5, width=7)


table_CountSS %>% 
  filter(NSS>4) %>%
  filter(NSS<300) %>%
  group_by(NSS) %>%
  summarise(CountTotal=sum(Count)) %>%
  ggplot(aes(NSS, CountTotal/1000))+
  geom_point()+
  geom_line() +
  labs(x = "NSS",
       y = "Count (10^3)") +
  # facet_grid(chr ~ ., scale="free")+
  scale_x_continuous(breaks=c(5, seq(50,500,50))) +
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor = element_line(color = "grey80"),
        axis.text = element_text(size=24), 
        axis.title = element_text(size=24))

ggsave("04_ed_distribution_pairwiseSS_per10kbWin_perWindow_all_5to300NSS.png", height=5, width=7)






table_CountSS %>%
  filter(NSS>2) %>%
  filter(NSS<50) %>%
  ggplot(aes(NSS, Count/1000))+
  geom_point()+
  geom_line() +
  labs(x = "NSS",
       y = "Count (10^3)") +
  facet_grid(chr ~ ., scale="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))

ggsave("05_distribution_pairwiseSS_per10kbWin_perWindow_allchr_2to50NSS.png", height=12, width=6)



