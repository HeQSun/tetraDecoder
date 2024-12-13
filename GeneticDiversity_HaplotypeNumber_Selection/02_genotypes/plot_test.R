# .libPaths("/dss/dsshome1/lxc03/di36guz2/R")
.libPaths()

setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes")

library(ggplot2)
library(tidyverse)


chr_size.txt

table_chr_size <- read_table("chr_size.txt", F)
names(table_chr_size) <- c("chr", "size")


table_AF <- read_table("AllSamplesPopPar_Win10000_AllCHR_AF.txt", F)
names(table_AF) <- c("chr", "NHaplotypes", "NVar")

head(table_AF)

table_AF %>%
  group_by(NHaplotypes) %>%
  summarise(total_NVar=sum(NVar)) %>%
  ggplot(aes(NHaplotypes, total_NVar))+
  geom_bar(stat="identity")+
  # scale_y_continuous(breaks=seq(0,10000000,500000)) +
  scale_x_continuous(breaks=seq(0,50,5)) +
  labs(x = "Num. Haplotypes",
       y = "N. Variant Sites") +
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor = element_line(color = "grey80"))


# ggsave("01_total_AFDis.png", height=7, width=8)




table_AF %>%
  left_join(table_chr_size, by="chr") %>%
  mutate(pro_NVar=NVar/size) %>%
  ggplot(aes(NHaplotypes, pro_NVar))+
  geom_bar(stat="identity")+
  # scale_y_continuous(breaks=seq(0,10000000,500000)) +
  scale_x_continuous(breaks=seq(0,50,5)) +
  labs(x = "Num. Haplotypes",
       y = "N. Variant Sites/bp") +
  facet_grid(chr ~ .)+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor = element_line(color = "grey80"))

ggsave("02_total_AFDis_perChr.png", height=10, width=8)




# population parameters along the genome:

test_popParameters.txt

table_pop <- read_table("test_popParameters.txt", F)
names(table_pop) <- c("chr", "win", "mean_pi", "mean_pi2", "NVar", "mean_pi_pw", "W_theta", "Taj_D")

head(table_pop)



# Pi:
table_pop %>%
  ggplot(aes(win/1000000, mean_pi_pw))+
  geom_line(alpha=0.2)+
  geom_point(alpha=0.2)+
  scale_y_continuous(breaks=seq(0,1,0.01)) +
  scale_x_continuous(breaks=seq(0,100,5)) +
  labs(x = "Pos. (Mb)",
       y = "Pi") +
  theme_classic()

ggsave("test_Pi.png", height=5, width=12)


# Watterson Theta
table_pop %>%
  ggplot(aes(win/1000000, W_theta))+
  geom_line(alpha=0.2)+
  geom_point(alpha=0.2)+
  scale_y_continuous(breaks=seq(0,1,0.01)) +
  scale_x_continuous(breaks=seq(0,100,5)) +
  labs(x = "Pos. (Mb)",
       y = "W_theta") +
  theme_classic()

ggsave("test_W_theta.png", height=5, width=12)






# Taj_D
table_pop %>% 
  ggplot(aes(win/1000000, Taj_D))+
  geom_line(alpha=0.2)+
  geom_point(alpha=0.2)+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  #scale_y_continuous(breaks=seq(0,1,0.01)) +
  scale_x_continuous(breaks=seq(0,100,5)) +
  labs(x = "Pos. (Mb)",
       y = "W_theta") +
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor = element_line(color = "grey80"))

ggsave("test_Taj_D", height=5, width=12)


