# .libPaths("/dss/dsshome1/lxc03/di36guz2/R")
.libPaths()

setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes")

library(ggplot2)
library(tidyverse)
library(car)
library(lme4)



# loading Pi table
table_pop_pi <- read_table("AllSamplesMultiAlleles_Win10000_ConsideringmissingData_AllCHR_popParameters.txt", F)
names(table_pop_pi) <- c("chr", "win", "sum_pi", "SS", "NVar", "mean_pi_pw", "W_theta", "Taj_D")

head(table_pop_pi)

winSize = 500000
table_pop_pi_winSize <- table_pop_pi %>%
  mutate(pos=floor(win/winSize)*winSize) %>%
  group_by(chr, pos) %>%
  summarise(meanPi=mean(mean_pi_pw )) 
head(table_pop_pi_winSize)


# loading haplotype Number table
table_pop_NHap <- read_csv("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/03_Nhaplotypes/merged_haplotypeGroups_min10SS_allChr_500KbWin.csv", T)
head(table_pop_NHap)



# loading LD values


LD_table <- read_table(paste0("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/04_LD/ed_within1.5Mb_MeanLDValue_all_maxDis100Mb_WinSize50000.txt"), F)
names(LD_table) <- c("chr", "pos1", "pos2", "r2", "Nvar")

winSize=500000
meanLD_table <- LD_table %>%
  mutate(pos=floor(pos1/winSize)*winSize) %>%
  mutate(pos2w=floor(pos2/winSize)*winSize) %>%
  filter(pos==pos2w) %>%
  group_by(chr, pos) %>%
  dplyr::summarise(mean_r2=mean(r2)) 


head(meanLD_table)



# Merge tables:
table_pop_merged <- table_pop_pi_winSize %>%
  left_join(table_pop_NHap, by=c("chr", "pos")) %>%
  left_join(meanLD_table, by=c("chr", "pos"))




# testing correlations:
# meanPi ~ MeanNHaplotypes

m1_Pi_NHap<-lm(meanPi ~ MeanNHaplotypes, data=table_pop_merged)
Anova(m1_Pi_NHap, Type="III") ; summary(m1_Pi_NHap)

chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

data.frame(chr="chr", 
           r2="r2", 
           adj_r2="adj_r2",
           p.value="p.value") %>%
  write.table(file=paste0("./summary_lm_Pi_MeanNHaplotypes.csv"), sep = ",", row.names = F, col.names = F, append=F)

for (chrID in chr_list) {
  chrUsed=paste0("chr",chrID)
  m1_Pi_NHap<-lm(meanPi ~ MeanNHaplotypes, data=table_pop_merged %>% filter(chr==chrUsed))
  #Anova(m1_Pi_NHap, Type="III") ; 
  res <- summary(m1_Pi_NHap)
  p.value<-res$coefficients[2,4]
  r2 <- res$r.squared
  adj_r2 <- res$adj.r.squared
  data.frame(chr=chrUsed, 
             r2=r2, 
             adj_r2=adj_r2,
             p.value=p.value) %>%
    write.table(file=paste0("./summary_lm_Pi_MeanNHaplotypes.csv"), sep = ",", row.names = F, col.names = F, append=T)
}


m1_Pi_NHap<-lm(meanPi ~ MeanNHaplotypes, data=table_pop_merged)
res <- summary(m1_Pi_NHap)
p.value<-res$coefficients[2,4]
r2 <- res$r.squared
adj_r2 <- res$adj.r.squared
data.frame(chr="All", 
           r2=r2, 
           adj_r2=adj_r2,
           p.value=p.value) %>%
  write.table(file=paste0("./summary_lm_Pi_MeanNHaplotypes.csv"), sep = ",", row.names = F, col.names = F, append=T)



stats_pi_NHap <- read.csv("summary_lm_Pi_MeanNHaplotypes.csv", T) %>%
  mutate(label = paste0("AdjR2=", format(round(adj_r2, 2), nsmall = 2), 
                        "\nP-Val=", format(signif(p.value, digits=3), scientific = T, nsmall = 2))) %>%
  mutate(meanPi=0.03, 
         MeanNHaplotypes=10)

table_pop_merged %>%
  filter(MeanNHaplotypes<20) %>%
  ggplot(aes(MeanNHaplotypes, meanPi))+
  geom_point(alpha=0.4)+
  geom_smooth(method='lm', formula= y~x, alpha=0.5, colour="black")+
  geom_text(data=stats_pi_NHap %>% filter(chr!="All"), 
             aes(label=label), size=3)+
  # facet_grid(chr ~ .)+
  facet_grid(. ~ chr)+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor.x = element_line(color = "grey80"),
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12)) 

ggsave("07_CorNHap_Pi_mean500kb_perChr.pdf", height=4, width=18)


stats_pi_NHap <- read.csv("summary_lm_Pi_MeanNHaplotypes.csv", T) %>%
  mutate(label = paste0("AdjR2=", format(round(adj_r2, 2), nsmall = 2), 
                        "\nP-Val=", format(signif(p.value, digits=3), scientific = T, nsmall = 2))) %>%
  mutate(meanPi=0.03, 
         MeanNHaplotypes=6)

table_pop_merged %>%
  filter(MeanNHaplotypes<20) %>%
  ggplot(aes(MeanNHaplotypes, meanPi))+
  geom_point(alpha=0.4)+
  geom_smooth(method='lm', formula= y~x, alpha=0.5, colour="black")+
  geom_text(data=stats_pi_NHap %>% filter(chr=="All"), 
            aes(label=label), size=3)+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor.x = element_line(color = "grey80"),
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12)) 

ggsave("08_CorNHap_Pi_mean500kb.pdf", height=4, width=5)







# meanPi ~ MeanMaxNShared

m2_Pi_MaxHap <-lm(meanPi ~ MeanMaxNShared, data=table_pop_merged)
Anova(m2_Pi_MaxHap, Type="III") ; summary(m2_Pi_MaxHap)

chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

data.frame(chr="chr", 
           r2="r2", 
           adj_r2="adj_r2",
           p.value="p.value") %>%
  write.table(file=paste0("./summary_lm_Pi_MeanMaxNShared.csv"), sep = ",", row.names = F, col.names = F, append=F)

for (chrID in chr_list) {
  chrUsed=paste0("chr",chrID)
  m2_Pi_MaxHap<-lm(meanPi ~ MeanMaxNShared, data=table_pop_merged %>% filter(chr==chrUsed))
  #Anova(m1_Pi_NHap, Type="III") ; 
  res <- summary(m2_Pi_MaxHap)
  p.value<-res$coefficients[2,4]
  r2 <- res$r.squared
  adj_r2 <- res$adj.r.squared
  data.frame(chr=chrUsed, 
             r2=r2, 
             adj_r2=adj_r2,
             p.value=p.value) %>%
    write.table(file=paste0("./summary_lm_Pi_MeanMaxNShared.csv"), sep = ",", row.names = F, col.names = F, append=T)
}


m2_Pi_MaxHap<-lm(meanPi ~ MeanMaxNShared, data=table_pop_merged)
res <- summary(m2_Pi_MaxHap)
p.value<-res$coefficients[2,4]
r2 <- res$r.squared
adj_r2 <- res$adj.r.squared
data.frame(chr="All", 
           r2=r2, 
           adj_r2=adj_r2,
           p.value=p.value) %>%
  write.table(file=paste0("./summary_lm_Pi_MeanMaxNShared.csv"), sep = ",", row.names = F, col.names = F, append=T)



stats_pi_MeanMaxNShared <- read.csv("./summary_lm_Pi_MeanMaxNShared.csv", T) %>%
  mutate(label = paste0("AdjR2=", format(round(adj_r2, 2), nsmall = 2), 
                        "\nP-Val=", format(signif(p.value, digits=3), scientific = T, nsmall = 2))) %>%
  mutate(meanPi=0.03, 
         MeanMaxNShared=25)

stats_pi_MeanMaxNShared %>% head()

table_pop_merged %>% 
  # filter(MeanNHaplotypes<20) %>%
  ggplot(aes(MeanMaxNShared, meanPi))+
  geom_point(alpha=0.4)+
  geom_smooth(method='lm', formula= y~x, alpha=0.5, colour="black")+
  geom_text(data=stats_pi_MeanMaxNShared %>% filter(chr!="All"), 
            aes(label=label), size=3)+
  # facet_grid(chr ~ .)+
  facet_grid(. ~ chr)+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor.x = element_line(color = "grey80"),
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12)) 

ggsave("09_CorMeanMaxNShared_Pi_mean500kb_perChr.pdf", height=4, width=18)


stats_pi_MeanMaxNShared <- read.csv("./summary_lm_Pi_MeanMaxNShared.csv", T) %>%
  mutate(label = paste0("AdjR2=", format(round(adj_r2, 2), nsmall = 2), 
                        "\nP-Val=", format(signif(p.value, digits=3), scientific = T, nsmall = 2))) %>%
  mutate(meanPi=0.03, 
         MeanMaxNShared=25)

table_pop_merged %>%
  ggplot(aes(MeanMaxNShared, meanPi))+
  geom_point(alpha=0.4)+
  geom_smooth(method='lm', formula= y~x, alpha=0.5, colour="black")+
  geom_text(data=stats_pi_MeanMaxNShared %>% filter(chr=="All"), 
            aes(label=label), size=3)+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor.x = element_line(color = "grey80"),
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12)) 

ggsave("10_CorMeanMaxNShared_Pi_mean500kb.pdf", height=4, width=5)


m3_NHap_MaxHap <-lm(MeanNHaplotypes ~ MeanMaxNShared, data=table_pop_merged)
Anova(m3_NHap_MaxHap, Type="III") ; summary(m3_NHap_MaxHap)







## LD Vs Pi

m4_Pi_LD <-lm(meanPi ~ mean_r2, data=table_pop_merged)
Anova(m4_Pi_LD, Type="III") ; summary(m4_Pi_LD)

chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

data.frame(chr="chr", 
           r2="r2", 
           adj_r2="adj_r2",
           p.value="p.value") %>%
  write.table(file=paste0("./summary_lm_Pi_R2.csv"), sep = ",", row.names = F, col.names = F, append=F)

for (chrID in chr_list) {
  chrUsed=paste0("chr",chrID)
  m4_Pi_LD<-lm(meanPi ~ mean_r2, data=table_pop_merged %>% filter(chr==chrUsed))
  #Anova(m4_Pi_LD, Type="III") ; 
  res <- summary(m4_Pi_LD)
  p.value<-res$coefficients[2,4]
  r2 <- res$r.squared
  adj_r2 <- res$adj.r.squared
  data.frame(chr=chrUsed, 
             r2=r2, 
             adj_r2=adj_r2,
             p.value=p.value) %>%
    write.table(file=paste0("./summary_lm_Pi_R2.csv"), sep = ",", row.names = F, col.names = F, append=T)
}


m4_Pi_LD<-lm(meanPi ~ mean_r2, data=table_pop_merged)
res <- summary(m4_Pi_LD)
p.value<-res$coefficients[2,4]
r2 <- res$r.squared
adj_r2 <- res$adj.r.squared
data.frame(chr="All", 
           r2=r2, 
           adj_r2=adj_r2,
           p.value=p.value) %>%
  write.table(file=paste0("./summary_lm_Pi_R2.csv"), sep = ",", row.names = F, col.names = F, append=T)



stats_pi_R2 <- read.csv("summary_lm_Pi_R2.csv", T) %>%
  mutate(label = paste0("AdjR2=", format(round(adj_r2, 2), nsmall = 2), 
                        "\nP-Val=", format(signif(p.value, digits=3), scientific = T, nsmall = 2))) %>%
  mutate(meanPi=0.03, 
         mean_r2=0.4)



table_pop_merged %>%
  ggplot(aes(mean_r2, meanPi))+
  geom_point(alpha=0.4)+
  geom_smooth(method='lm', formula= y~x, alpha=0.5, colour="black")+
  geom_text(data=stats_pi_R2 %>% filter(chr!="All"), 
            aes(label=label), size=3)+
  ylim(0, NA)+
  facet_grid(. ~ chr)+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor.x = element_line(color = "grey80"),
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12)) 

ggsave("11_CorR2_Pi_mean500kb_perChr.pdf", height=4, width=18)



stats_pi_NHap <- read.csv("summary_lm_Pi_R2.csv", T) %>%
  mutate(label = paste0("AdjR2=", format(round(adj_r2, 2), nsmall = 2), 
                        "\nP-Val=", format(signif(p.value, digits=3), scientific = T, nsmall = 2))) %>%
  mutate(meanPi=0.03, 
         mean_r2=0.4)

table_pop_merged %>%
  ggplot(aes(mean_r2, meanPi))+
  geom_point(alpha=0.4)+
  geom_smooth(method='lm', formula= y~x, alpha=0.5, colour="black")+
  geom_text(data=stats_pi_R2 %>% filter(chr=="All"), 
            aes(label=label), size=3)+
  ylim(0, NA)+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor.x = element_line(color = "grey80"),
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12)) 

ggsave("12_CorR2_Pi_mean500kb.pdf", height=4, width=5)







