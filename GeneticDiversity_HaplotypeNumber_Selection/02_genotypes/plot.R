# .libPaths("/dss/dsshome1/lxc03/di36guz2/R")
.libPaths()

setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes")

library(ggplot2)
library(tidyverse)


table_chr_size <- read_table("chr_size.txt", F)
names(table_chr_size) <- c("chr", "size")

centromere_table <- read.csv("centromeres.csv", header = T)


table_AF <- read_table("AllSamplesMultiAlleles_Win10000_ConsideringmissingData_AllCHR_AF.txt", F)
names(table_AF) <- c("chr", "fq", "NVar")

head(table_AF)

table_AF %>%
  mutate(interval = cut(fq,
                        seq(0,1,0.1),
                        include.lowest = TRUE,
                        right = FALSE)) %>%
  group_by(interval) %>%
  summarise(total_NVar=sum(NVar)) %>%
  ggplot(aes(interval, total_NVar))+
  geom_bar(stat="identity") +
  # scale_y_continuous(breaks=seq(0,10000000,500000)) +
  # scale_x_continuous(breaks=seq(0,50,5)) +
  labs(x = "Allele Fq.",
       y = "N. Variant Sites") +
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"),
        axis.text.x = element_text(angle = 45, hjust=1))


ggsave("01_total_AFDis.png", height=7, width=8)




table_AF %>%
  mutate(interval = cut(fq,
                        seq(0,1,0.1),
                        include.lowest = TRUE,
                        right = FALSE)) %>%
  group_by(chr, interval) %>%
  summarise(total_NVar=sum(NVar)) %>%
  left_join(table_chr_size, by="chr") %>%
  mutate(pro_NVar=total_NVar/size) %>%
  ggplot(aes(interval, pro_NVar))+
  geom_bar(stat="identity")+
  # scale_y_continuous(breaks=seq(0,10000000,500000)) +
  #scale_x_continuous(breaks=seq(0,50,5)) +
  labs(x = "Allele Fq.",
       y = "N. Variant Sites/bp") +
  facet_grid(chr ~ .)+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"),
        axis.text.x = element_text(angle = 45, hjust=1))

ggsave("02_total_AFDis_perChr.png", height=10, width=8)




# population parameters along the genome:
# chr     win     sum_mean_pi     sum_mean_pi2    SS      meanNsamples    mean_pi W_theta Taj_D

table_pop <- read_table("AllSamplesMultiAlleles_Win10000_ConsideringmissingData_AllCHR_popParameters.txt", F)
names(table_pop) <- c("chr", "win", "sum_pi", "SS", "NVar", "mean_pi_pw", "W_theta", "Taj_D", "Mean_NSamples")

# table_pop$Taj_D <- as.numeric(table_pop$Taj_D)

head(table_pop)
str(table_pop)

# Pi:
table_pop %>%
  filter(mean_pi_pw<0.06) %>%
  ggplot(aes(win/1000000, mean_pi_pw))+
  # geom_line(alpha=0.1)+
  geom_point(alpha=0.04)+
  geom_hline(yintercept=mean(table_pop$mean_pi_pw), linetype="dashed", color = "red")+
  scale_y_continuous(breaks=seq(0,1,0.02)) +
  scale_x_continuous(breaks=seq(0,100,5)) +
  facet_grid(chr ~ .)+
  labs(x = "Pos. (Mb)",
       y = "Pi") +
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"))

# ggsave("03_Pi_AllSamples.png", height=10, width=12)
ggsave("03_Pi_AllSamples_max06.png", height=10, width=12)


# edited version for manuscript:
winSize = 500000
table_pop %>%
  mutate(win_ed=floor(win/winSize)*winSize) %>%
  group_by(chr, win_ed) %>%
  summarise(meanPi=mean(mean_pi_pw, na.rm = T)) %>%
  mutate(low_pi=(meanPi<mean(table_pop$mean_pi_pw))*1) %>%
  filter(meanPi<0.04) %>%
  ggplot(aes(win_ed/1000000, meanPi))+
  # geom_area(aes(x = ifelse(meanPi<mean(table_pop$mean_pi_pw) , win_ed/1000000, win_ed/1000000), y = ifelse(meanPi<mean(table_pop$mean_pi_pw) , meanPi, 0)), fill = "red") +
  # geom_area(aes(x = ifelse(meanPi<mean(table_pop$mean_pi_pw) , win_ed/1000000, win_ed/1000000), y = ifelse(meanPi<0.01 , meanPi, 0)), fill = "red") +
  geom_ribbon(aes(ymin=ifelse(meanPi<mean(table_pop$mean_pi_pw), meanPi,mean(table_pop$mean_pi_pw)), ymax=mean(table_pop$mean_pi_pw)), fill="#dc8f95") +
  # geom_area(aes(fill=factor(low_pi))) +
  geom_line(alpha=0.8)+
  geom_line(data=centromere_table %>%
              gather(key="posB", value = "pos", -chr) %>%
              mutate(posy = 0), 
            aes(pos/1000000,posy), colour="darkgreen", size=1.5, alpha=0.5) +
  #geom_point(alpha=0.04)+
  geom_hline(yintercept=mean(table_pop$mean_pi_pw), linetype="dashed", color = "red")+
  scale_y_continuous(breaks=seq(0,1,0.01)) +
  scale_x_continuous(breaks=seq(0,100,15)) +
  facet_grid(. ~ chr, scale="free", space="free")+
  labs(x = "Pos. (Mb)",
       y = "Pi") +
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12), 
        axis.title = element_text(size=14))

ggsave("03_ed_Pi_AllSamples_mean500kb.png", height=2, width=15)
ggsave("03_ed_Pi_AllSamples_mean500kb.pdf", height=2, width=15)




# test for the association between Pi and centromeres location:

winSize = 500000
CenRegion_table <-  table_pop %>%
  mutate(win_ed=floor(win/winSize)*winSize) %>%
  group_by(chr, win_ed) %>%
  summarise(meanPi=mean(mean_pi_pw)) %>%
  left_join(centromere_table, by="chr") %>%
  mutate(CenRegion=factor(if_else((win_ed>start & win_ed<end), "PC", "No-PC")))


CenRegion_table %>% 
  ggplot(aes(CenRegion, meanPi))+
  geom_boxplot(fill="black", alpha=0.3)+
  facet_grid(. ~ chr)+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        strip.text = element_text(size=12), 
        axis.text.y = element_text(size=12), 
        axis.text.x = element_text(size=12, hjust = 1, angle = 45), 
        axis.title = element_text(size=14))


ggsave("03_2_ed_ComparisonPi_Pericentromeres_mean500kb.pdf", height=3, width=15)


data.frame(chr="chr", 
           W_val="W_val", 
           p.value="p.value", 
           diff_mean="diff_mean") %>%
  write.table(file=paste0("./summary_wilcox.test_pericentromeres_Pi.csv"), sep = ",", row.names = F, col.names = F)

chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
chrID <- "01"

for (chrID in chr_list) {
  chrUsed=paste0("chr",chrID)
  res <- wilcox.test(meanPi ~ CenRegion,
                     paired=F,
                     conf.int=T,
                     data = CenRegion_table %>% filter(chr==chrUsed),
                     exact = FALSE)
  W <- as.vector(res$statistic[1])
  PVal <- res$p.value
  diff_mean <- as.vector(res$estimate[1])
  
  data.frame(chr=chrUsed, 
             W_val=W, 
             p.value=PVal, 
             diff_mean=diff_mean) %>%
    write.table(file=paste0("./summary_wilcox.test_pericentromeres_Pi.csv"), sep = ",", row.names = F, col.names = F, append=T)
}

chrUsed="All"
res <- wilcox.test(meanPi ~ CenRegion,
                   paired=F,
                   conf.int=T,
                   data = CenRegion_table,
                   exact = FALSE)
W <- as.vector(res$statistic[1])
PVal <- res$p.value
diff_mean <- as.vector(res$estimate[1])

data.frame(chr=chrUsed, 
           W_val=W, 
           p.value=PVal, 
           diff_mean=diff_mean) %>%
  write.table(file=paste0("./summary_wilcox.test_pericentromeres_Pi.csv"), sep = ",", row.names = F, col.names = F, append=T)






# Watterson Theta
table_pop %>%
  filter(W_theta<0.06) %>%
  ggplot(aes(win/1000000, W_theta))+
  # geom_line(alpha=0.1)+
  geom_point(alpha=0.04)+
  geom_hline(yintercept=mean(table_pop$W_theta), linetype="dashed", color = "red")+
  scale_y_continuous(breaks=seq(0,1,0.02)) +
  scale_x_continuous(breaks=seq(0,100,5)) +
  labs(x = "Pos. (Mb)",
       y = "W_theta") +
  facet_grid(chr ~ .)+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"))

# ggsave("04_W_theta.png", height=10, width=12)
ggsave("04_W_theta_max06.png", height=10, width=12)






# Taj_D
table_pop %>% 
  ggplot(aes(win/1000000, Taj_D))+
  # geom_line(alpha=0.1)+
  geom_point(alpha=0.04)+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  #scale_y_continuous(breaks=seq(0,1,0.01)) +
  scale_x_continuous(breaks=seq(0,100,5)) +
  labs(x = "Pos. (Mb)",
       y = "Taj_D") +
  facet_grid(chr ~ .)+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))

ggsave("05_Taj_D.png", height=10, width=12)


# eddited version for manuscript:
winSize = 500000
table_pop %>%
  mutate(win_ed=floor(win/winSize)*winSize) %>%
  group_by(chr, win_ed) %>%
  summarise(meanTaj_D=mean(Taj_D, na.rm = T)) %>%
  ggplot(aes(win_ed/1000000, meanTaj_D))+
  geom_ribbon(aes(ymin=ifelse(meanTaj_D<0, meanTaj_D,0), ymax=0), alpha = 0.9, fill="#f0bd9d") +
  geom_line(alpha=0.8)+
  geom_line(data=centromere_table %>%
              gather(key="posB", value = "pos", -chr) %>%
              mutate(posy = (-2)), 
            aes(pos/1000000,posy), colour="darkgreen", size=1.5, alpha=0.5) +
  #geom_point(alpha=0.04)+
  geom_hline(yintercept=mean(table_pop$mean_pi_pw), linetype="dashed", color = "red")+
  # scale_y_continuous(breaks=seq(0,1,0.01)) +
  scale_x_continuous(breaks=seq(0,100,15)) +
  facet_grid(. ~ chr, scale="free", space="free")+
  labs(x = "Pos. (Mb)",
       y = "Taj_D") +
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12), 
        axis.title = element_text(size=14))

ggsave("05_ed_Taj_D_mean500kb.png", height=2, width=15)
ggsave("05_ed_Taj_D_mean500kb.pdf", height=2, width=15)





# Mean Number of Samples:
table_pop %>%
  ggplot(aes(win/1000000, Mean_NSamples))+
  # geom_line(alpha=0.1)+
  geom_point(alpha=0.04)+
  geom_hline(yintercept=mean(table_pop$Mean_NSamples), linetype="dashed", color = "red")+
  geom_line(data=centromere_table %>%
              gather(key="posB", value = "pos", -chr) %>%
              mutate(posy = (-2)), 
            aes(pos/1000000,posy), colour="darkgreen", linewidth=1.5, alpha=0.5) +
  # scale_y_continuous(breaks=seq(0,1,0.02)) +
  scale_x_continuous(breaks=seq(0,100,5)) +
  facet_grid(chr ~ .)+
  labs(x = "Pos. (Mb)",
       y = "Mean NSamples") +
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"))

ggsave("06_MeanNSamples.png", height=10, width=12)
