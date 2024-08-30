# Genotype file / Genetic diversity and Taj's D



## Job all Chromosomes



python script: geneticDiversity_genotyping.py

```python
import numpy as np
from collections import defaultdict
import math
import sys
import os

#file="merged_genotypes.chr01.txt"
#list_samples_file = str("list_haplotypesNames.chr01.txt")
#winSize = int(10000)
#output = "test"

file = str(sys.argv[1])
list_samples_file = str(sys.argv[2])
winSize = int(sys.argv[3])
output = str(sys.argv[4])

def round_up_to_nearest_win(num, win_size):
	return math.ceil(int(num) / win_size) * win_size

# load samples names
with open(list_samples_file, "r") as f:
	list_samples = f.read().splitlines()


Nsamples= len(list_samples)
samples = np.array(list_samples)


#output file 1: table with genotypes:
output_genotype = open(f'{output}_genotypeTable.txt', "w")

samples_join = '\t'.join(map(str, samples))
output_genotype.write(f'CHR\tStart\tEnd\tRef\tAlt\t{samples_join}\n')

# output file 2: table qieh allele fq:
output_AF = open(f'{output}_AF.txt', "w")
AF_dis = defaultdict(int)

#output file 3: population parameters per Window
output_popParameters = open(f'{output}_popParameters.txt', "w")

# Calculate parametes:
current_win = 0
mean_pi = 0
mean_pi2 = 0
NVar = 0

# calculate harmonic mean for according to N to calculate Watterson Theta and Tajima's D
total_sum_harmonicMean = 0
total_sum_harmonicMean2 = 0
for sampleN in range(1,Nsamples):
	#print(sampleN)
	total_sum_harmonicMean += 1 / sampleN
	total_sum_harmonicMean2 += 1 / (sampleN**2)

# parameters for Tajima's D
a1_pop1 = float(total_sum_harmonicMean)
a2_pop1 = float(total_sum_harmonicMean2)
b1_pop1 = (Nsamples+1)/(3*(Nsamples-1))
b2_pop1 = (2*((Nsamples**2)+Nsamples+3))/(9*Nsamples*(Nsamples-1))
c1_pop1 = b1_pop1-(1/a1_pop1)
c2_pop1 = b2_pop1-((Nsamples+2)/(a1_pop1*Nsamples))+(a2_pop1/(a1_pop1**2))
e1_pop1 = c1_pop1/a1_pop1
e2_pop1 = c2_pop1/((a1_pop1**2)+a2_pop1)

for line in open(file, "r"):
	chr, start, end, ref, alt, haplotypes  = line.strip().split("\t")
	win = round_up_to_nearest_win(start, winSize)
	# produce genotype table
	index_haplotypes = [np.where(samples == haplotype)[0][0] for haplotype in haplotypes.split(",")]
	genotypes = np.zeros((Nsamples), dtype=int)
	genotypes[index_haplotypes] = 1
	genotypes_join = '\t'.join(map(str, genotypes))
	output_genotype.write(f'{chr}\t{start}\t{end}\t{ref}\t{alt}\t{genotypes_join}\n')
	# calculate allele frequency. values are stored in a dictionary (AF_dis)
	AF_dis[sum(genotypes)] += 1
	# calculate pi (pw) per window:
	p = sum(genotypes)/Nsamples
	q = 1- p
	pi = (Nsamples/(Nsamples-1))*2*p*q
	pi2 = sum(genotypes) * (Nsample-sum(genotypes))
	pi2_unco = pi2/(((Nsamples-1)*Nsamples)/2)
	if win!=current_win:
		if current_win == 0:
			current_win = win
			mean_pi = pi
			mean_pi2 = pi2_unco
			NVar = 1
		else:
			# Watterson Theta:
			W_theta = (NVar/winSize) / total_sum_harmonicMean
			W_theta_unco = NVar / a1_pop1
			# Tajima's D:
			Taj_D = (mean_pi2 - W_theta_unco)/(math.sqrt((e1_pop1*(NVar))+(e2_pop1*(NVar)*(NVar-1))))
			# produce output table
			output_popParameters.write(f'{chr}\t{current_win}\t{mean_pi}\t{mean_pi2}\t{NVar}\t{mean_pi/winSize}\t{W_theta}\t{Taj_D}\n')
			current_win = win
			mean_pi = pi
			mean_pi2 = pi2_unco
			NVar = 1		
	else:
		mean_pi += pi
		mean_pi2 += pi2_unco
		NVar += 1


# produce allele fq table output:
for N, totalVar in AF_dis.items():
	output_AF.write(f'{N}\t{totalVar}\n')

output_genotype.close()
output_AF.close()
output_popParameters.close()


```

Script: geneticDiversity_genotyping.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=12:00:00
#SBATCH -J Pop_Par

list_samples=$1
chr_used=$2
win_size=$3
prefix=$4

#ls -1 ../01_nucmer_syri_vs_DM/nucmer_syri_*_vs_DM/refDM_$chr_used""*/syri.out > list_haplotypes_$chr_used.txt

# produce genotype table per haplotype/sample
for haplotype_file in $(cat $list_samples )
do echo $haplotype_file
haplotype_name=$(echo $haplotype_file | sed 's@.*/01_nucmer_syri_vs_DM/nucmer_syri_@@g' | sed 's@vs_DM/refDM_chr.*_vs_qry_@@g' | sed 's@/syri.out@@g' )
echo $haplotype_name
echo $chr_used
#produce bed file per haplotype with chr, start-1, end, ref, alt, sample
awk -v del="DEL" -v ins="INS" -v snp="SNP" -v haplotype=$haplotype_name '{if($11==del || $11==ins || $11==snp) print $1"\t"$2-1"\t"$3"\t"$4"\t"$5"\t"haplotype}' $haplotype_file | sort | uniq > genotypes.$prefix.$chr_used.$haplotype_name.txt
# produce file with list og haplotype's names:
echo $haplotype_name >> list_haplotypesNames.$prefix.$chr_used.txt
done

# merged file with genotypes with all samples:
cat genotypes.$prefix.$chr_used.* | sort -nk2,2 -k4,4 -k5,5 | awk -v chr=$chr_used 'BEGIN {sample_list=""; startpos=0; endpos=0; ref=""; alt=""} {if($2==startpos && $3==endpos && $4==ref && $5==alt) sample_list=sample_list","$6 ; else {{if(sample_list!="") print chr"\t"startpos"\t"endpos"\t"ref"\t"alt"\t"sample_list}; sample_list=$6; startpos=$2; endpos=$3; ref=$4; alt=$5 }} END {print chr"\t"startpos"\t"endpos"\t"ref"\t"alt"\t"sample_list}' > merged_genotypes.$prefix.$chr_used.txt

# remove files:
rm genotypes.$prefix.$chr_used.*

# calculate population parameters:
python3 geneticDiversity_genotyping.py merged_genotypes.$prefix.$chr_used.txt  list_haplotypesNames.$prefix.$chr_used.txt $win_size $prefix"_"Win$win_size"_"$chr_used

```

Script: geneticDiversity_pythonStep.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=12:00:00
#SBATCH -J Pop_Par

list_samples=$1
chr_used=$2
win_size=$3
prefix=$4

# calculate population parameters:
python3 geneticDiversity_genotyping.py merged_genotypes.$prefix.$chr_used.txt  list_haplotypesNames.$prefix.$chr_used.txt $win_size $prefix"_"Win$win_size"_"$chr_used

```



Job:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes

mamba activate mypython3

# window size 10000
for i in $(seq -w 1 12)
do echo $i
ls -1 /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/01_nucmer_syri_vs_DM/nucmer_syri_*_vs_DM/refDM_chr$i""*/syri.out | grep -v  nucmer_syri_Y_vs_DM > list_haplotypes_chr$i.txt
sbatch geneticDiversity_genotyping.sh list_haplotypes_chr$i.txt chr$i 10000 AllSamplesPopPar
done
# old: 1939011..1939022
# 1949892..1949903

for i in $(seq -w 1 12)
do echo $i
sbatch geneticDiversity_pythonStep.sh list_haplotypes_chr$i.txt chr$i 10000 AllSamplesPopPar
done
# 1949922..1949933


# merge tables:
grep "" AllSamplesPopPar_Win10000_chr*_AF.txt | sed 's@AllSamplesPopPar_Win10000_@@g' | sed 's@_AF.txt:@\t@g' > AllSamplesPopPar_Win10000_AllCHR_AF.txt

cat AllSamplesPopPar_Win10000_chr*_popParameters.txt > AllSamplesPopPar_Win10000_AllCHR_popParameters.txt

grep -v scaffold_ ../../../02_reference/DM_1-3_516_R44_potato_genome_assembly.v6.1/scatter.interval_list | grep "@SQ" | awk '{print $2"\t"$3}' | sed 's@SN:@@g' | sed 's@LN:@@g' > chr_size.txt

# Using other window sizes:
# 50000 and 100000
for i in $(seq -w 1 12)
do echo $i
sbatch geneticDiversity_pythonStep.sh list_haplotypes_chr$i.txt chr$i 50000 AllSamplesPopPar
sbatch geneticDiversity_pythonStep.sh list_haplotypes_chr$i.txt chr$i 100000 AllSamplesPopPar
done
# old: 1939811..1939834

```

Job per Sample:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes

mamba activate mypython3

ls -1 /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/01_nucmer_syri_vs_DM/nucmer_syri_*_vs_DM/refDM_chr01""*/syri.out | sed 's@/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/01_nucmer_syri_vs_DM/nucmer_syri_@@g' | sed 's@_vs_DM/refDM_chr01_vs_qry_ha.*@@g' | sort | uniq > list_Samples_Cultivar.txt

# window size 10000
for sample in $(cat list_Samples_Cultivar.txt)
do echo $sample
for i in $(seq -w 1 12)
do echo $i
ls -1 /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/01_nucmer_syri_vs_DM/nucmer_syri_$sample"_"vs_DM/refDM_chr$i""*/syri.out > list_perSamplehaplotypes_$sample"_"chr$i.txt
sbatch geneticDiversity_genotyping.sh list_perSamplehaplotypes_$sample"_"chr$i.txt chr$i 10000 perSample_PopPar_$sample
done
done
# old: 1939238..1939381

```



### Plots:



```
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

table_pop <- read_table("AllSamplesPopPar_Win10000_AllCHR_popParameters.txt", F)
names(table_pop) <- c("chr", "win", "mean_pi", "mean_pi2", "NVar", "mean_pi_pw", "W_theta", "Taj_D")

head(table_pop)

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



```



## Including missing data:



### Produce genotype table multiple alleles:

python script: genotypeTable_multipleAlleles.py

```{python}
import numpy as np
from collections import defaultdict
import math
import sys
import os

# file="merged_genotypes.AllSamplesMultiAlleles.chr01.txt"
# list_samples_file = str("list_haplotypesNames.AllSamplesMultiAlleles.chr01.txt")
# output = "test"

file = str(sys.argv[1])
list_samples_file = str(sys.argv[2])
output = str(sys.argv[3])


# load samples names
with open(list_samples_file, "r") as f:
	list_samples = f.read().splitlines()


Nsamples= len(list_samples)
samples = np.array(list_samples)


#output file 1: table with genotypes:
output_genotype = open(f'{output}_genotypeTable_miltipleAlleles.txt', "w")

samples_join = '\t'.join(map(str, samples))
output_genotype.write(f'CHR\tStart\tEnd\tRef\tAlt\t{samples_join}\n')

current_start_pos = 0
current_end_pos = 0


for line in open(file, "r"):
	chr, start, end, ref, alt, haplotypes  = line.strip().split("\t")
	#win = round_up_to_nearest_win(start, winSize)
	if (current_start_pos != int(start)) and (current_end_pos != int(end)):
		if current_start_pos != 0:
			genotypes_join = '\t'.join(map(str, genotypes))
			output_genotype.write(f'{chr}\t{current_start_pos}\t{current_end_pos}\t{current_ref}\t{current_alt}\t{genotypes_join}\n')
		# produce genotype table
		index_haplotypes = [np.where(samples == haplotype)[0][0] for haplotype in haplotypes.split(",")]
		genotypes = np.zeros((Nsamples), dtype=int)
		genotypes[index_haplotypes] = 1
		current_haplotype = 2
		current_start_pos = int(start)
		current_end_pos = int(end)
		current_ref = ref
		current_alt = alt
	else:
		index_haplotypes = [np.where(samples == haplotype)[0][0] for haplotype in haplotypes.split(",")]
		genotypes[index_haplotypes] = current_haplotype
		current_haplotype += 1
		current_alt = ','.join(map(str, [current_alt, alt]))


genotypes_join = '\t'.join(map(str, genotypes))
output_genotype.write(f'{chr}\t{current_start_pos}\t{current_end_pos}\t{current_ref}\t{current_alt}\t{genotypes_join}\n')

output_genotype.close()


```



Script: genotypeTable_multipleAlleles.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=12:00:00
#SBATCH -J Genotype_table

list_samples=$1
chr_used=$2
prefix=$3

#ls -1 ../01_nucmer_syri_vs_DM/nucmer_syri_*_vs_DM/refDM_$chr_used""*/syri.out > list_haplotypes_$chr_used.txt

if test -f "list_haplotypesNames.$prefix.$chr_used.txt"; then
    rm list_haplotypesNames.$prefix.$chr_used.txt
fi

# produce genotype table per haplotype/sample
for haplotype_file in $(cat $list_samples )
do echo $haplotype_file
haplotype_name=$(echo $haplotype_file | sed 's@.*/01_nucmer_syri_vs_DM/nucmer_syri_@@g' | sed 's@vs_DM/refDM_chr.*_vs_qry_@@g' | sed 's@/syri.out@@g' )
echo $haplotype_name
echo $chr_used
#produce bed file per haplotype with chr, start-1, end, ref, alt, sample
awk -v del="DEL" -v ins="INS" -v snp="SNP" -v haplotype=$haplotype_name '{if($11==del || $11==ins || $11==snp) print $1"\t"$2-1"\t"$3"\t"$4"\t"$5"\t"haplotype}' $haplotype_file | sort | uniq > genotypes.$prefix.$chr_used.$haplotype_name.txt
# produce file with list og haplotype's names:
echo $haplotype_name >> list_haplotypesNames.$prefix.$chr_used.txt
done

# merged file with genotypes with all samples:
cat genotypes.$prefix.$chr_used.* | sort -nk2,2 -k4,4 -k5,5 | awk -v chr=$chr_used 'BEGIN {sample_list=""; startpos=0; endpos=0; ref=""; alt=""} {if($2==startpos && $3==endpos && $4==ref && $5==alt) sample_list=sample_list","$6 ; else {{if(sample_list!="") print chr"\t"startpos"\t"endpos"\t"ref"\t"alt"\t"sample_list}; sample_list=$6; startpos=$2; endpos=$3; ref=$4; alt=$5 }} END {print chr"\t"startpos"\t"endpos"\t"ref"\t"alt"\t"sample_list}' > merged_genotypes.$prefix.$chr_used.txt

# remove files:
rm genotypes.$prefix.$chr_used.*

# calculate population parameters:
python3 genotypeTable_multipleAlleles.py merged_genotypes.$prefix.$chr_used.txt  list_haplotypesNames.$prefix.$chr_used.txt $prefix"_"$chr_used


```



Job:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes

mamba activate mypython3

# C88 and Ovata were removed (Samples Y and O).
for i in $(seq -w 1 12)
do echo $i
ls -1 /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/01_nucmer_syri_vs_DM/nucmer_syri_*_vs_DM/refDM_chr$i""*/syri.out | grep -v  nucmer_syri_[Y,O]_vs_DM > list_haplotypes_chr$i.txt
sbatch genotypeTable_multipleAlleles.sh list_haplotypes_chr$i.txt chr$i AllSamplesMultiAlleles
done
# 22.11: 2243078..2243089

#for i in $(seq -w 1 12)
#do echo $i
##cp merged_genotypes.AllSamplesPopPar.chr$i.txt merged_genotypes.AllSamplesMultiAlleles.chr$i.txt
#sbatch genotypeTable_multipleAllele_pythonstep.sh list_haplotypes_chr$i.txt chr$i AllSamplesMultiAlleles
#done
# 

for i in {2243002..2243013}; do stop $i ; done

```









### Produce genotype table with missing data:

From the SYRI output we took regions which were identified as "deletions" in the reference, as well as non-aligned regions to add missing data to the genotype table.



python script: geneticDiversity_genotyping_addmissingData.py

```python
import sys
import os
from collections import defaultdict

#python geneticDiversity_genotyping_missingData.py AllSamplesPopPar_Win10000_chr01_genotypeTable.txt test_missing_AllSamplesPopPar_Win10000_chr01_genotypeTable.txt 1 list_test_missing_AllSamplesPopPar_Win10000_chr01_del_chr01_A_hap1.txt A_hap1
# $prefix"_"genotypeTable.txt $index $del_file 

#genotype_table="AllSamplesPopPar_Win10000_chr01_genotypeTable.txt"
#sample_index=int("1")
#del_regions="list_test_missing_AllSamplesPopPar_Win10000_chr01_del_chr01_A_hap1.txt"
#haplotype_name="A_hap1"

genotype_table = str(sys.argv[1])
sample_index = int(sys.argv[2])
del_regions=str(sys.argv[3])
haplotype_name=str(sys.argv[4])

# load position with deletions:

deletions_dic = defaultdict(int)
for line in open(del_regions, "r"):
	list_start, list_end = line.strip().split("\t")
	deletions_dic[list_start] = list_end 

#output file 1: table with genotypes:
output_genotype = open(f'tem_{genotype_table}', "w")

for line in open(genotype_table, "r"):
	if line[0:3]=="CHR":
		output_genotype.write(f'{line}')
	else:
		line_elements  = line.strip().split("\t")
		if line_elements[1] in deletions_dic.keys():
			if line_elements[2] == deletions_dic[line_elements[1]]:
				line_elements[sample_index+4] = "nan"
				line_join = '\t'.join(map(str, line_elements))
				output_genotype.write(f'{line_join}\n')
			else:
				output_genotype.write(f'{line}')
		else:
			output_genotype.write(f'{line}')


output_genotype.close()

```

script: geneticDiversity_genotyping_addmissingData.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00
#SBATCH -J addMissingGT

#python geneticDiversity_genotyping_missingData.sh AllSamplesPopPar_Win10000_chr01_genotypeTable.txt list_haplotypes_chr01.txt AllSamplesPopPar_Win10000_chr01_missing chr01
#genotype_table=AllSamplesPopPar_Win10000_chr01_genotypeTable.txt
#list_samples=list_haplotypes_chr01.txt
#prefix=test_missing_AllSamplesPopPar_Win10000_chr01
#chr=chr01

genotype_table=$1
list_samples=$2
prefix=$3
chr=$4

head -1 $genotype_table | sed 's/\t/\n/'g | grep hap | cat -n | sed 's/ //g' > sample_order.$prefix.txt
cp $genotype_table $prefix"_"genotypeTable.txt
awk '{if($1!="CHR") print $1"\t"$2"\t"$3}' $genotype_table > $prefix"_"genotype.bed 

for haplotype_file in $(cat $list_samples )
do echo $haplotype_file

haplotype_name=$(echo $haplotype_file | sed 's@.*/01_nucmer_syri_vs_DM/nucmer_syri_@@g' | sed 's@vs_DM/refDM_chr.*_vs_qry_@@g' | sed 's@/syri.out@@g' )
echo $haplotype_name
echo $chr
index=$(grep -w $haplotype_name sample_order.$prefix.txt | awk '{print $1}')
echo $index

awk '{if($11=="DEL" && $3-($2-1)>100) print $1"\t"$2-1"\t"$3}' $haplotype_file > $prefix"_"del_$chr"_"$haplotype_name.txt

awk '{if($11=="HDR" && $3-($2-1)>100) print $1"\t"$2-1"\t"$3}' $haplotype_file > $prefix"_"hdr_$chr"_"$haplotype_name.txt

awk '{if($11=="CPL" && $3-($2-1)>100) print $1"\t"$2-1"\t"$3}' $haplotype_file > $prefix"_"cpl_$chr"_"$haplotype_name.txt

awk '{if($1!="-" && $11=="NOTAL" && $3-($2-1)>100) print $1"\t"$2-1"\t"$3}' $haplotype_file > $prefix"_"NOTAL_$chr"_"$haplotype_name.txt

cat $prefix"_"del_$chr"_"$haplotype_name.txt $prefix"_"NOTAL_$chr"_"$haplotype_name.txt $prefix"_"hdr_$chr"_"$haplotype_name.txt $prefix"_"cpl_$chr"_"$haplotype_name.txt  | sort -k1,1 -k2,2n > $prefix"_"ALLMissing_$chr"_"$haplotype_name.txt

bedtools merge -i $prefix"_"ALLMissing_$chr"_"$haplotype_name.txt > $prefix"_"ALLMissing_merged_$chr"_"$haplotype_name.txt

#bedtools intersect  -wa -wb -a $prefix"_"genotype.bed -b $prefix"_"del_$chr"_"$haplotype_name.txt | awk '{varSize=$3-$2 ; delSize=$6-$5 ; if (varSize!=delSize) print $2"\t"$3}' > list_$prefix"_"del_$chr"_"$haplotype_name.txt
#del_file=list_$prefix"_"del_$chr"_"$haplotype_name.txt

#python geneticDiversity_genotyping_addmissingData.py $prefix"_"genotypeTable.txt $index $del_file $haplotype_name

bedtools intersect  -wa -wb -a $prefix"_"genotype.bed -b $prefix"_"ALLMissing_merged_$chr"_"$haplotype_name.txt | awk '{varSize=$3-$2 ; delSize=$6-$5 ; if (varSize!=delSize) print $2"\t"$3}' > list_$prefix"_"ALLMissing_merged_$chr"_"$haplotype_name.txt
ALLMissing_file=list_$prefix"_"ALLMissing_merged_$chr"_"$haplotype_name.txt

python geneticDiversity_genotyping_addmissingData.py $prefix"_"genotypeTable.txt $index $ALLMissing_file $haplotype_name

mv tem_$prefix"_"genotypeTable.txt $prefix"_"genotypeTable.txt

#rm $prefix"_"del_$chr"_"$haplotype_name.txt list_$prefix"_"del_$chr"_"$haplotype_name.txt

rm $prefix"_"del_$chr"_"$haplotype_name.txt $prefix"_"NOTAL_$chr"_"$haplotype_name.txt $prefix"_"ALLMissing_$chr"_"$haplotype_name.txt $prefix"_"ALLMissing_merged_$chr"_"$haplotype_name.txt list_$prefix"_"ALLMissing_merged_$chr"_"$haplotype_name.txt $prefix"_"hdr_$chr"_"$haplotype_name.txt $prefix"_"cpl_$chr"_"$haplotype_name.txt

done

rm sample_order.$prefix.txt $prefix"_"genotype.bed 


```



Job adding missing data to genotype table:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes

mamba activate env_others

for i in $(seq -w 1 12)
do echo $i
genotype_table=AllSamplesMultiAlleles_chr$i""_genotypeTable_miltipleAlleles.txt
list_samples=list_haplotypes_chr$i"".txt
prefix=AllSamplesMultiAlleles_chr$i""_missing
chr=chr$i""
sbatch geneticDiversity_genotyping_addmissingData.sh $genotype_table $list_samples $prefix $chr
done
# 05.04.2024:2452814..2452825

#for i in {2452770..2452781}; do stop $i ; done


```





### Calculate pop parameters with missing data:



Script: geneticDiversityCal_addmissingData.py

```{python}
import numpy as np
from collections import defaultdict
from collections import Counter
import math
import sys
import os

#file="AllSamplesMultiAlleles_chr01_missing_genotypeTable.txt"
#chr_used="chr01"
#Nsamples=int(48)
#winSize=int(10000)
#output="test"

file = str(sys.argv[1])
chr_used = str(sys.argv[2])
Nsamples=int(sys.argv[3])
winSize=int(sys.argv[4])
output=str(sys.argv[5])

def round_up_to_nearest_win(num, win_size):
	return math.ceil(int(num) / win_size) * win_size


# output file 1: table with allele fq:
output_AF = open(f'{output}_AF.txt', "w")
AF_dis = defaultdict(int)


#output file 2: population parameters per Window
output_popParameters = open(f'{output}_popParameters.txt', "w")


# calculate harmonic mean according to N to calculate Watterson Theta and Tajima's D
total_sum_harmonicMean = 0
total_sum_harmonicMean2 = 0
for sampleN in range(1,Nsamples):
	#print(sampleN)
	total_sum_harmonicMean += 1 / sampleN
	total_sum_harmonicMean2 += 1 / (sampleN**2)

# parameters for Tajima's D
a1_pop1 = float(total_sum_harmonicMean)
a2_pop1 = float(total_sum_harmonicMean2)
b1_pop1 = (Nsamples+1)/(3*(Nsamples-1))
b2_pop1 = (2*((Nsamples**2)+Nsamples+3))/(9*Nsamples*(Nsamples-1))
c1_pop1 = b1_pop1-(1/a1_pop1)
c2_pop1 = b2_pop1-((Nsamples+2)/(a1_pop1*Nsamples))+(a2_pop1/(a1_pop1**2))
total_e1_pop1 = c1_pop1/a1_pop1
total_e2_pop1 = c2_pop1/((a1_pop1**2)+a2_pop1)


current_win = 0

for line in open(file, "r"):
	#chr, start, end, ref, alt, haplotypes  = line.strip().split("\t")
	if line[0:3] == "CHR":
		output_popParameters.write(f'chr\twin\tsum_pi\tSS\tNvar\tmean_pi\tW_theta\tTaj_D\tmean_Nsamples\n')
	else:
		line_values  = line.strip().split("\t")
		chr, start, end, ref, alt = line_values[0:5]
		win = round_up_to_nearest_win(start, winSize)
		#genotypes = np.array(list(map(int, line_values[5:])))
		genotypes = [int(x) for x in line_values[5:] if x!='nan']
		Nsamples = len(genotypes)
		# calculate Pi per base:
		if Nsamples > 1:
			total_sum_harmonicMean = 0
			for sampleN in range(1,Nsamples):
				total_sum_harmonicMean += 1 / sampleN
			# parameters for Tajima's D
			a1_pop1 = float(total_sum_harmonicMean)
		if len(set(genotypes)) > 1:
			allele_count = list(Counter(genotypes).values())
			NumSS = 0
			for i in range(0,len(allele_count)):
				for j in range(0,len(allele_count)):
					if j > i:
						NumSS += allele_count[i]*allele_count[j]
			NCop = (Nsamples*(Nsamples-1))/2
			pi = NumSS/NCop
			SS = 1
			w_theta_local = 1 / a1_pop1
			AF = round(((max(allele_count)/(Nsamples))*100))
			if AF < 50:
				AF_dis[AF/100] += 1
			else:
				AF_dis[(100-AF)/100] += 1
		else:
			pi = 0
			SS = 0
			w_theta_local = 0
			AF_dis[0] += 1
		# add value per window:
		if win!=current_win:
			if current_win == 0:
				current_win = win
				sum_pi = pi
				total_SS = SS
				W_theta = w_theta_local
				NVar = 1
				mean_Nsamples = Nsamples
			else:
				# Watterson Theta:
				total_W_theta = W_theta/winSize
				# Tajima's D:
				if total_SS > 0:
					Taj_D = (sum_pi - W_theta)/(math.sqrt((total_e1_pop1*(total_SS))+(total_e2_pop1*(total_SS)*(total_SS-1))))
				else:
					Taj_D = 'nan'
				# produce output table
				output_popParameters.write(f'{chr}\t{current_win}\t{sum_pi}\t{total_SS}\t{NVar}\t{sum_pi/winSize}\t{total_W_theta}\t{Taj_D}\t{mean_Nsamples/NVar}\n')
				current_win = win
				sum_pi = pi
				total_SS = SS
				W_theta = w_theta_local
				NVar = 1
				mean_Nsamples = Nsamples
		else:
			sum_pi += pi
			total_SS += SS
			W_theta += w_theta_local
			NVar += 1
			mean_Nsamples += Nsamples



# produce allele fq table output:
for fq, totalVar in AF_dis.items():
	output_AF.write(f'{fq}\t{totalVar}\n')

output_AF.close()
output_popParameters.close()


```



Bash script: geneticDiversityCal_addmissingData.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=12:00:00
#SBATCH -J Pop_Par

genotype_table=$1
chr_used=$2
Nsamples=$3
win_size=$4
prefix=$5

#genotype_table="AllSamplesMultiAlleles_chr12_missing_genotypeTable.txt"
#chr_used="chr01"
#Nsamples=int(48)
#win_size=int(10000)
#prefix_output="test"

# calculate population parameters:
python3 geneticDiversityCal_addmissingData.py $genotype_table $chr_used $Nsamples $win_size $prefix_output



```



job:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes

mamba activate env_others


# window size 10000
for i in $(seq -w 1 12)
do echo $i
Nsamples=$(grep CHR AllSamplesMultiAlleles_chr$i""_missing_genotypeTable.txt | awk -F"Alt\t" '{print $2}' | sed 's/\t/\n/g' | wc -l )
sbatch geneticDiversityCal_addmissingData.sh AllSamplesMultiAlleles_chr$i""_missing_genotypeTable.txt chr$i $Nsamples 10000 AllSamplesMultiAlleles_Win10000_ConsideringmissingData_chr$i
done
# 05.04.2024: 2452905..2452916


# window size 200000
for i in $(seq -w 1 12)
do echo $i
Nsamples=$(grep CHR AllSamplesMultiAlleles_chr$i""_missing_genotypeTable.txt | awk -F"Alt\t" '{print $2}' | sed 's/\t/\n/g' | wc -l )
sbatch geneticDiversityCal_addmissingData.sh AllSamplesMultiAlleles_chr$i""_missing_genotypeTable.txt chr$i $Nsamples 200000 AllSamplesMultiAlleles_Win200000_ConsideringmissingData_chr$i
done
# old 22.11:2243222..2243233


```

Merge Tables:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes

# 10kb Window size:
# merge tables:
grep "" AllSamplesMultiAlleles_Win10000_ConsideringmissingData_chr*_AF.txt | sed 's@AllSamplesMultiAlleles_Win10000_ConsideringmissingData_@@g' | sed 's@_AF.txt:@\t@g' > AllSamplesMultiAlleles_Win10000_ConsideringmissingData_AllCHR_AF.txt

cat AllSamplesMultiAlleles_Win10000_ConsideringmissingData_chr*_popParameters.txt | grep -vw win > AllSamplesMultiAlleles_Win10000_ConsideringmissingData_AllCHR_popParameters.txt

grep -v scaffold_ ../../../02_reference/DM_1-3_516_R44_potato_genome_assembly.v6.1/scatter.interval_list | grep "@SQ" | awk '{print $2"\t"$3}' | sed 's@SN:@@g' | sed 's@LN:@@g' > chr_size.txt

sed -i 's/nan/NA/g' AllSamplesMultiAlleles_Win10000_ConsideringmissingData_AllCHR_popParameters.txt


# 200kb Window size:
# merge tables:
grep "" AllSamplesMultiAlleles_Win200000_ConsideringmissingData_chr*_AF.txt | sed 's@AllSamplesMultiAlleles_Win200000_ConsideringmissingData_@@g' | sed 's@_AF.txt:@\t@g' > AllSamplesMultiAlleles_Win200000_ConsideringmissingData_AllCHR_AF.txt

cat AllSamplesMultiAlleles_Win200000_ConsideringmissingData_chr*_popParameters.txt | grep -vw win > AllSamplesMultiAlleles_Win200000_ConsideringmissingData_AllCHR_popParameters.txt

sed -i 's/nan/NA/g' AllSamplesMultiAlleles_Win200000_ConsideringmissingData_AllCHR_popParameters.txt

```



### Plots 10 Kb Windows:



```R
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


# ggsave("01_total_AFDis.png", height=7, width=8)




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

# ggsave("02_total_AFDis_perChr.png", height=10, width=8)




# population parameters along the genome:
# chr     win     sum_mean_pi     sum_mean_pi2    SS      meanNsamples    mean_pi W_theta Taj_D

table_pop <- read_table("AllSamplesMultiAlleles_Win10000_ConsideringmissingData_AllCHR_popParameters.txt", F)
names(table_pop) <- c("chr", "win", "sum_pi", "SS", "NVar", "mean_pi_pw", "W_theta", "Taj_D", "Mean_NSamples")

# table_pop$Taj_D <- as.numeric(table_pop$Taj_D)

head(table_pop)
str(table_pop)

# Pi:
table_pop %>%
  # filter(mean_pi_pw<0.06) %>%
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
ggsave("03_Pi_AllSamples_max06.pdf", height=10, width=12)




table_pop %>% 
  # filter(mean_pi_pw<0.06) %>%
  ggplot(aes(mean_pi_pw))+
  geom_histogram(colour="#dc8f95", fill="#dc8f95", binwidth = 0.002)+
  geom_vline(xintercept = mean(table_pop$mean_pi_pw), color = "red") +
  scale_x_continuous(breaks=seq(0,1,0.02)) +
  scale_y_continuous(breaks=seq(0,10000,1000)) +
  labs(x = "Pi per 10kb Window",
       y = "Window Count") +
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.title = element_text(size=16), 
        axis.text = element_text(size=14))
  
ggsave("03_v2_Pi_AllSamples_max06.png", height=2.6, width=3.9)
ggsave("03_v2_Pi_AllSamples_max06.pdf", height=2.6, width=3.9)






# edited version for manuscript:

winSize = 100000
quantile_5tem <- table_pop %>%
  mutate(win_ed=floor(win/winSize)*winSize) %>%
  group_by(chr, win_ed) %>%
  summarise(meanPi=mean(mean_pi_pw, na.rm = T)) %>%
  mutate(low_pi=(meanPi<mean(table_pop$mean_pi_pw))*1) %>%
  filter(meanPi<0.04) %>%
  pull(meanPi) %>%
  quantile(0.05)
table_pop %>%
  mutate(win_ed=floor(win/winSize)*winSize) %>%
  group_by(chr, win_ed) %>%
  summarise(meanPi=mean(mean_pi_pw, na.rm = T)) %>%
  mutate(low_pi=(meanPi<mean(table_pop$mean_pi_pw))*1) %>%
  filter(meanPi<quantile_5tem) %>%
  dim()

quantile_5tem
366*100000

winSize = 500000

quantile_5 <- table_pop %>%
  mutate(win_ed=floor(win/winSize)*winSize) %>%
  group_by(chr, win_ed) %>%
  summarise(meanPi=mean(mean_pi_pw, na.rm = T)) %>%
  mutate(low_pi=(meanPi<mean(table_pop$mean_pi_pw))*1) %>%
  filter(meanPi<0.04) %>%
  pull(meanPi) %>%
  quantile(0.05)


table_pop %>%
  mutate(win_ed=floor(win/winSize)*winSize) %>%
  group_by(chr, win_ed) %>%
  summarise(meanPi=mean(mean_pi_pw, na.rm = T)) %>%
  mutate(low_pi=(meanPi<mean(table_pop$mean_pi_pw))*1) %>%
  filter(meanPi<0.35) %>%
  ggplot(aes(win_ed/1000000, meanPi))+
  # geom_area(aes(x = ifelse(meanPi<mean(table_pop$mean_pi_pw) , win_ed/1000000, win_ed/1000000), y = ifelse(meanPi<mean(table_pop$mean_pi_pw) , meanPi, 0)), fill = "red") +
  # geom_area(aes(x = ifelse(meanPi<mean(table_pop$mean_pi_pw) , win_ed/1000000, win_ed/1000000), y = ifelse(meanPi<0.01 , meanPi, 0)), fill = "red") +
  geom_ribbon(aes(ymin=ifelse(meanPi<mean(table_pop$mean_pi_pw), meanPi,mean(table_pop$mean_pi_pw)), ymax=mean(table_pop$mean_pi_pw)), fill="#dc8f95") +
  # geom_ribbon(aes(ymin=ifelse(meanPi<quantile_5, meanPi,quantile_5), ymax=quantile_5), fill="#b01e29") +
  geom_point(data=table_pop %>%
               mutate(win_ed=floor(win/winSize)*winSize) %>%
               group_by(chr, win_ed) %>%
               summarise(meanPi=mean(mean_pi_pw, na.rm = T)) %>%
               mutate(low_pi=(meanPi<mean(table_pop$mean_pi_pw))*1) %>%
               filter(meanPi<quantile_5) %>%
               mutate(posy=0.002), aes(win_ed/1000000, posy), colour="#b01e29") +
  # geom_area(aes(fill=factor(low_pi))) +
  geom_line(alpha=0.8)+
  geom_line(data=centromere_table %>%
              gather(key="posB", value = "pos", -chr) %>%
              mutate(posy = 0), 
            aes(pos/1000000,posy), colour="darkgreen", linewidth=2, alpha=0.5) +
  #geom_point(alpha=0.04)+
  geom_hline(yintercept=mean(table_pop$mean_pi_pw), linetype="dashed", color = "red")+
  scale_y_continuous(breaks=seq(0,1,0.01), limits = c(0,0.035)) +
  scale_x_continuous(breaks=seq(0,100,15)) +
  facet_grid(. ~ chr, scale="free", space="free")+
  labs(x = "Pos. (Mb)",
       y = "Pi") +
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12), 
        axis.title = element_text(size=14))

ggsave("03_ed_Pi_AllSamples_mean500kb.png", height=2.5, width=15)
ggsave("03_ed_Pi_AllSamples_mean500kb.pdf", height=2.5, width=15)



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
       y = "Watterson Theta") +
  facet_grid(chr ~ .)+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"))

# ggsave("04_W_theta.png", height=10, width=12)
ggsave("04_W_theta_max06.png", height=10, width=12)
ggsave("04_W_theta_max06.pdf", height=10, width=12)






# Taj_D
table_pop %>% 
  ggplot(aes(win/1000000, Taj_D))+
  # geom_line(alpha=0.1)+
  geom_point(alpha=0.04)+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  #scale_y_continuous(breaks=seq(0,1,0.01)) +
  scale_x_continuous(breaks=seq(0,100,5)) +
  labs(x = "Pos. (Mb)",
       y = "Tajima's D") +
  facet_grid(chr ~ .)+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))

ggsave("05_Taj_D.png", height=10, width=12)
ggsave("05_Taj_D.pdf", height=10, width=12)


# eddited version for manuscript:

winSize = 500000

quantile_5 <- table_pop %>%
  mutate(win_ed=floor(win/winSize)*winSize) %>%
  group_by(chr, win_ed) %>%
  summarise(meanTaj_D=mean(Taj_D, na.rm = T)) %>%
  pull(meanTaj_D) %>%
  quantile(0.04)

quantile_5

table_pop %>%
  mutate(win_ed=floor(win/winSize)*winSize) %>%
  group_by(chr, win_ed) %>%
  summarise(meanTaj_D=mean(Taj_D, na.rm = T)) %>%
  ggplot(aes(win_ed/1000000, meanTaj_D))+
  geom_ribbon(aes(ymin=ifelse(meanTaj_D<0, meanTaj_D,0), ymax=0), alpha = 0.9, fill="#f0bd9d") +
  geom_point(data=table_pop %>%
               mutate(win_ed=floor(win/winSize)*winSize) %>%
               group_by(chr, win_ed) %>%
               summarise(meanTaj_D=mean(Taj_D, na.rm = T)) %>%
               filter(meanTaj_D<quantile_5) %>%
               mutate(posy=(-2.2)), aes(win_ed/1000000, posy), colour="#d86e2c") +
  geom_line(alpha=0.8)+
  geom_line(data=centromere_table %>%
              gather(key="posB", value = "pos", -chr) %>%
              mutate(posy = (-2.4)), 
            aes(pos/1000000,posy), colour="darkgreen", linewidth=1.5, alpha=0.5) +
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

ggsave("05_ed_Taj_D_mean500kb.png", height=2.5, width=15)
ggsave("05_ed_Taj_D_mean500kb.pdf", height=2.5, width=15)






winSize = 100000
table_pop %>%
  mutate(win_ed=floor(win/winSize)*winSize) %>%
  group_by(chr, win_ed) %>%
  summarise(meanTaj_D=mean(Taj_D, na.rm = T)) %>%
  filter(meanTaj_D<quantile_5) %>%
  filter(chr=="chr12") %>%
  # filter(Taj_D<quantile_5) %>%
  data.frame()




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
       y = "Mean N. Samples") +
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"))

ggsave("06_MeanNSamples.png", height=10, width=12)
ggsave("06_MeanNSamples.pdf", height=10, width=12)




# eddited version for manuscript:
winSize = 500000
table_pop %>%
  mutate(win_ed=floor(win/winSize)*winSize) %>%
  group_by(chr, win_ed) %>%
  summarise(MeanNSamples=mean(Mean_NSamples, na.rm = T)) %>%
  ggplot(aes(win_ed/1000000, MeanNSamples))+
  # geom_ribbon(aes(ymin=ifelse(meanTaj_D<0, meanTaj_D,0), ymax=0), alpha = 0.9, fill="#f0bd9d") +
  geom_line(alpha=0.8)+
  geom_line(data=centromere_table %>%
              gather(key="posB", value = "pos", -chr) %>%
              mutate(posy = (10)), 
            aes(pos/1000000,posy), colour="darkgreen", size=1.5, alpha=0.5) +
  #geom_point(alpha=0.04)+
  geom_hline(yintercept=mean(table_pop$Mean_NSamples), linetype="dashed", color = "red")+
  # scale_y_continuous(breaks=seq(0,1,0.01)) +
  scale_x_continuous(breaks=seq(0,100,15)) +
  facet_grid(. ~ chr, scale="free", space="free")+
  labs(x = "Pos. (Mb)",
       y = "Mean N. Samples") +
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12), 
        axis.title = element_text(size=14))

ggsave("06_ed_MeanNSamples_mean500kb.png", height=2, width=15)
ggsave("06_ed_MeanNSamples_mean500kb.pdf", height=2, width=15)



```



### Diversity within individuals (Pi):

Script: PiperInd.py

```python
import numpy as np
from collections import defaultdict
from collections import Counter
import math
import sys
import os

#file="/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr01_missing_genotypeTable.txt"
#chr_used="chr01"
#winSize=int(10000)
#output="test"
#haplotypes =str("1,2,3,4")

file = str(sys.argv[1])
chr_used = str(sys.argv[2])
winSize=int(sys.argv[3])
output=str(sys.argv[4])
haplotypes = str(sys.argv[5])



haplotype_list= [int(x) for x in haplotypes.split(',')]

def round_down_to_nearest_win(num, win_size):
	return math.floor(int(num) / win_size) * win_size

#output file 1: population parameters per Window for a single individual using the index of the haplotypes
output_popParameters = open(f'{output}_Pi.txt', "w")

current_win = -1

for line in open(file, "r"):
	#chr, start, end, ref, alt, haplotypes  = line.strip().split("\t")
	if line[0:3] == "CHR":
		output_popParameters.write(f'chr\twin\tsum_pi\tSS\tNvar\tmean_pi\tmean_Nsamples\n')
	else:
		line_values  = line.strip().split("\t")
		chr, start, end, ref, alt = line_values[0:5]
		win = round_down_to_nearest_win(start, winSize)
		#genotypes = np.array(list(map(int, line_values[5:])))
		# extract genotypes only for the haplotypes given. (From one individual)
		genotypes =[int(line_values[i+4]) for i in haplotype_list if line_values[i+4] != 'nan']
		Nsamples = len(genotypes)
		# calculate Pi per base:
		if len(set(genotypes)) > 1:
			allele_count = list(Counter(genotypes).values())
			NumSS = 0
			for i in range(0,len(allele_count)):
				for j in range(0,len(allele_count)):
					if j > i:
						NumSS += allele_count[i]*allele_count[j]
			NCop = (Nsamples*(Nsamples-1))/2
			pi = NumSS/NCop
			SS = 1
		else:
			pi = 0
			SS = 0
		# add value per window:
		if win!=current_win:
			if current_win == -1:
				current_win = win
				sum_pi = pi
				total_SS = SS
				NVar = 1
				mean_Nsamples = Nsamples
			else:
				# produce output table
				output_popParameters.write(f'{chr}\t{current_win}\t{sum_pi}\t{total_SS}\t{NVar}\t{sum_pi/winSize}\t{mean_Nsamples/NVar}\n')
				current_win = win
				sum_pi = pi
				total_SS = SS
				NVar = 1
				mean_Nsamples = Nsamples
		else:
			sum_pi += pi
			total_SS += SS
			NVar += 1
			mean_Nsamples += Nsamples



output_popParameters.close()


```



Bash script: PiperInd.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=12:00:00
#SBATCH -J Pop_Par

genotype_table=$1
chr_used=$2
win_size=$3
prefix=$4
haplotypes=$5

#file="/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr01_missing_genotypeTable.txt"
#chr_used="chr01"
#winSize=int(10000)
#output="test"
#haplotypes =str("1,2,3,4")

# calculate population parameters:
python3 PiperInd.py $genotype_table $chr_used $win_size $prefix $haplotypes


```



job:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/05_diversity_withinInd

mamba activate mypython3

# window size 10000
for i in $(seq 1 4 40)
do echo $i
sample=$(grep CHR /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr01_missing_genotypeTable.txt | awk -F"Alt\t" '{print $2}' | sed 's/\t/\n/g' | sed -n "$i"p | sed 's/_.*//g')
for chr in $(seq -w 1 12)
do echo $chr 
hap2=$(echo $i | awk '{print $1+1}')
hap3=$(echo $i | awk '{print $1+2}')
hap4=$(echo $i | awk '{print $1+3}')
# run job:
sbatch PiperInd.sh /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr$chr""_missing_genotypeTable.txt chr$chr 10000 TablePi_Win10000_ConsideringmissingData_chr$chr"_Sample"$sample $i","$hap2","$hap3","$hap4
done
done
# 05.04.2024: 2452919..2453038

```



Randomising haplotypes:

```bash

cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/05_diversity_withinInd/RandomSampling


mamba activate mypython3

for sample in $(seq -w 1 50)
do echo $sample
for chr in $(seq -w 1 12)
do echo $chr 
hap1=$(shuf -i 1-40 -n 1 )
hap2=$(shuf -i 1-40 -n 1 )
hap3=$(shuf -i 1-40 -n 1 )
hap4=$(shuf -i 1-40 -n 1 )
# run job:
sbatch PiperInd.sh /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr$chr""_missing_genotypeTable.txt chr$chr 10000 TablePi_Win10000_ConsideringmissingData_chr$chr"_Sample"$sample $hap1","$hap2","$hap3","$hap4
done
done
# 05.04.2024: 2453939..2454538


cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/05_diversity_withinInd/

grep -v "mean_Nsamples" RandomSampling/TablePi_Win10000_ConsideringmissingData_chr*txt | sed 's@RandomSampling/TablePi_Win10000_ConsideringmissingData_chr.._Sample@@g' | sed 's@_Pi.txt:@\t@g' > RandomSampling_TablePi_Win10000_ConsideringmissingData_AllChr_AllSamples.txt

grep -v "mean_Nsamples" TablePi_Win10000_ConsideringmissingData_chr*txt | sed 's@TablePi_Win10000_ConsideringmissingData_chr.._Sample@@g' | sed 's@_Pi.txt:@\t@g' > Cultivars_TablePi_Win10000_ConsideringmissingData_AllChr_AllSamples.txt
# done


# for i in {2215213..2215555}; do stop $i ; rm slurm-$i* ; done

```



#### Plots:



```R
# .libPaths("/dss/dsshome1/lxc03/di36guz2/R")
# mamba activate env_others
.libPaths()

setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/05_diversity_withinInd")

library(ggplot2)
library(tidyverse)

centromere_table <- read.csv("centromeres.csv", header = T)

chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

chrID <- "01"

# for (chrID in chr_list) {

table_RandomSampling <- read_table("RandomSampling_TablePi_Win10000_ConsideringmissingData_AllChr_AllSamples.txt", F)
names(table_RandomSampling) <- c("sample", "chr", "pos", "sum_pi", "SS", "NVar", "meanPi", "meanNSamples")

table_cultivars <- read_table("Cultivars_TablePi_Win10000_ConsideringmissingData_AllChr_AllSamples.txt", F)
names(table_cultivars) <- c("sample", "chr", "pos", "sum_pi", "SS", "NVar", "meanPi", "meanNSamples")




table_RandomSampling %>%
  group_by(sample, chr) %>%
  summarise(totalMeanPi=mean(meanPi)) %>%
  ggplot(aes(totalMeanPi))+
  geom_histogram(fill="black") +
  geom_vline(data=table_cultivars %>%
               group_by(sample, chr) %>%
               summarise(totalMeanPi=mean(meanPi)) , 
             aes(xintercept = totalMeanPi), color = "red", alpha=0.5) +
  labs(x = "Mean Pi/bp",
       y = "Count") +
  facet_grid(. ~ chr)+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))


ggsave(paste0("01_MeanPi_perSample_perChr.png"), height=2, width=12)
ggsave(paste0("01_MeanPi_perSample_perChr.pdf"), height=2, width=12)
  



winSize = 500000
table_RandomSampling %>%
  group_by(chr, pos) %>%
  summarise(maxPi=max(meanPi), 
            totalMeanPi=mean(meanPi, na.rm = T), 
            minPi=min(meanPi)) %>%
  mutate(win_ed=floor(pos/winSize)*winSize) %>%
  group_by(chr, win_ed) %>%
  summarise(MeanmaxPi=mean(maxPi, na.rm = T), 
            MeantotalMeanPi=mean(totalMeanPi, na.rm = T), 
            MeanminPi=mean(minPi, na.rm = T)) %>%
  filter(MeanmaxPi<0.06)%>%
  ggplot(aes(win_ed/1000000, MeantotalMeanPi))+
  geom_line(aes(win_ed/1000000, MeanminPi), alpha=0.4)+
  geom_line(aes(win_ed/1000000, MeanmaxPi), alpha=0.5)+
  geom_line(data=table_cultivars %>%
              mutate(win_ed=floor(pos/winSize)*winSize) %>%
              group_by(sample, chr, win_ed) %>%
              summarise(MeantotalMeanPi=mean(meanPi, na.rm = T)), 
            aes(group=sample), colour="red", alpha=0.4)+
  geom_line(alpha=0.4, colour="blue")+
  geom_line(data=centromere_table %>%
              gather(key="posB", value = "pos", -chr) %>%
              mutate(posy = 0), 
            aes(pos/1000000,posy), colour="darkgreen", size=1.5, alpha=0.5) +
  scale_x_continuous(breaks=seq(0,100,15)) +
  labs(x = "Pos. (Mb)",
       y = "Pi") +
  facet_grid(. ~ chr, scale="free", space="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12), 
        axis.title = element_text(size=14))


ggsave("02_meanPi_alongGenome_mean500kb.png", height=2, width=17)
ggsave("02_meanPi_alongGenome_mean500kb.pdf", height=2, width=17)



```



### Diversity within individuals (Theta):



Script: ThetaperInd.py

```python
import numpy as np
from collections import defaultdict
from collections import Counter
import math
import sys
import os


#file="/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr01_missing_genotypeTable.txt"
#chr_used="chr01"
#winSize=int(10000)
#output="test"
#haplotypes =str("1,2,3,4")

file = str(sys.argv[1])
chr_used = str(sys.argv[2])
winSize=int(sys.argv[3])
output=str(sys.argv[4])
haplotypes = str(sys.argv[5])


haplotype_list= [int(x) for x in haplotypes.split(',')]

def round_up_to_nearest_win(num, win_size):
  return math.floor(int(num) / win_size) * win_size


#output file 1: population parameters per Window
output_popParameters = open(f'{output}.Theta', "w")

current_win = -1

for line in open(file, "r"):
  #chr, start, end, ref, alt, haplotypes  = line.strip().split("\t")
  if line[0:3] == "CHR":
    output_popParameters.write(f'chr\twin\tsum_theta\tSS\tNvar\tW_theta\tmean_Nsamples\n')
  else:
    line_values  = line.strip().split("\t")
    chr, start, end, ref, alt = line_values[0:5]
    win = round_up_to_nearest_win(start, winSize)
    genotypes =[int(line_values[i+4]) for i in haplotype_list if line_values[i+4] != 'nan']
    Nsamples = len(genotypes)
    # calculate theta per base:
    if len(set(genotypes)) > 1:
      total_sum_harmonicMean = 0
      for sampleN in range(1,Nsamples):
        total_sum_harmonicMean += 1 / sampleN
      a1_pop1 = float(total_sum_harmonicMean)
      w_theta_local = 1 / a1_pop1
      SS = 1
    else:
      w_theta_local = 0
      SS = 0
    # add value per window:
    if win!=current_win:
      if current_win == -1:
        current_win = win
        total_SS = SS
        W_theta = w_theta_local
        NVar = 1
        mean_Nsamples = Nsamples
      else:
        # Watterson Theta:
        total_W_theta = W_theta/winSize
        # produce output table
        output_popParameters.write(f'{chr}\t{current_win}\t{W_theta}\t{total_SS}\t{NVar}\t{total_W_theta}\t{mean_Nsamples/NVar}\n')
        current_win = win
        total_SS = SS
        W_theta = w_theta_local
        NVar = 1
        mean_Nsamples = Nsamples
    else:
      total_SS += SS
      W_theta += w_theta_local
      NVar += 1
      mean_Nsamples += Nsamples


output_popParameters.close()


```



Bash script: ThetaperInd.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=12:00:00
#SBATCH -J Pop_Par

genotype_table=$1
chr_used=$2
win_size=$3
prefix=$4
haplotypes=$5

#file="/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr01_missing_genotypeTable.txt"
#chr_used="chr01"
#winSize=int(10000)
#output="test"
#haplotypes =str("1,2,3,4")

# calculate population parameters:
python3 ThetaperInd.py $genotype_table $chr_used $win_size $prefix $haplotypes


```



job:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/05_diversity_withinInd/theta

mamba activate mypython3

# window size 10000
for i in $(seq 1 4 40)
do echo $i
sample=$(grep CHR /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr01_missing_genotypeTable.txt | awk -F"Alt\t" '{print $2}' | sed 's/\t/\n/g' | sed -n "$i"p | sed 's/_.*//g')
for chr in $(seq -w 1 12)
do echo $chr 
hap2=$(echo $i | awk '{print $1+1}')
hap3=$(echo $i | awk '{print $1+2}')
hap4=$(echo $i | awk '{print $1+3}')
# run job:
sbatch ThetaperInd.sh /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr$chr""_missing_genotypeTable.txt chr$chr 10000 TableTheta_Win10000_ConsideringmissingData_chr$chr"_Sample"$sample $i","$hap2","$hap3","$hap4
done
done
# 05.04.2024: 2453039..2453158

```



Randomising haplotypes:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/05_diversity_withinInd/theta/RandomSampling


mamba activate mypython3

for sample in $(seq -w 1 100)
do echo $sample
for chr in $(seq -w 1 12)
do echo $chr 
hap1=$(shuf -i 1-40 -n 1 )
hap2=$(shuf -i 1-40 -n 1 )
hap3=$(shuf -i 1-40 -n 1 )
hap4=$(shuf -i 1-40 -n 1 )
# run job:
sbatch ThetaperInd.sh /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr$chr""_missing_genotypeTable.txt chr$chr 10000 TableTheta_Win10000_ConsideringmissingData_chr$chr"_Sample"$sample $hap1","$hap2","$hap3","$hap4
done
done
# 05.04.2024: 2453159..2453926: actually only 64 replicates


cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/05_diversity_withinInd/theta

grep -v "mean_Nsamples" RandomSampling/TableTheta_Win10000_ConsideringmissingData_chr*Theta | sed 's@RandomSampling/TableTheta_Win10000_ConsideringmissingData_chr.._Sample@@g' | sed 's@.Theta:@\t@g' > RandomSampling_TableTheta_Win10000_ConsideringmissingData_AllChr_AllSamples.txt

grep -v "mean_Nsamples" TableTheta_Win10000_ConsideringmissingData_chr*Theta | sed 's@TableTheta_Win10000_ConsideringmissingData_chr.._Sample@@g' | sed 's@.Theta:@\t@g' > Cultivars_TableTheta_Win10000_ConsideringmissingData_AllChr_AllSamples.txt
# done

# for i in {2215213..2215555}; do stop $i ; rm slurm-$i* ; done

```



#### Plots:



```R
# .libPaths("/dss/dsshome1/lxc03/di36guz2/R")
# mamba activate env_others
.libPaths()

setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/05_diversity_withinInd/theta")

library(ggplot2)
library(tidyverse)

centromere_table <- read.csv("centromeres.csv", header = T)

chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

chrID <- "01"

# for (chrID in chr_list) {

table_RandomSampling <- read_table("RandomSampling_TableTheta_Win10000_ConsideringmissingData_AllChr_AllSamples.txt", F)
names(table_RandomSampling) <- c("sample", "chr", "pos", "sum_Theta", "SS", "NVar", "meanTheta", "meanNSamples")

table_cultivars <- read_table("Cultivars_TableTheta_Win10000_ConsideringmissingData_AllChr_AllSamples.txt", F)
names(table_cultivars) <- c("sample", "chr", "pos", "sum_Theta", "SS", "NVar", "meanTheta", "meanNSamples")




table_RandomSampling %>%
  group_by(sample, chr) %>%
  summarise(totalMeanTheta=mean(meanTheta)) %>%
  ggplot(aes(totalMeanTheta))+
  geom_histogram(fill="black") +
  geom_vline(data=table_cultivars %>%
               group_by(sample, chr) %>%
               summarise(totalMeanTheta=mean(meanTheta)) , 
             aes(xintercept = totalMeanTheta), color = "red", alpha=0.5) +
  labs(x = "Mean Theta/bp",
       y = "Count") +
  facet_grid(. ~ chr)+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))


ggsave(paste0("01_MeanTheta_perSample_perChr.png"), height=2, width=12)
ggsave(paste0("01_MeanTheta_perSample_perChr.pdf"), height=2, width=12)




winSize = 500000
table_RandomSampling %>%
  group_by(chr, pos) %>%
  summarise(maxTheta=max(meanTheta), 
            totalMeanTheta=mean(meanTheta, na.rm = T), 
            minTheta=min(meanTheta)) %>%
  mutate(win_ed=floor(pos/winSize)*winSize) %>%
  group_by(chr, win_ed) %>%
  summarise(MeanmaxTheta=mean(maxTheta, na.rm = T), 
            MeantotalMeanTheta=mean(totalMeanTheta, na.rm = T), 
            MeanminTheta=mean(minTheta, na.rm = T)) %>%
  filter(MeanmaxTheta<0.06)%>%
  ggplot(aes(win_ed/1000000, MeantotalMeanTheta))+
  geom_line(aes(win_ed/1000000, MeanminTheta), alpha=0.4)+
  geom_line(aes(win_ed/1000000, MeanmaxTheta), alpha=0.5)+
  geom_line(data=table_cultivars %>%
              mutate(win_ed=floor(pos/winSize)*winSize) %>%
              group_by(sample, chr, win_ed) %>%
              summarise(MeantotalMeanTheta=mean(meanTheta, na.rm = T))%>%
              filter(MeanmaxTheta<0.06), 
            aes(group=sample), colour="red", alpha=0.4)+
  geom_line(alpha=0.4, colour="blue")+
  geom_line(data=centromere_table %>%
              gather(key="posB", value = "pos", -chr) %>%
              mutate(posy = 0), 
            aes(pos/1000000,posy), colour="darkgreen", size=1.5, alpha=0.5) +
  scale_x_continuous(breaks=seq(0,100,15)) +
  labs(x = "Pos. (Mb)",
       y = "Theta") +
  facet_grid(. ~ chr, scale="free", space="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12), 
        axis.title = element_text(size=14))


ggsave("02_meanTheta_alongGenome_mean500kb.png", height=2, width=17)
ggsave("02_meanTheta_alongGenome_mean500kb.pdf", height=2, width=17)



```



#### Analyses removing singletons:



##### Remove singletons and subsampling genotype table:



python script: removeSingletons.py

```python
import numpy as np
import sys
import os

#file="/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr01_missing_genotypeTable.txt"
#output="filterSingletons_AllSamplesMultiAlleles_chr01_missing_genotypeTable.txt"


file = str(sys.argv[1])
output=str(sys.argv[2])

#output file 1:
output_table = open(f'{output}', "w")

def mac(List):
    unique, counts = np.unique(List, return_counts=True)
    return [min(counts), sum(counts)]

for line in open(file, "r"):
	if line[0:3] == "CHR":
		output_table.write(f'{line}')
	else:
		line_values  = line.strip().split("\t")
		chr, start, end, ref, alt = line_values[0:5]
		genotypes = [x for x in line_values[5:]]
		Non_nan_genotypes = [x != 'nan' for x in genotypes] 
		index_Non_nan_genotypes = [x for x in range(len(Non_nan_genotypes)) if Non_nan_genotypes[x]] 
		filtered_genotypes = [genotypes[i] for i in index_Non_nan_genotypes]
		if len(set(filtered_genotypes)) > 1:
			MAC, totalSample = mac(filtered_genotypes)
			if (MAC > 2) and (totalSample > 10):
				output_table.write(f'{line}')


output_table.close()

```



bash script: removeSingletons.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=4:00:00
#SBATCH -J removeSingletons

file=$1
output=$2

python3 removeSingletons.py $file $output 


```

job:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/05_diversity_withinInd/withoutSingletons

mamba activate mypython3

# /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr01_missing_genotypeTable.txt

# excluding mac < 3 and min10 haplotypes with no missing data:
for i in $(seq -w 1 12)
do chr="chr"$i ; echo $chr
sbatch removeSingletons.sh /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr$i"_"missing_genotypeTable.txt AllSamplesMultiAlleles_chr$i"_"missing_genotypeTable_NoSingletons.txt 
done
# 05.04.2024: 2453927..2453938

# for i in {2261814..2261825}; do stop $i ; rm slurm-$i* ; done
```





##### Calculate Theta



Script: ThetaperInd.py

```python
import numpy as np
from collections import defaultdict
from collections import Counter
import math
import sys
import os


#file="/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr01_missing_genotypeTable.txt"
#chr_used="chr01"
#winSize=int(10000)
#output="test"
#haplotypes =str("1,2,3,4")

file = str(sys.argv[1])
chr_used = str(sys.argv[2])
winSize=int(sys.argv[3])
output=str(sys.argv[4])
haplotypes = str(sys.argv[5])


haplotype_list= [int(x) for x in haplotypes.split(',')]

def round_up_to_nearest_win(num, win_size):
  return math.floor(int(num) / win_size) * win_size


#output file 1: population parameters per Window
output_popParameters = open(f'{output}.Theta', "w")

current_win = -1

for line in open(file, "r"):
  #chr, start, end, ref, alt, haplotypes  = line.strip().split("\t")
  if line[0:3] == "CHR":
    output_popParameters.write(f'chr\twin\tsum_theta\tSS\tNvar\tW_theta\tmean_Nsamples\n')
  else:
    line_values  = line.strip().split("\t")
    chr, start, end, ref, alt = line_values[0:5]
    win = round_up_to_nearest_win(start, winSize)
    genotypes =[int(line_values[i+4]) for i in haplotype_list if line_values[i+4] != 'nan']
    Nsamples = len(genotypes)
    # calculate theta per base:
    if len(set(genotypes)) > 1:
      total_sum_harmonicMean = 0
      for sampleN in range(1,Nsamples):
        total_sum_harmonicMean += 1 / sampleN
      #a1_pop1 = float(total_sum_harmonicMean)
      w_theta_local = 1 / total_sum_harmonicMean
      SS = 1
    else:
      w_theta_local = 0
      SS = 0
    # add value per window:
    if win!=current_win:
      if current_win == -1:
        current_win = win
        total_SS = SS
        W_theta = w_theta_local
        NVar = 1
        mean_Nsamples = Nsamples
      else:
        # Watterson Theta:
        total_W_theta = W_theta/winSize
        # produce output table
        output_popParameters.write(f'{chr}\t{current_win}\t{W_theta}\t{total_SS}\t{NVar}\t{total_W_theta}\t{mean_Nsamples/NVar}\n')
        current_win = win
        total_SS = SS
        W_theta = w_theta_local
        NVar = 1
        mean_Nsamples = Nsamples
    else:
      total_SS += SS
      W_theta += w_theta_local
      NVar += 1
      mean_Nsamples += Nsamples


output_popParameters.close()


```



Bash script: ThetaperInd.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=12:00:00
#SBATCH -J Pop_Par

genotype_table=$1
chr_used=$2
win_size=$3
prefix=$4
haplotypes=$5

#file="/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr01_missing_genotypeTable.txt"
#chr_used="chr01"
#winSize=int(10000)
#output="test"
#haplotypes =str("1,2,3,4")

# calculate population parameters:
python3 ThetaperInd.py $genotype_table $chr_used $win_size $prefix $haplotypes


```



job:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/05_diversity_withinInd/withoutSingletons

mamba activate mypython3

# window size 10000
for i in $(seq 1 4 40)
do echo $i
sample=$(grep CHR AllSamplesMultiAlleles_chr01_missing_genotypeTable_NoSingletons.txt | awk -F"Alt\t" '{print $2}' | sed 's/\t/\n/g' | sed -n "$i"p | sed 's/_.*//g')
for chr in $(seq -w 1 12)
do echo $chr 
hap2=$(echo $i | awk '{print $1+1}')
hap3=$(echo $i | awk '{print $1+2}')
hap4=$(echo $i | awk '{print $1+3}')
# run job:
sbatch ThetaperInd.sh AllSamplesMultiAlleles_chr$chr""_missing_genotypeTable_NoSingletons.txt chr$chr 10000 TableTheta_Win10000_ConsideringmissingData_NoSingletons_chr$chr"_Sample"$sample $i","$hap2","$hap3","$hap4
done
done
# 05.04.2024: 2454539..2454658

```



Randomising haplotypes:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/05_diversity_withinInd/withoutSingletons/RandomSampling

mamba activate mypython3

for sample in $(seq -w 1 50)
do echo $sample
for chr in $(seq -w 1 12)
do echo $chr 
hap1=$(shuf -i 1-40 -n 1 )
hap2=$(shuf -i 1-40 -n 1 )
hap3=$(shuf -i 1-40 -n 1 )
hap4=$(shuf -i 1-40 -n 1 )
# run job:
sbatch ThetaperInd.sh ../AllSamplesMultiAlleles_chr$chr""_missing_genotypeTable_NoSingletons.txt chr$chr 10000 TableTheta_Win10000_ConsideringmissingData_chr$chr"_Sample"$sample $hap1","$hap2","$hap3","$hap4
done
done
# 05.04.2024: 2454659..2455258


# sampling whole genome onces:
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/05_diversity_withinInd/withoutSingletons/RandomSampling_WG

mamba activate mypython3

for sample in $(seq -w 1 50)
do echo $sample
hap1=$(shuf -i 1-40 -n 1 )
hap2=$(shuf -i 1-40 -n 1 )
hap3=$(shuf -i 1-40 -n 1 )
hap4=$(shuf -i 1-40 -n 1 )
for chr in $(seq -w 1 12)
do echo $chr 
# run job:
sbatch ThetaperInd.sh ../AllSamplesMultiAlleles_chr$chr""_missing_genotypeTable_NoSingletons.txt chr$chr 10000 TableTheta_Win10000_ConsideringmissingData_chr$chr"_Sample"$sample $hap1","$hap2","$hap3","$hap4
done
done
# 05.04.2024: 2455259..2455858
#for i in {2215213..2215555}; do stop $i ; rm slurm-$i* ; done


cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/05_diversity_withinInd/withoutSingletons

grep -v "mean_Nsamples" RandomSampling_WG/TableTheta_Win10000_ConsideringmissingData_chr* | sed 's@RandomSampling_WG/TableTheta_Win10000_ConsideringmissingData_chr.._Sample@@g' | sed 's@.Theta:@\t@g' > RandomSampling_TableTheta_Win10000_ConsideringmissingData_AllChr_AllSamples.txt

grep -v "mean_Nsamples" TableTheta_Win10000_ConsideringmissingData_NoSingletons_chr* | sed 's@TableTheta_Win10000_ConsideringmissingData_NoSingletons_chr.._Sample@@g' | sed 's@.Theta:@\t@g' > Cultivars_TableTheta_Win10000_ConsideringmissingData_AllChr_AllSamples.txt
# done

# for i in {2215213..2215555}; do stop $i ; rm slurm-$i* ; done


```



##### Plots:



```R
# .libPaths("/dss/dsshome1/lxc03/di36guz2/R")
# mamba activate env_others
.libPaths()

setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/05_diversity_withinInd/withoutSingletons")

library(ggplot2)
library(tidyverse)

centromere_table <- read.csv("centromeres.csv", header = T)

chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

chrID <- "01"

# for (chrID in chr_list) {

table_RandomSampling <- read_table("RandomSampling_TableTheta_Win10000_ConsideringmissingData_AllChr_AllSamples.txt", F)
names(table_RandomSampling) <- c("sample", "chr", "pos", "sum_Theta", "SS", "NVar", "meanTheta", "meanNSamples")

table_cultivars <- read_table("Cultivars_TableTheta_Win10000_ConsideringmissingData_AllChr_AllSamples.txt", F)
names(table_cultivars) <- c("sample", "chr", "pos", "sum_Theta", "SS", "NVar", "meanTheta", "meanNSamples")




table_RandomSampling %>%
  group_by(sample, chr) %>%
  summarise(totalMeanTheta=mean(meanTheta)) %>%
  ggplot(aes(totalMeanTheta))+
  geom_histogram(fill="black") +
  geom_vline(data=table_cultivars %>%
               group_by(sample, chr) %>%
               summarise(totalMeanTheta=mean(meanTheta)) , 
             aes(xintercept = totalMeanTheta), color = "red", alpha=0.5) +
  labs(x = "Mean Theta/bp",
       y = "Count") +
  facet_grid(. ~ chr)+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"), 
        axis.text.x = element_text(angle=45, hjust=1))


ggsave(paste0("01_MeanTheta_perSample_perChr.png"), height=2, width=12)
ggsave(paste0("01_MeanTheta_perSample_perChr.pdf"), height=2, width=12)




winSize = 500000
table_RandomSampling %>%
  group_by(chr, pos) %>%
  summarise(maxTheta=max(meanTheta), 
            totalMeanTheta=mean(meanTheta, na.rm = T), 
            minTheta=min(meanTheta)) %>%
  mutate(win_ed=floor(pos/winSize)*winSize) %>%
  group_by(chr, win_ed) %>%
  summarise(MeanmaxTheta=mean(maxTheta, na.rm = T), 
            MeantotalMeanTheta=mean(totalMeanTheta, na.rm = T), 
            MeanminTheta=mean(minTheta, na.rm = T)) %>%
  filter(MeanmaxTheta<0.06)%>%
  ggplot(aes(win_ed/1000000, MeantotalMeanTheta))+
  geom_line(aes(win_ed/1000000, MeanminTheta), alpha=0.4)+
  geom_line(aes(win_ed/1000000, MeanmaxTheta), alpha=0.5)+
  geom_line(data=table_cultivars %>%
              mutate(win_ed=floor(pos/winSize)*winSize) %>%
              group_by(sample, chr, win_ed) %>%
              summarise(MeantotalMeanTheta=mean(meanTheta, na.rm = T)), 
            aes(group=sample), colour="red", alpha=0.4)+
  geom_line(alpha=0.4, colour="blue")+
  geom_line(data=centromere_table %>%
              gather(key="posB", value = "pos", -chr) %>%
              mutate(posy = 0), 
            aes(pos/1000000,posy), colour="darkgreen", size=1.5, alpha=0.5) +
  scale_x_continuous(breaks=seq(0,100,15)) +
  labs(x = "Pos. (Mb)",
       y = "Theta") +
  facet_grid(. ~ chr, scale="free", space="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12), 
        axis.title = element_text(size=14))


ggsave("02_meanTheta_alongGenome_mean500kb.png", height=2, width=17)
ggsave("02_meanTheta_alongGenome_mean500kb.pdf", height=2, width=17)





table_cultivars %>%
  group_by(sample) %>%
  summarise(totalMeanTheta=mean(meanTheta)) 


table_RandomSampling %>%
  group_by(sample) %>%
  summarise(totalMeanTheta=mean(meanTheta)) %>%
  ggplot(aes(totalMeanTheta))+
  geom_histogram(fill="black") +
  geom_vline(data=table_cultivars %>%
               group_by(sample) %>%
               summarise(totalMeanTheta=mean(meanTheta)) , 
             aes(xintercept = totalMeanTheta), color = "red", alpha=0.5) +
  labs(x = "Mean Theta/bp",
       y = "Count") +
  # facet_grid(. ~ chr)+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"), 
        axis.text.x = element_text(angle=45, hjust=1))


ggsave(paste0("03_MeanTheta_perSample.png"), height=4, width=4)
ggsave(paste0("03_MeanTheta_perSample.pdf"), height=4, width=4)







table_RandomSampling %>%
  group_by(sample, chr) %>%
  summarise(totalSS=sum(SS)) %>%
  ggplot(aes(totalSS))+
  geom_histogram(fill="black") +
  geom_vline(data=table_cultivars %>%
               group_by(sample, chr) %>%
               summarise(totalSS=sum(SS)), 
             aes(xintercept = totalSS), color = "red", alpha=0.5) +
  labs(x = "Total SS",
       y = "Count") +
  facet_grid(. ~ chr)+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"), 
        axis.text.x = element_text(angle=45, hjust=1))


ggsave(paste0("04_SS_perSample_perChr.png"), height=2, width=12)
ggsave(paste0("04_SS_perSample_perChr.pdf"), height=2, width=12)







table_RandomSampling %>%
  group_by(sample) %>%
  summarise(totalSS=sum(SS)) %>%
  ggplot(aes(totalSS))+
  geom_histogram(fill="black") +
  geom_vline(data=table_cultivars %>%
               group_by(sample) %>%
               summarise(totalSS=sum(SS)), 
             aes(xintercept = totalSS), color = "red", alpha=0.5) +
  labs(x = "Total SS",
       y = "Count") +
  # facet_grid(. ~ chr)+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"), 
        axis.text.x = element_text(angle=45, hjust=1))


ggsave(paste0("05_SS_perSample.png"), height=2, width=12)
ggsave(paste0("05_SS_perSample.pdf"), height=2, width=12)



```







# Number of variant sites per window



```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes

# file to work with:
grep -vw CHR AllSamplesMultiAlleles_chr*_genotypeTable_miltipleAlleles.txt | sed 's/:/\t/g' | awk '{Win=int($3/10000) ; print $2"\t"Win}' | sort | uniq -c | sed 's/^[ \t]*//' | awk -F' ' '{print $2"\t"$3"\t"$1}' > VarCount_AllSamplesMultiAlleles_Allchr_genotypeTable_miltipleAlleles.txt

```

Plot:

```R
# .libPaths("/dss/dsshome1/lxc03/di36guz2/R")
.libPaths()

setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes")

library(ggplot2)
library(tidyverse)


table_VarN <- read_table("VarCount_AllSamplesMultiAlleles_Allchr_genotypeTable_miltipleAlleles.txt", F)
names(table_VarN) <- c("chr", "pos", "NVar")

centromere_table <- read.csv("centromeres.csv", header = T)


head(table_VarN)
centromere_table


# edited version for manuscript:
winSize = 500000
table_VarN %>%
  mutate(win_ed=floor(pos*10000/winSize)*winSize) %>%
  group_by(chr, win_ed) %>%
  summarise(meanNVar=mean(NVar)) %>%
  filter(meanNVar<1300) %>%
  ggplot(aes(win_ed/1000000, meanNVar))+
  geom_line(alpha=0.8)+
  geom_line(data=centromere_table %>%
              gather(key="posB", value = "pos", -chr) %>%
              mutate(posy = 0), 
            aes(pos/1000000,posy), colour="darkgreen", size=1.5, alpha=0.5) +
  #geom_point(alpha=0.04)+
  geom_hline(yintercept=mean(table_VarN$NVar), linetype="dashed", color = "red")+
  scale_y_continuous(breaks=seq(0,1500,250)) +
  scale_x_continuous(breaks=seq(0,100,15)) +
  facet_grid(. ~ chr, scale="free", space="free")+
  labs(x = "Pos. (Mb)",
       y = "N. Var") +
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12), 
        axis.title = element_text(size=14))

ggsave("07_ed_NVar_AllSamples_mean500kb.png", height=2, width=15)
ggsave("07_ed_NVar_AllSamples_mean500kb.pdf", height=2, width=15)






table_VarN_mean <- table_VarN %>%
  group_by(chr) %>%
  summarise(mean_NVar=mean(NVar))

table_VarN %>%
  ggplot(aes(NVar))+
  geom_histogram(fill="black")+
  geom_vline(data=table_VarN_mean, aes(xintercept = mean_NVar), 
             #linetype="dotted", 
             color = "red", alpha=0.5) +
  geom_vline(xintercept = mean(table_VarN$NVar), 
             #linetype="dotted", 
             color = "blue", alpha=0.5) +
  scale_x_continuous(breaks=seq(0,5000,1000)) +
  labs(x = "N. Var/10kb Win",
       y = "Count") +
  facet_grid(.~chr)+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"), 
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12, hjust = 1, angle = 45))

ggsave(paste0("08_Distribution_NVar_10kbWin_perChr.png"), height=2.5, width=15)
ggsave(paste0("08_Distribution_NVar_10kbWin_perChr.pdf"), height=2.5, width=15)



table_VarN %>%
  ggplot(aes(NVar))+
  geom_histogram(fill="black", binwidth = 50)+
  geom_vline(xintercept = mean(table_VarN$NVar), 
             #linetype="dotted", 
             color = "red") +
  labs(x = "N. Var/10kb Win",
       y = "Count") +
  # facet_grid(chr~., scale="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12, hjust = 1, angle = 45))

ggsave(paste0("09_Distribution_NVar_10kbWin_allChr.png"), height=4, width=5)
ggsave(paste0("09_Distribution_NVar_10kbWin_allChr.pdf"), height=4, width=5)


table_VarN %>%
  filter(NVar<2500) %>%
  ggplot(aes(NVar))+
  geom_histogram(colour="black", fill="black", binwidth = 20)+
  geom_vline(xintercept = mean(table_VarN$NVar), 
             #linetype="dotted", 
             color = "red") +
  labs(x = "N. Var/10kb Window",
       y = "Count") +
  # facet_grid(chr~., scale="free")+
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.title = element_text(size=16), 
        axis.text = element_text(size=14))


ggsave(paste0("09_v2_Distribution_NVar_10kbWin_allChr.png"), height=2.6, width=3.9)
ggsave(paste0("09_v2_Distribution_NVar_10kbWin_allChr.pdf"), height=2.6, width=3.9)

```



# Distribution non-aligned regions



script: meanCovNOTAL_reg.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00
#SBATCH -J covNOTAL_reg

winSize=$1
list_samples=$2
prefix=$3
chr=$4

for haplotype_file in $(cat $list_samples )
do echo $haplotype_file

haplotype_name=$(echo $haplotype_file | sed 's@.*/01_nucmer_syri_vs_DM/nucmer_syri_@@g' | sed 's@vs_DM/refDM_chr.*_vs_qry_@@g' | sed 's@/syri.out@@g' )
echo $haplotype_name
echo $chr

awk '{if($11=="HDR" && $3-($2-1)>100) print $1"\t"$2-1"\t"$3}' $haplotype_file > $prefix"_"hdr_$chr"_"$haplotype_name.txt

awk '{if($11=="CPL" && $3-($2-1)>100) print $1"\t"$2-1"\t"$3}' $haplotype_file > $prefix"_"cpl_$chr"_"$haplotype_name.txt

awk '{if($1!="-" && $11=="NOTAL" && $3-($2-1)>100) print $1"\t"$2-1"\t"$3}' $haplotype_file > $prefix"_"NOTAL_$chr"_"$haplotype_name.txt

cat $prefix"_"NOTAL_$chr"_"$haplotype_name.txt $prefix"_"hdr_$chr"_"$haplotype_name.txt $prefix"_"cpl_$chr"_"$haplotype_name.txt  | sort -k1,1 -k2,2n > $prefix"_"ALLMissing_$chr"_"$haplotype_name.txt

bedtools merge -i $prefix"_"ALLMissing_$chr"_"$haplotype_name.txt > $prefix"_"$haplotype_name.bed

rm $prefix"_"hdr_$chr"_"$haplotype_name.txt $prefix"_"cpl_$chr"_"$haplotype_name.txt $prefix"_"NOTAL_$chr"_"$haplotype_name.txt $prefix"_"ALLMissing_$chr"_"$haplotype_name.txt

#awk '{if($1!="-" && $11=="NOTAL" && $3-($2-1)>100) print $1"\t"$2-1"\t"$3}' $haplotype_file | sort -k1,1 -k2,2n > $prefix"_"$haplotype_name.bed

done

cat $prefix"_"*hap*.bed | sort -k1,1 -k2,2n > $prefix"_All".bed

chr_size=$(awk 'BEGIN{maxVal=0} {if($2>$3) {line_max=$2} else {line_max=$3} ; if (line_max>maxVal) maxVal=line_max} END {print maxVal}' $prefix"_All".bed )

rm Win_$chr.bed
for start in $(seq 0 $winSize $chr_size )
do end=$(echo $start | awk -v WS=$winSize '{print $1+WS-1}' )
echo -e $chr'\t'$start'\t'$end >> Win_$chr.bed
done


bedtools coverage -a Win_$chr.bed -b $prefix"_All".bed -d | awk -v CHR=$chr 'BEGIN {win=0 ; SumCov=0 ; NLines=0} {if($2==win) {SumCov+=$5 ; NLines+=1} else {print CHR"\t"win"\t"SumCov/NLines ; win=$2 ; SumCov=$5 ; NLines=1} } END {print CHR"\t"win"\t"SumCov/NLines}' > meanCoverage_$prefix.txt


rm $prefix"_"*hap*.bed $prefix"_All".bed Win_$chr.bed

```



Job calculating coverage of non-aligned regions along the genome:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes

mamba activate env_others

winSize=500000
for i in $(seq -w 1 12)
do echo $i
list_samples=list_haplotypes_chr$i"".txt
prefix=AllSamplesMultiAlleles_chr$i""_NOTAL
chr=chr$i""
sbatch meanCovNOTAL_reg.sh $winSize $list_samples $prefix $chr
done
# 05.04.2024: 2455871..2455882

cat meanCoverage_AllSamplesMultiAlleles_chr*_NOTAL.txt > all_meanCoverage_AllSamplesMultiAlleles_NOTAL.txt
# done

#for i in {2455859..2455870}; do stop $i ; done


```



## Plot

```R
# .libPaths("/dss/dsshome1/lxc03/di36guz2/R")
.libPaths()

setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes")


library(ggplot2)
library(tidyverse)


table_NOTAL <- read_table("all_meanCoverage_AllSamplesMultiAlleles_NOTAL.txt", F)
names(table_NOTAL) <- c("chr", "pos", "MeanNOTALSamples")

centromere_table <- read.csv("centromeres.csv", header = T)


head(table_NOTAL)
centromere_table


table_NOTAL %>%
  ggplot(aes(pos/1000000, MeanNOTALSamples))+
  geom_line(alpha=0.8)+
  geom_line(data=centromere_table %>%
              gather(key="posB", value = "pos", -chr) %>%
              mutate(posy = 0), 
            aes(pos/1000000,posy), colour="darkgreen", size=1.5, alpha=0.5) +
  #geom_point(alpha=0.04)+
  geom_hline(yintercept=median(table_NOTAL$MeanNOTALSamples), linetype="dashed", color = "red")+
  # scale_y_continuous(breaks=seq(0,1500,250)) +
  scale_x_continuous(breaks=seq(0,100,15)) +
  facet_grid(. ~ chr, scale="free", space="free")+
  labs(x = "Pos. (Mb)",
       y = "N. Not_Al. Samples") +
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12), 
        axis.title = element_text(size=14))

ggsave("10_ed_NNotAlSamples_mean500kb.png", height=2, width=15)
ggsave("10_ed_NNotAlSamples_mean500kb.pdf", height=2, width=15)

```



# Number of segregating sites and haplotypes per window:



## Including nan as different genotypes:

This version uses missing data as a different genotype, which is coded as *-1* genotype.

script: calculateSS.py

```python
import os
import sys
import numpy as np
import math

#genotype_file = str("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr01_missing_genotypeTable.txt")
#win_size = int(10000)
#output_name = str("test")

genotype_file = str(sys.argv[1])
win_size = int(sys.argv[2])
output_name = str(sys.argv[3])

output_file = open(f'{output_name}.txt', "w")

def str_to_int(x):
        if x!='nan':
                return int(x)
        else:
                return int(-1)

for line in open(genotype_file, "r"):
        if line[0:3] == "CHR":
                haplotypes = line.strip().split("\t")[5:]
                SS_array = np.zeros((len(haplotypes), len(haplotypes)))
                NVar = 0
                win_current = 0
        else:
                line_values = line.strip().split("\t")
                pos = int(line_values[1])
                win = math.ceil(int(pos) / win_size) * win_size
                genotypes = np.asarray([str_to_int(x) for x in line_values[5:]])
                SS_var = genotypes[:,None] != genotypes
                if win == win_current:
                        SS_array += SS_var
                        NVar += 1
                else:
                        for i,x in enumerate(haplotypes):
                                for j,y in enumerate(haplotypes):
                                        if i > j:
                                                output_file.write(f'{win_current}\t{x}\t{y}\t{int(SS_array[i,j])}\n')
                        SS_array = SS_var*1
                        NVar = 0
                        win_current = win

output_file.close()

```

script: calculateSS.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00
#SBATCH -J calculateSS

genotype_file=$1
win_size=$2
output_name=$3

# python3 /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr01_missing_genotypeTable.txt 10000 NSS_chr01
python calculateSS.py $genotype_file $win_size $output_name

```

job:

This script produces a table with pairwise comparisons between haplotypes and windows and report the number of segregating sites. 

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/03_Nhaplotypes

mamba activate env_others

# 10 kb WIndow size 
for chN in $( seq -w 1 12); do echo $chN
sbatch calculateSS.sh /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr$chN"_"missing_genotypeTable.txt 10000 NSS_10kbW_chr$chN
done
# 05.04.2024: 2455883..2455894

for chN in $( seq -w 1 12); do echo $chN
awk '{print $4}' NSS_10kbW_chr$chN.txt | sort | uniq -c | sed 's/^[ \t]*//' | awk -F' ' '{print $2"\t"$1}' > Total_NSS_10kbW_chr$chN.txt
done

grep "" Total_NSS_10kbW_chr*.txt | sed 's@Total_NSS_10kbW_chr\(..\).txt:@\1\t@g' > all_Total_NSS_10kbW.txt
# done

# 100 kb WIndow size 
for chN in $( seq -w 1 12); do echo $chN
sbatch calculateSS.sh /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr$chN"_"missing_genotypeTable.txt 100000 NSS_100kbW_chr$chN
done
# 05.04.2024: 2455895..2455906

for chN in $( seq -w 1 12); do echo $chN
awk '{print $4}' NSS_100kbW_chr$chN.txt | sort | uniq -c | sed 's/^[ \t]*//' | awk -F' ' '{print $2"\t"$1}' > Total_NSS_100kbW_chr$chN.txt
done

grep "" Total_NSS_100kbW_chr*.txt | sed 's@Total_NSS_100kbW_chr\(..\).txt:@\1\t@g' > all_Total_NSS_100kbW.txt
# done

# 200 kb WIndow size 
for chN in $( seq -w 1 12); do echo $chN
sbatch calculateSS.sh /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr$chN"_"missing_genotypeTable.txt 200000 NSS_200kbW_chr$chN
done
# 05.04.2024: 2455907..2455918

for chN in $( seq -w 1 12); do echo $chN
awk '{print $4}' NSS_200kbW_chr$chN.txt | sort | uniq -c | sed 's/^[ \t]*//' | awk -F' ' '{print $2"\t"$1}' > Total_NSS_200kbW_chr$chN.txt
done

grep "" Total_NSS_200kbW_chr*.txt | sed 's@Total_NSS_200kbW_chr\(..\).txt:@\1\t@g' > all_Total_NSS_200kbW.txt
# done

# 500 kb WIndow size 
for chN in $( seq -w 1 12); do echo $chN
sbatch calculateSS.sh /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr$chN"_"missing_genotypeTable.txt 500000 NSS_500kbW_chr$chN
done
# 05.04.2024: 2455919..2455930

for chN in $( seq -w 1 12); do echo $chN
awk '{print $4}' NSS_500kbW_chr$chN.txt | sort | uniq -c | sed 's/^[ \t]*//' | awk -F' ' '{print $2"\t"$1}' > Total_NSS_500kbW_chr$chN.txt
done

grep "" Total_NSS_500kbW_chr*.txt | sed 's@Total_NSS_500kbW_chr\(..\).txt:@\1\t@g' > all_Total_NSS_500kbW.txt
# done


# 1 Mb WIndow size 
for chN in $( seq -w 1 12); do echo $chN
sbatch calculateSS.sh /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr$chN"_"missing_genotypeTable.txt 1000000 NSS_1MbW_chr$chN
done
# 05.04.2024: 2455931..2455942

for chN in $( seq -w 1 12); do echo $chN
awk '{print $4}' NSS_1MbW_chr$chN.txt | sort | uniq -c | sed 's/^[ \t]*//' | awk -F' ' '{print $2"\t"$1}' > Total_NSS_1MbW_chr$chN.txt
done

grep "" Total_NSS_1MbW_chr*.txt | sed 's@Total_NSS_1MbW_chr\(..\).txt:@\1\t@g' > all_Total_NSS_1MbW.txt
# done

```



### Plots 10 kb Windows



```R
# .libPaths("/dss/dsshome1/lxc03/di36guz2/R")
# mamba activate env_others
.libPaths()

setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/03_Nhaplotypes")

library(ggplot2)
library(tidyverse)

chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

chr <- "03"

for (chr in chr_list) {
  table_CountSS <- read_table(paste0("Total_NSS_10kbW_chr",chr,".txt"), F)
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


table_CountSS <- read_table("all_Total_NSS_10kbW.txt", F)
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

ggsave("04_distribution_pairwiseSS_per10kbWin_perWindow_allchr_allNSS.pdf", height=12, width=6)
ggsave("04_distribution_pairwiseSS_per10kbWin_perWindow_allchr_2to500NSS.png", height=12, width=6)
ggsave("04_distribution_pairwiseSS_per10kbWin_perWindow_allchr_2to500NSS.pdf", height=12, width=6)




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




```



## Excluding nan from the comparison:

In this version nan are not considered in the comparison between haplotypes. 

script: calculateSS_removingNAN.py

```python
import os
import sys
import numpy as np
import math

#genotype_file = str("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr01_missing_genotypeTable.txt")
#win_size = int(10000)
#output_name = str("test")

genotype_file = str(sys.argv[1])
win_size = int(sys.argv[2])
output_name = str(sys.argv[3])

output_file = open(f'{output_name}.txt', "w")

def str_to_int(x):
        if x!='nan':
                return int(x)
        else:
                return int(-1)

for line in open(genotype_file, "r"):
        if line[0:3] == "CHR":
                haplotypes = line.strip().split("\t")[5:]
                SS_array = np.zeros((len(haplotypes), len(haplotypes)))
                NVar = 0
                win_current = 0
        else:
                line_values = line.strip().split("\t")
                pos = int(line_values[1])
                win = math.ceil(int(pos) / win_size) * win_size
                genotypes = np.asarray([str_to_int(x) for x in line_values[5:]])
                SS_var = genotypes[:,None] != genotypes
                NA_genotypes_boolean = (genotypes==(-1))
                SS_var[NA_genotypes_boolean, :] = False
                SS_var[:, NA_genotypes_boolean] = False
                if win == win_current:
                        SS_array += SS_var
                        NVar += 1
                else:
                        for i,x in enumerate(haplotypes):
                                for j,y in enumerate(haplotypes):
                                        if i > j:
                                                output_file.write(f'{win_current}\t{x}\t{y}\t{int(SS_array[i,j])}\n')
                        SS_array = SS_var*1
                        NVar = 0
                        win_current = win

output_file.close()


```

script: calculateSS_removingNAN.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00
#SBATCH -J calculateSS

genotype_file=$1
win_size=$2
output_name=$3

# python3 /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr01_missing_genotypeTable.txt 10000 NSS_chr01
python calculateSS_removingNAN.py $genotype_file $win_size $output_name

```

job:

This script produces a table with pairwise comparisons between haplotypes and windows and report the number of segregating sites. 

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/03_Nhaplotypes

mamba activate env_others

# 10 kb WIndow size 
for chN in $( seq -w 1 12); do echo $chN
sbatch calculateSS_removingNAN.sh /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr$chN"_"missing_genotypeTable.txt 10000 removingNAN_NSS_10kbW_chr$chN
done
# 05.04.2024: 2455943..2455954

# 100 kb WIndow size 
for chN in $( seq -w 1 12); do echo $chN
sbatch calculateSS_removingNAN.sh /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr$chN"_"missing_genotypeTable.txt 100000 removingNAN_NSS_100kbW_chr$chN
done
# 

# 200 kb WIndow size 
for chN in $( seq -w 1 12); do echo $chN
sbatch calculateSS_removingNAN.sh /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr$chN"_"missing_genotypeTable.txt 200000 removingNAN_NSS_200kbW_chr$chN
done
# 

# 500 kb WIndow size 
for chN in $( seq -w 1 12); do echo $chN
sbatch calculateSS_removingNAN.sh /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr$chN"_"missing_genotypeTable.txt 500000 removingNAN_NSS_500kbW_chr$chN
done
# 

# 1 Mb WIndow size 
for chN in $( seq -w 1 12); do echo $chN
sbatch calculateSS_removingNAN.sh /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr$chN"_"missing_genotypeTable.txt 1000000 removingNAN_NSS_1MbW_chr$chN
done
# 

# 06.04.2024: 2456647..2456647 all the above


for chN in $( seq -w 1 12); do echo $chN
awk '{print $4}' removingNAN_NSS_10kbW_chr$chN.txt | sort | uniq -c | sed 's/^[ \t]*//' | awk -F' ' '{print $2"\t"$1}' > Total_removingNAN_NSS_10kbW_chr$chN.txt
done

grep "" Total_removingNAN_NSS_10kbW_chr*.txt | sed 's@Total_removingNAN_NSS_10kbW_chr\(..\).txt:@\1\t@g' > all_Total_removingNAN_NSS_10kbW.txt


for chN in $( seq -w 1 12); do echo $chN
awk '{print $4}' removingNAN_NSS_100kbW_chr$chN.txt | sort | uniq -c | sed 's/^[ \t]*//' | awk -F' ' '{print $2"\t"$1}' > Total_removingNAN_NSS_100kbW_chr$chN.txt
done

grep "" Total_removingNAN_NSS_100kbW_chr*.txt | sed 's@Total_removingNAN_NSS_100kbW_chr\(..\).txt:@\1\t@g' > all_Total_removingNAN_NSS_100kbW.txt


for chN in $( seq -w 1 12); do echo $chN
awk '{print $4}' removingNAN_NSS_200kbW_chr$chN.txt | sort | uniq -c | sed 's/^[ \t]*//' | awk -F' ' '{print $2"\t"$1}' > Total_removingNAN_NSS_200kbW_chr$chN.txt
done

grep "" Total_removingNAN_NSS_200kbW_chr*.txt | sed 's@Total_removingNAN_NSS_200kbW_chr\(..\).txt:@\1\t@g' > all_Total_removingNAN_NSS_200kbW.txt


for chN in $( seq -w 1 12); do echo $chN
awk '{print $4}' removingNAN_NSS_500kbW_chr$chN.txt | sort | uniq -c | sed 's/^[ \t]*//' | awk -F' ' '{print $2"\t"$1}' > Total_removingNAN_NSS_500kbW_chr$chN.txt
done

grep "" Total_removingNAN_NSS_500kbW_chr*.txt | sed 's@Total_removingNAN_NSS_500kbW_chr\(..\).txt:@\1\t@g' > all_Total_removingNAN_NSS_500kbW.txt


for chN in $( seq -w 1 12); do echo $chN
awk '{print $4}' removingNAN_NSS_1MbW_chr$chN.txt | sort | uniq -c | sed 's/^[ \t]*//' | awk -F' ' '{print $2"\t"$1}' > Total_removingNAN_NSS_1MbW_chr$chN.txt
done

grep "" Total_removingNAN_NSS_1MbW_chr*.txt | sed 's@Total_removingNAN_NSS_1MbW_chr\(..\).txt:@\1\t@g' > all_Total_removingNAN_NSS_1MbW.txt


```



### Other window sizes plots:



```R
# .libPaths("/dss/dsshome1/lxc03/di36guz2/R")
# mamba activate env_others
.libPaths()

setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/03_Nhaplotypes")

library(ggplot2)
library(tidyverse)

chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

chr <- "01"

# WinS <- "500kbW"
WinS <- "1MbW"

for (chr in chr_list) {
  table_CountSS <- read_table(paste0("Total_NSS_", WinS, "_chr",chr,".txt"), F)
  names(table_CountSS) <- c("NSS", "Count")
  
  table_CountSS %>%
    ggplot(aes(NSS, Count))+
    geom_point()+
    theme_classic()+
    theme(panel.grid.major = element_line(color = "grey80"),
          panel.grid.minor = element_line(color = "grey80"))
  
  ggsave(paste0("01_distribution_pairwiseSS_per", WinS, "_perWindow_chr", chr,".png"), height=7, width=8)
  
  table_CountSS %>%
    filter(NSS>2) %>%
    filter(NSS<500) %>%
    ggplot(aes(NSS, Count))+
    geom_point()+
    geom_line() +
    theme_classic()+
    theme(panel.grid.major = element_line(color = "grey80"),
          panel.grid.minor = element_line(color = "grey80"))
  
  ggsave(paste0("02_distribution_pairwiseSS_per", WinS, "_perWindow_chr", chr,"_2to500NSS.png"), height=7, width=8)
  
  # table_CountSS %>%
  #   filter(NSS>2) %>%
  #   filter(NSS<300) %>%
  #   ggplot(aes(NSS, Count))+
  #   geom_point()+
  #   geom_line() +
  #   theme_classic()+
  #   theme(panel.grid.major = element_line(color = "grey80"),
  #         panel.grid.minor = element_line(color = "grey80"))
  # 
  # ggsave(paste0("03_distribution_pairwiseSS_per", WinS, "_perWindow_chr", chr,"_2to100NSS.png"), height=7, width=8)
}


table_CountSS <- read_table(paste0("all_Total_NSS_", WinS, ".txt"), F)
names(table_CountSS) <- c("chr", "NSS", "Count")

table_CountSS %>%
  ggplot(aes(NSS, Count/1000))+
  geom_point()+
  geom_line() +
  labs(x = "NSS",
       y = "Count (10^3)") +
  facet_grid(chr ~ ., scale="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))

ggsave(paste0("04_distribution_pairwiseSS_per", WinS, "_perWindow_allchr_allNSS.pdf"), height=12, width=6)


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

ggsave(paste0("04_distribution_pairwiseSS_per", WinS, "_perWindow_allchr_2to500NSS.png"), height=12, width=6)
ggsave(paste0("04_distribution_pairwiseSS_per", WinS, "_perWindow_allchr_2to500NSS.pdf"), height=12, width=6)




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

ggsave(paste0("04_ed_distribution_pairwiseSS_per", WinS, "_perWindow_all.png"), height=5, width=7)


table_CountSS %>% 
  filter(NSS>4) %>%
  filter(NSS<500) %>%
  group_by(NSS) %>%
  summarise(CountTotal=sum(Count)) %>%
  ggplot(aes(NSS, CountTotal/1000))+
  geom_point()+
  geom_line() +
  labs(x = "NSS",
       y = "Count (10^3)") +
  # facet_grid(chr ~ ., scale="free")+
  scale_x_continuous(breaks=c(5, seq(50,1500,50))) +
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor = element_line(color = "grey80"),
        axis.text = element_text(size=16), 
        axis.title = element_text(size=16))

ggsave(paste0("04_ed_distribution_pairwiseSS_per", WinS, "_perWindow_all_5to300NSS.png"), height=5, width=7)



table_CountSS %>% 
  filter(NSS>4) %>%
  filter(NSS<50) %>%
  group_by(NSS) %>%
  summarise(CountTotal=sum(Count)) %>%
  ggplot(aes(NSS, CountTotal/1000))+
  geom_point()+
  geom_line() +
  labs(x = "NSS",
       y = "Count (10^3)") +
  # facet_grid(chr ~ ., scale="free")+
  scale_x_continuous(breaks=c(5, seq(0,500,5))) +
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor = element_line(color = "grey80"),
        axis.text = element_text(size=24), 
        axis.title = element_text(size=24))






# table_CountSS %>%
#   filter(NSS>2) %>%
#   filter(NSS<50) %>%
#   ggplot(aes(NSS, Count/1000))+
#   geom_point()+
#   geom_line() +
#   labs(x = "NSS",
#        y = "Count (10^3)") +
#   facet_grid(chr ~ ., scale="free")+
#   theme_classic()+
#   theme(panel.grid.major = element_line(color = "grey80"),
#         panel.grid.minor = element_line(color = "grey80"))
# 
# ggsave(paste0("05_distribution_pairwiseSS_per", WinS, "_perWindow_allchr_2to50NSS.png"), height=12, width=6)
# 





# plots removing missing data:


WinS <- "500kbW"
WinS <- "1MbW"


table_CountSS <- read_table(paste0("all_Total_removingNAN_NSS_", WinS, ".txt"), F)
names(table_CountSS) <- c("chr", "NSS", "Count")

table_CountSS %>%
  ggplot(aes(NSS, Count/1000))+
  geom_point()+
  geom_line() +
  labs(x = "NSS",
       y = "Count (10^3)") +
  facet_grid(chr ~ ., scale="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))

ggsave(paste0("04_distribution_pairwiseSS_per", WinS, "_perWindow_allchr_allNSS_removingNAN.pdf"), height=12, width=6)


# table_CountSS %>%
#   filter(NSS>2) %>%
#   filter(NSS<500) %>%
#   ggplot(aes(NSS, Count/1000))+
#   geom_point()+
#   geom_line() +
#   labs(x = "NSS",
#        y = "Count (10^3)") +
#   facet_grid(chr ~ ., scale="free")+
#   theme_classic()+
#   theme(panel.grid.major = element_line(color = "grey80"),
#         panel.grid.minor = element_line(color = "grey80"))
# 
# ggsave(paste0("04_distribution_pairwiseSS_per", WinS, "_perWindow_allchr_2to500NSS_removingNAN.png"), height=12, width=6)
# ggsave(paste0("04_distribution_pairwiseSS_per", WinS, "_perWindow_allchr_2to500NSS_removingNAN.pdf"), height=12, width=6)
# 



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

ggsave(paste0("04_ed_distribution_pairwiseSS_per", WinS, "_perWindow_all_removingNAN.png"), height=5, width=7)


table_CountSS %>% 
  filter(NSS>4) %>%
  filter(NSS<500) %>%
  group_by(NSS) %>%
  summarise(CountTotal=sum(Count)) %>%
  ggplot(aes(NSS, CountTotal/1000))+
  geom_point()+
  geom_line() +
  labs(x = "NSS",
       y = "Count (10^3)") +
  # facet_grid(chr ~ ., scale="free")+
  scale_x_continuous(breaks=c(5, seq(50,1500,50))) +
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor = element_line(color = "grey80"),
        axis.text = element_text(size=16), 
        axis.title = element_text(size=16))

ggsave(paste0("04_ed_distribution_pairwiseSS_per", WinS, "_perWindow_all_5to300NSS_removingNAN.png"), height=5, width=7)



```



## Clustering haplotypes



script: ClusteringHaplotypes.py

```python
import os
import sys
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as shc
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
import math

# genotype_file = str("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr01_missing_genotypeTable.txt")
# win_size = int(10000)
# output_name = str("test")
# max_d = int(10)

# this script takes a genotype table, a window size and and output name and: produce a distance matrix and cluster haplotypes with a maximum of around 10SNPs (max_d parameter).
# the output is a table with window start position, sample and haplotype ID after clustering.
# The scripts uses the haplotype groups produced from the previous window to allocate haplotypes ID in the current window. 

genotype_file = str(sys.argv[1])
win_size = int(sys.argv[2])
output_name = str(sys.argv[3])
max_d = int(sys.argv[4])
#max_d = 10

def str_to_int(x):
    if x!='nan':
        return int(x)
    else:
        return int(-1)


for line in open(genotype_file, "r"):
    if line[0:3] == "CHR":
        haplotypes = line.strip().split("\t")[5:]
        SS_array = np.zeros((len(haplotypes), len(haplotypes)))
        NA_genotypes_mean = np.zeros((len(haplotypes)))
        NVar = 0
        win_current = 0
        table_hap_groups_pre = pd.DataFrame()
        if os.path.exists(f"{output_name}.txt"):
            os.remove(f"{output_name}.txt")
    else:
        line_values = line.strip().split("\t")
        pos = int(line_values[1])
        win = int(pos / win_size)
        genotypes = np.asarray([str_to_int(x) for x in line_values[5:]])
        SS_var = genotypes[:,None] != genotypes
        NA_genotypes = (genotypes==(-1))*1
        if win == win_current:
            SS_array += SS_var
            NA_genotypes_mean += NA_genotypes
            NVar += 1
        else:
            #Z = shc.ward(pdist(SS_array))
            #clusters = shc.fcluster(Z, max_d, criterion='distance')
            condensed_dist = squareform(SS_array)
            linkresult = shc.linkage(condensed_dist, method='complete')
            clusters = shc.fcluster(linkresult, max_d, criterion='distance')
            if NVar > 0:
                frac_NA_genotypes_mean = NA_genotypes_mean/NVar
            else:
                frac_NA_genotypes_mean = NA_genotypes_mean
            if not table_hap_groups_pre.empty:
                max_cluster_pre = table_hap_groups_pre['haplotype_group_pre'].max()
                table_hap_groups = pd.DataFrame(data={'labs': list(haplotypes), 'haplotype_group': list(clusters+max_cluster_pre), 'Mean_NAN_haplo': list(frac_NA_genotypes_mean)})
                comparisonDF = table_hap_groups.merge(table_hap_groups_pre, on='labs', how='left').dropna().groupby(['haplotype_group', 'haplotype_group_pre']).agg({'count'})['labs'].reset_index()
                idx = comparisonDF.groupby('haplotype_group')['count'].transform(max) == comparisonDF['count']
                idx2 = comparisonDF[idx].groupby('haplotype_group_pre')['count'].transform(max) == comparisonDF[idx]['count']
                previous_hap_groups = comparisonDF[idx][idx2].groupby('haplotype_group_pre').first().reset_index().groupby('haplotype_group').first().reset_index()
                table_hap_groups_merge = table_hap_groups.merge( previous_hap_groups[['haplotype_group', 'haplotype_group_pre']], on='haplotype_group', how='left')
                table_hap_groups_merge['haplotype_group_ed'] = [int(table_hap_groups_merge['haplotype_group'][x]) if math.isnan(table_hap_groups_merge['haplotype_group_pre'][x]) else int(table_hap_groups_merge['haplotype_group_pre'][x]) for x in range(0, table_hap_groups_merge.shape[0])]
                table_hap_groups_merge['win'] = win_current*win_size
                table_hap_groups_merge[['win', 'labs', 'haplotype_group_ed', 'Mean_NAN_haplo']].to_csv(f"{output_name}.txt", sep="\t", index=False, header=False, mode='a')
                table_hap_groups_pre = table_hap_groups_merge[['labs', 'haplotype_group_ed']].rename(columns={"labs": "labs", "haplotype_group_ed": "haplotype_group_pre"})
                SS_array = SS_var*1
                NA_genotypes_mean = NA_genotypes
                NVar = 1
                win_current = win
            else:
                table_hap_groups = pd.DataFrame(data={'labs': list(haplotypes), 'haplotype_group': list(clusters), 'Mean_NAN_haplo': list(frac_NA_genotypes_mean)})
                table_hap_groups['win'] = win_current*win_size
                table_hap_groups[['win', 'labs', 'haplotype_group', 'Mean_NAN_haplo']].to_csv(f"{output_name}.txt", sep="\t", index=False, header=False, mode='a')
                table_hap_groups_pre = table_hap_groups[['labs', 'haplotype_group']].rename(columns={"labs": "labs", "haplotype_group": "haplotype_group_pre"})
                SS_array = SS_var*1
                NA_genotypes_mean = NA_genotypes
                NVar = 1
                win_current = win



```



bash script: ClusteringHaplotypes.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00
#SBATCH -J clusterHaplotypes

genotype_file=$1
win_size=$2
output_name=$3
maxdif=$4

# python3 /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesPopPar_Win10000_chr01_missing_genotypeTable.txt 10000 Cluster_hapIDs_chr01
python3 ClusteringHaplotypes.py $genotype_file $win_size $output_name $maxdif

```



Job:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/03_Nhaplotypes

mamba activate env_others

#/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr$chN"_"genotypeTable_miltipleAlleles.txt

# using window size of 10 kb and min number of SS 10
for chN in $( seq -w 1 12); do echo $chN
sbatch ClusteringHaplotypes.sh /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr$chN"_"missing_genotypeTable.txt 10000 haplotypeGroups_min10SS_chr$chN 10
done
# 05.04.2024: 2456003..2456014

# using window size of 100 kb and min number of SS 100
for chN in $( seq -w 1 12); do echo $chN
sbatch ClusteringHaplotypes.sh /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr$chN"_"missing_genotypeTable.txt 100000 haplotypeGroups_100kbW_min100SS_chr$chN 100
done
# old: 2105406..2105417
# old: 2150961..2150972

# using window size of 100 kb and min number of SS 20
for chN in $( seq -w 1 12); do echo $chN
sbatch ClusteringHaplotypes.sh /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr$chN"_"missing_genotypeTable.txt 100000 haplotypeGroups_100kbW_min20SS_chr$chN 20
done
# old: 2105418..2105429
# old: 2150973..2150984

#for i in {2223055..2223066}; do stop $i ; rm slurm-$i* ; done

```

### Plots 10 kb:

```R
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

chrID <- "10"

for (chrID in chr_list) {
  
  table_hapG <- read_table(paste0("haplotypeGroups_min10SS_chr",chrID,".txt"), F)
  names(table_hapG) <- c("pos", "sample", "hapGroup", "Mean_NAN_haplo")
  
  table_hapG_sharedVal <- table_hapG %>%
    filter(Mean_NAN_haplo<0.9) %>%
    group_by(pos, hapGroup) %>%
    summarise(NShared=n()) 
  
  
  table_hapG %>% 
    left_join(table_hapG_sharedVal, by=c("pos", "hapGroup")) %>%
    ggplot(aes(pos/1000000, sample)) +
    geom_tile(aes(fill=NShared, colour=NShared))+
    #scico::scale_fill_scico(palette = "lajolla")+
    scale_fill_viridis_c(option = "magma", direction = -1, limits = c(1, 40))+
    scale_colour_viridis_c(option = "magma", direction = -1, limits = c(1, 40))+
    scale_x_continuous(breaks=seq(0,100,5), name = "Position (Mb)", expand = c(0,0))+
    scale_y_discrete(name = "Haplotype", expand = c(0,0))+
    ggtitle(paste0("Chr:", chrID))+
    theme_classic()+
    theme(panel.grid.major = element_line(color = "grey80"),
          panel.grid.minor = element_line(color = "grey80"))
  
  ggsave(paste0("06_haplotypeSharing_per10kbWin_chr", chrID,".png"), height=7, width=14)
  ggsave(paste0("06_haplotypeSharing_per10kbWin_chr", chrID,".pdf"), height=7, width=14)
  
  #
  WinSizeUsed = 100000
  table_hapG %>%
    left_join(table_hapG_sharedVal, by=c("pos", "hapGroup")) %>%
    mutate(win100=floor(pos/WinSizeUsed)) %>%
    group_by(sample, win100) %>%
    summarise(MedianNShared=mean(NShared, na.rm = T)) %>%
    ggplot(aes(win100*WinSizeUsed/1000000, sample, fill=MedianNShared, colour=MedianNShared)) +
    geom_tile()+
    #scico::scale_fill_scico(palette = "lajolla")+
    scale_fill_viridis_c(option = "magma", direction = -1, limits = c(1, 40))+
    scale_colour_viridis_c(option = "magma", direction = -1, limits = c(1, 40))+
    scale_x_continuous(breaks=seq(0,100,5), name = "Position (Mb)", expand = c(0,0))+
    scale_y_discrete(name = "Haplotype", expand = c(0,0))+
    ggtitle(paste0("Chr:", chrID))+
    theme_classic()+
    theme(panel.grid.major = element_line(color = "grey80"),
          panel.grid.minor = element_line(color = "grey80"))
  
  
  ggsave(paste0("06b_100kbmedian_haplotypeSharing_per10kbWin_chr", chrID,".png"), height=12, width=15)
  ggsave(paste0("06b_100kbmedian_haplotypeSharing_per10kbWin_chr", chrID,".pdf"), height=12, width=15)
  # #
  # #
  WinSizeUsed = 1000000
  
  table_hapG %>%
    left_join(table_hapG_sharedVal, by=c("pos", "hapGroup")) %>%
    mutate(win100=floor(pos/WinSizeUsed)) %>%
    group_by(sample, win100) %>%
    summarise(MedianNShared=mean(NShared, na.rm = T)) %>%
    ggplot(aes(win100*WinSizeUsed/1000000, sample, fill=MedianNShared, colour=MedianNShared)) +
    geom_tile()+
    #scico::scale_fill_scico(palette = "lajolla")+
    scale_fill_viridis_c(option = "magma", direction = -1, limits = c(1, 40))+
    scale_colour_viridis_c(option = "magma", direction = -1, limits = c(1, 40))+
    scale_x_continuous(breaks=seq(0,100,5), name = "Position (Mb)", expand = c(0,0))+
    scale_y_discrete(name = "Haplotype", expand = c(0,0))+
    ggtitle(paste0("Chr:", chrID))+
    theme_classic()+
    theme_classic()+
    theme(panel.grid.major = element_line(color = "grey80"),
          panel.grid.minor = element_line(color = "grey80"),
          axis.text = element_text(size=18),
          axis.title = element_text(size=18))
  
  
  # ggsave(paste0("06c_1Mbmedian_haplotypeSharing_per10kbWin_chr", chrID,".png"), height=12, width=15*((max(table_hapG$pos)/1000000)/100))
  ggsave(paste0("06c_1Mbmedian_haplotypeSharing_per10kbWin_chr", chrID,".png"), height=12, width=15)
  ggsave(paste0("06c_1Mbmedian_haplotypeSharing_per10kbWin_chr", chrID,".pdf"), height=12, width=15)
  #
  WinSizeUsed = 1000000
  
  table_hapG %>%
    left_join(table_hapG_sharedVal, by=c("pos", "hapGroup")) %>%
    mutate(win100=floor(pos/WinSizeUsed)) %>%
    group_by(sample, win100) %>%
    summarise(MedianNShared=mean(NShared, na.rm = T)) %>%
    ggplot(aes(win100*WinSizeUsed/1000000, sample)) +
    geom_tile(aes(fill=MedianNShared))+
    #scico::scale_fill_scico(palette = "lajolla")+
    scale_fill_viridis_c(option = "magma", direction = -1, limits = c(1, 40))+
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
  
  
  ggsave(paste0("06d_1Mbmedian_haplotypeSharing_per10kbWin_chr", chrID,".png"), height=12, width=15*((max(table_hapG$pos)/1000000)/100))
  
  table_hapG %>%
    group_by(pos, hapGroup) %>%
    summarise(NShared=n()) %>%
    group_by(pos) %>%
    summarise(MaxNShared=max(NShared, na.rm = T)) %>%
    ggplot(aes(pos/1000000, MaxNShared))+
    geom_line(alpha=0.05)+
    geom_point(alpha=0.1)+
    scale_x_continuous(breaks=seq(0,100,5), name = "Position (Mb)", expand = c(0,0))+
    ggtitle(paste0("Chr:", chrID))+
    theme_classic()+
    theme(panel.grid.major = element_line(color = "grey80"),
          panel.grid.minor = element_line(color = "grey80"))
  
  ggsave(paste0("07_MaxhaplotypeSharing_per10kbWin_chr", chrID,".png"), height=5, width=12)
  
  WinSizeU=500000
  if (chrID == "01"){
    merged_table <- table_hapG %>% 
      group_by(pos, hapGroup) %>%
      summarise(NShared=n(),
                Mean_NAN_haplo = mean(Mean_NAN_haplo)) %>%
      mutate(NAN_haplo=if_else(Mean_NAN_haplo>0.9, 1,0)) %>%
      # filter(NAN_haplo==0) %>%
      group_by(pos) %>%
      summarise(MaxNShared=max(NShared*(NAN_haplo==0)),
                NHaplotypes = length(unique(hapGroup*(NAN_haplo==0)))-1) %>%
      mutate(win_ed=floor(pos/WinSizeU)) %>%
      # filter(win_ed==21500000/WinSizeU) %>%
      group_by(win_ed) %>%
      summarise(MeanMaxNShared=mean(MaxNShared),
                MeanNHaplotypes=mean(NHaplotypes)) %>%
      mutate(chr=paste0("chr", chrID),
             pos=win_ed*WinSizeU) %>%
      ungroup() %>% 
      select(chr, pos, MeanMaxNShared, MeanNHaplotypes) 
    merged_tableNA <- table_hapG %>% 
      group_by(pos, hapGroup) %>%
      summarise(NShared=n(),
                Mean_NAN_haplo = mean(Mean_NAN_haplo)) %>%
      mutate(NAN_haplo=if_else(Mean_NAN_haplo>0.9, 1,0)) %>%
      # filter(NAN_haplo==0) %>%
      group_by(pos) %>%
      summarise(TotalNASamples=sum(NShared*(NAN_haplo==1))) %>%
      mutate(win_ed=floor(pos/WinSizeU)) %>%
      group_by(win_ed) %>%
      summarise(MeanTotalNASamples=mean(TotalNASamples)) %>%
      mutate(chr=paste0("chr", chrID),
             pos=win_ed*WinSizeU) %>%
      ungroup() %>% 
      select(chr, pos, MeanTotalNASamples)
  } else{
    merged_table <- rbind(merged_table,
                          table_hapG %>% 
                            group_by(pos, hapGroup) %>%
                            summarise(NShared=n(),
                                      Mean_NAN_haplo = mean(Mean_NAN_haplo)) %>%
                            mutate(NAN_haplo=if_else(Mean_NAN_haplo>0.9, 1,0)) %>%
                            # filter(NAN_haplo==0) %>%
                            group_by(pos) %>%
                            summarise(MaxNShared=max(NShared*(NAN_haplo==0)),
                                      NHaplotypes = length(unique(hapGroup*(NAN_haplo==0)))-1) %>%
                            mutate(win_ed=floor(pos/WinSizeU)) %>%
                            # filter(win_ed==21500000/WinSizeU) %>%
                            group_by(win_ed) %>%
                            summarise(MeanMaxNShared=mean(MaxNShared),
                                      MeanNHaplotypes=mean(NHaplotypes)) %>%
                            mutate(chr=paste0("chr", chrID),
                                   pos=win_ed*WinSizeU) %>%
                            ungroup() %>% 
                            select(chr, pos, MeanMaxNShared, MeanNHaplotypes))
    merged_tableNA <- rbind(merged_tableNA, table_hapG %>% 
                              group_by(pos, hapGroup) %>%
                              summarise(NShared=n(),
                                        Mean_NAN_haplo = mean(Mean_NAN_haplo)) %>%
                              mutate(NAN_haplo=if_else(Mean_NAN_haplo>0.9, 1,0)) %>%
                              # filter(NAN_haplo==0) %>%
                              group_by(pos) %>%
                              summarise(TotalNASamples=sum(NShared*(NAN_haplo==1))) %>%
                              mutate(win_ed=floor(pos/WinSizeU)) %>%
                              group_by(win_ed) %>%
                              summarise(MeanTotalNASamples=mean(TotalNASamples)) %>%
                              mutate(chr=paste0("chr", chrID),
                                     pos=win_ed*WinSizeU) %>%
                              ungroup() %>% 
                              select(chr, pos, MeanTotalNASamples))
  }
  
  
  table_hapG %>%
    group_by(pos, hapGroup) %>%
    summarise(NShared=n(),
              Mean_NAN_haplo = mean(Mean_NAN_haplo)) %>%
    mutate(NAN_haplo=if_else(Mean_NAN_haplo>0.9, "NA", "Hap")) %>%
    group_by(pos, NAN_haplo) %>%
    summarise(MaxNShared=max(NShared),
              NHaplotypes = length(unique(hapGroup))) %>%
    mutate(win100=floor(pos/100000)) %>%
    group_by(win100, NAN_haplo) %>%
    summarise(MeanMaxNShared=mean(MaxNShared),
              MeanNHaplotypes=mean(NHaplotypes)) %>%
    ggplot(aes(win100/10, MeanMaxNShared))+
    geom_line(alpha=0.5)+
    geom_point(alpha=0.1)+
    geom_line(aes(win100/10, MeanNHaplotypes), colour="darkgreen", alpha=0.5)+
    geom_point(aes(win100/10, MeanNHaplotypes), colour="darkgreen", alpha=0.1)+
    scale_x_continuous(breaks=seq(0,100,5), name = "Position (Mb)", expand = c(0,0))+
    ggtitle(paste0("Chr:", chrID))+
    facet_grid(NAN_haplo ~ .)+
    theme_classic()+
    theme(panel.grid.major = element_line(color = "grey80"),
          panel.grid.minor = element_line(color = "grey80"))
  
  
  ggsave(paste0("08_MeanMaxhaplotypeSharing_per100kbWin_chr", chrID,".png"), height=5, width=12)
  
  
  table_hapG %>% 
    ggplot(aes(pos/1000000, sample)) +
    geom_tile(aes(fill=Mean_NAN_haplo, colour=Mean_NAN_haplo))+
    #scico::scale_fill_scico(palette = "lajolla")+
    scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, 1))+
    scale_colour_viridis_c(option = "magma", direction = -1, limits = c(0, 1))+
    scale_x_continuous(breaks=seq(0,100,5), name = "Position (Mb)", expand = c(0,0))+
    scale_y_discrete(name = "Haplotype", expand = c(0,0))+
    ggtitle(paste0("Chr:", chrID))+
    theme_classic()+
    theme(panel.grid.major = element_line(color = "grey80"),
          panel.grid.minor = element_line(color = "grey80"))
  
  ggsave(paste0("09_haplotype_Missing_per10kbWin_chr", chrID,".png"), height=7, width=14)
  ggsave(paste0("09_haplotype_Missing_per10kbWin_chr", chrID,".pdf"), height=7, width=14)
  
}




# merged_table %>%
#   write.table(file=paste0("./merged_haplotypeGroups_min10SS_allChr_500KbWin.csv"), sep = ",", row.names = F, col.names = T)
# 
# merged_tableNA %>%
#   write.table(file=paste0("./merged_tableNA_haplotypeGroups_min10SS_allChr_500KbWin.csv"), sep = ",", row.names = F, col.names = T)

merged_table <- read.csv(paste0("./merged_haplotypeGroups_min10SS_allChr_500KbWin.csv"), T)
merged_tableNA <- read.csv(paste0("./merged_tableNA_haplotypeGroups_min10SS_allChr_500KbWin.csv"), T)


MeanNHaplotypes_GW_noNA <- merged_table %>%
  ungroup() %>%
  pull(MeanNHaplotypes) %>%
  mean()


MeanMaxNShared_GW_noNA <- merged_table %>%
  ungroup() %>%
  pull(MeanMaxNShared) %>%
  mean()


MeanMaxNShared_GW_NA <- merged_tableNA %>%
  ungroup() %>%
  pull(MeanTotalNASamples) %>%
  mean()


head(merged_tableNA)

# merged_table_MeanMaxNShared <- merged_table %>% 
#   mutate(NAN_haplo=if_else(NAN_haplo == 1, "nan", "Hap")) %>%
#   select(-MeanNHaplotypes) %>% 
#   spread(NAN_haplo, MeanMaxNShared, fill = 0) %>% 
#   mutate(total_NotNa=40-nan) %>%
#   gather(key="NAN_haplo", value = "MeanMaxNShared", Hap, nan, total_NotNa) 

merged_table_MeanMaxNShared

quantile_5_fq <- merged_table %>%
  left_join(merged_tableNA %>%
              mutate(TotalMeanNHaplotypes=40-MeanTotalNASamples), by=c("chr", "pos")) %>%
  mutate(prop_haplotypes=MeanNHaplotypes/TotalMeanNHaplotypes) %>%
  filter(MeanNHaplotypes<24) %>%
  pull(prop_haplotypes) %>%
  quantile(0.055)
# head()
# ggplot(aes(prop_haplotypes)) +
# geom_histogram()
# head()


# merged_table %>% 
#   left_join(merged_tableNA %>%
#               mutate(TotalMeanNHaplotypes=40-MeanTotalNASamples), by=c("chr", "pos")) %>%
#   mutate(prop_haplotypes=MeanNHaplotypes/TotalMeanNHaplotypes) %>%
#   filter(prop_haplotypes<24) %>%
#   pull(prop_haplotypes) %>%
#   quantile(0.055)
# 
quantile_5_count <-merged_table  %>%
  filter(MeanNHaplotypes<24) %>%
  pull(MeanNHaplotypes) %>%
  quantile(0.055)




merged_table %>% 
  filter(MeanNHaplotypes<24) %>%
  ggplot(aes(pos/1000000, MeanNHaplotypes))+
  geom_ribbon(aes(ymin=ifelse(MeanNHaplotypes<MeanNHaplotypes_GW_noNA, MeanNHaplotypes, MeanNHaplotypes_GW_noNA), 
                  ymax=MeanNHaplotypes_GW_noNA), fill="red") +
  geom_point(data=merged_table %>% 
               left_join(merged_tableNA %>%
                           mutate(TotalMeanNHaplotypes=40-MeanTotalNASamples), by=c("chr", "pos")) %>%
               mutate(prop_haplotypes=MeanNHaplotypes/TotalMeanNHaplotypes) %>%
               filter(MeanNHaplotypes<24) %>%
               filter(prop_haplotypes<quantile_5_fq) %>%
               # filter(MeanNHaplotypes<quantile_5_count) %>%
               mutate(posy=(1)), aes(pos/1000000, posy), colour="darkred") +
  geom_line(alpha=0.8)+
  geom_line(data=centromere_table %>%
              gather(key="posB", value = "pos", -chr) %>%
              mutate(posy = (-1)), 
            aes(pos/1000000,posy), colour="darkgreen", size=1.5) +
  geom_hline(yintercept=MeanNHaplotypes_GW_noNA, linetype="dashed", color = "red")+
  geom_line(data=merged_tableNA %>%
              mutate(MeanNHaplotypes=40-MeanTotalNASamples), 
            aes(pos/1000000,MeanNHaplotypes), colour="grey30", alpha=0.5) +
  scale_y_continuous(breaks=seq(0,50,5), limits = c((-1), 40)) +
  scale_x_continuous(breaks=seq(0,100,15)) +
  labs(x = "Pos. (Mb)",
       y = "Mean N. Haplo.") +
  facet_grid(. ~ chr, scale="free", space="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12), 
        axis.title = element_text(size=14))


ggsave("10_ed_MeanNHaplot_per500kbWin.png", height=2.5, width=15)
ggsave("10_ed_MeanNHaplot_per500kbWin.pdf", height=2.5, width=15)


merged_table %>% 
  # filter(MeanNHaplotypes<24) %>%
  pull(MeanNHaplotypes) %>%
  # pull(MeanMaxNShared) %>%
  summary()

merged_table  %>%
  # filter(MeanNHaplotypes<24) %>%
  pull(MeanNHaplotypes) %>%
  # quantile(0.95) 
  quantile(0.05)

merged_table %>% 
  # filter(chr=="chr01") %>%
  # filter(pos> 38000000) %>%
  # filter(pos< 41000000) 
  filter(chr=="chr10") %>%
  filter(pos> 7000000) %>%
  filter(pos<40000000) 
# filter(MeanNHaplotypes<24) %>%
# pull(MeanMaxNShared) %>%
# summary()


quantile_95 <- merged_table %>% 
  pull(MeanMaxNShared) %>%
  quantile(0.94)

quantile_95



quantile_94_fq <-merged_table %>%
  left_join(merged_tableNA %>%
              mutate(TotalMeanNHaplotypes=40-MeanTotalNASamples), by=c("chr", "pos")) %>%
  mutate(prop_haplotypes=(MeanMaxNShared/TotalMeanNHaplotypes)) %>%
  pull(prop_haplotypes) %>%
  quantile(0.9)


quantile_94_count <- merged_table %>% 
  pull(MeanMaxNShared) %>%
  quantile(0.94)

quantile_94_fq
quantile_94_count

merged_table %>% 
  ggplot(aes(pos/1000000, MeanMaxNShared))+
  geom_ribbon(aes(ymax=ifelse(MeanMaxNShared>MeanMaxNShared_GW_noNA, MeanMaxNShared, MeanMaxNShared_GW_noNA), 
                  ymin=MeanMaxNShared_GW_noNA), fill="#7f97ca") +
  geom_point(data=merged_table %>%
               left_join(merged_tableNA %>%
                           mutate(TotalMeanNHaplotypes=40-MeanTotalNASamples), by=c("chr", "pos")) %>%
               mutate(prop_haplotypes=(MeanMaxNShared/TotalMeanNHaplotypes)) %>%
               filter(prop_haplotypes>quantile_94_fq) %>%
               filter(MeanMaxNShared>quantile_94_count) %>%
               # merged_table %>% 
               # filter(MeanMaxNShared>quantile_95) %>%
               mutate(posy=(5)), aes(pos/1000000, posy), colour="#2859c2") +
  geom_line(alpha=0.8)+
  geom_line(data=centromere_table %>%
              gather(key="posB", value = "pos", -chr) %>%
              mutate(posy = 3), 
            aes(pos/1000000,posy), colour="darkgreen", size=1.5) +
  geom_hline(yintercept=MeanMaxNShared_GW_noNA, linetype="dashed", color = "red")+
  geom_line(data=merged_tableNA %>%
              mutate(MeanNHaplotypes=40-MeanTotalNASamples), 
            aes(pos/1000000,MeanNHaplotypes), colour="grey30", alpha=0.5) +
  scale_y_continuous(breaks=seq(0,50,5), limits = c(0, 40)) +
  scale_x_continuous(breaks=seq(0,100,15)) +
  labs(x = "Pos. (Mb)",
       y = "Max N. Shared Hap.") +
  facet_grid(. ~ chr, scale="free", space="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12), 
        axis.title = element_text(size=14))


ggsave("11_ed_MeanMaxNShared_per500kbWin.png", height=2.5, width=15)
ggsave("11_ed_MeanMaxNShared_per500kbWin.pdf", height=2.5, width=15)


```



### Plots 100 kb:

```R
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

chr <- "01"

for (chr in chr_list) {
table_hapG <- read_table(paste0("haplotypeGroups_100kbW_min100SS_chr",chr,".txt"), F)
names(table_hapG) <- c("pos", "sample", "hapGroup")

head(table_hapG)

# calculate number of shared haplotypes per window
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
  ggtitle(paste0("Chr:", chr))+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))


ggsave(paste0("06_haplotypeSharing_per100kbWin_min100SS_chr", chr,".png"), height=7, width=7)


table_hapG %>%
  group_by(pos, hapGroup) %>%
  summarise(NShared=n()) %>%
  group_by(pos) %>%
  summarise(MaxNShared=max(NShared)) %>%
  ggplot(aes(pos/1000000, MaxNShared))+
  geom_line(alpha=0.5)+
  geom_point(alpha=0.2)+
  scale_x_continuous(breaks=seq(0,100,5), name = "Position (Mb)", expand = c(0,0))+
  scale_y_continuous(limits = c(0, 44))+
  ggtitle(paste0("Chr:", chr))+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))

ggsave(paste0("07_MaxhaplotypeSharing_per100kbWin_min100SS_chr", chr,".png"), height=5, width=12)


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
   ggplot(aes(pos/1000000, MaxNShared))+
  geom_line(alpha=0.5)+
  geom_point(alpha=0.1)+
  geom_line(aes(pos/1000000, SecondMaxNShared), colour="red", alpha=0.5)+
  geom_point(aes(pos/1000000, SecondMaxNShared), colour="red", alpha=0.1)+
  scale_x_continuous(breaks=seq(0,100,5), name = "Position (Mb)", expand = c(0,0))+
  scale_y_continuous(limits = c(0, 44))+
  ggtitle(paste0("Chr:", chr))+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))
  

ggsave(paste0("08_MeanMaxhaplotypeSharing_top2_per1000kbWin_min100SS_chr", chr,".png"), height=5, width=12)

}





# plots using min number of SS of 22 sites
for (chr in chr_list) {
table_hapG <- read_table(paste0("haplotypeGroups_100kbW_min20SS_chr",chr,".txt"), F)
names(table_hapG) <- c("pos", "sample", "hapGroup")

head(table_hapG)

# calculate number of shared haplotypes per window
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
  ggtitle(paste0("Chr:", chr))+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))


ggsave(paste0("06_haplotypeSharing_per100kbWin_min20SS_chr", chr,".png"), height=7, width=7)


table_hapG %>%
  group_by(pos, hapGroup) %>%
  summarise(NShared=n()) %>%
  group_by(pos) %>%
  summarise(MaxNShared=max(NShared)) %>%
  ggplot(aes(pos/1000000, MaxNShared))+
  geom_line(alpha=0.5)+
  geom_point(alpha=0.2)+
  scale_x_continuous(breaks=seq(0,100,5), name = "Position (Mb)", expand = c(0,0))+
  scale_y_continuous(limits = c(0, 44))+
  ggtitle(paste0("Chr:", chr))+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))

ggsave(paste0("07_MaxhaplotypeSharing_per100kbWin_min20SS_chr", chr,".png"), height=5, width=12)


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
   ggplot(aes(pos/1000000, MaxNShared))+
  geom_line(alpha=0.5)+
  geom_point(alpha=0.1)+
  geom_line(aes(pos/1000000, SecondMaxNShared), colour="red", alpha=0.5)+
  geom_point(aes(pos/1000000, SecondMaxNShared), colour="red", alpha=0.1)+
  scale_x_continuous(breaks=seq(0,100,5), name = "Position (Mb)", expand = c(0,0))+
  scale_y_continuous(limits = c(0, 44))+
  ggtitle(paste0("Chr:", chr))+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))
  

ggsave(paste0("08_MeanMaxhaplotypeSharing_top2_per1000kbWin_min20SS_chr", chr,".png"), height=5, width=12)

}

```



## Clustering using overlapping Windows and remove missing data:

Version using overlaping windows:

script: ClusteringHaplotypes_overlappingWin.py

```python
import os
import sys
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as shc
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
import math

# genotype_file = str("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr01_missing_genotypeTable.txt")
# win_size = int(10000)
# output_name = str("test")
# max_d = int(10)

# this script takes a genotype table, a window size and and output name and: produce a distance matrix and cluster haplotypes with a maximum of around 10SNPs (max_d parameter).
# the output is a table with window start position, sample and haplotype ID after clustering.
# The scripts uses the haplotype groups produced from the previous window to allocate haplotypes ID in the current window. 

# genotype_file = str("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr01_missing_genotypeTable.txt")
# win_size = int(1000000)
# step_size= int(50000)
# max_d = int(2000)
# output_name = str("test")


genotype_file = str(sys.argv[1])
win_size = int(sys.argv[2])
step_size=int(sys.argv[3])
max_d = int(sys.argv[4])
output_name = str(sys.argv[5])

def str_to_int(x):
    if x!='nan':
        return int(x)
    else:
        return int(-1)


dic_win_SS_array = {}
dic_win_NA_genotypes_mean = {}
dic_win_NVar = {}

for line in open(genotype_file, "r"):
    if line[0:3] == "CHR":
        haplotypes = line.strip().split("\t")[5:]
        list_win_pre = [0] 
        dic_win_SS_array[0] = np.zeros((len(haplotypes), len(haplotypes)))
        dic_win_NA_genotypes_mean[0] = np.zeros((len(haplotypes)))
        dic_win_NVar[0] = 0
        table_hap_groups_pre = pd.DataFrame()
        if os.path.exists(f"{output_name}.txt"):
            os.remove(f"{output_name}.txt")
    else:
        line_values = line.strip().split("\t")
        pos = int(line_values[1])
        list_win = [i for i in range(int((pos-win_size)/step_size),int(pos/step_size)+1) if i>=0]
        list_win_withPre = set(list_win_pre + list_win)
        genotypes = np.asarray([str_to_int(x) for x in line_values[5:]])
        SS_var = genotypes[:,None] != genotypes
        NA_genotypes = (genotypes==(-1))*1
        NA_genotypes_boolean = (genotypes==(-1))
        SS_var[NA_genotypes_boolean, :] = False
        SS_var[:, NA_genotypes_boolean] = False
        #SS_var[NA_genotypes_boolean, :][:, NA_genotypes_boolean] = False
        #SS_var[NA_genotypes_boolean,] = False
        #SS_var[:,NA_genotypes_boolean] = False
        for win_current in list_win_withPre:
            if not win_current in dic_win_SS_array:
                dic_win_SS_array[win_current] = np.zeros((len(haplotypes), len(haplotypes)))
                dic_win_NA_genotypes_mean[win_current] = np.zeros((len(haplotypes)))
                dic_win_NVar[win_current] = 0
            if win_current in list_win:
                dic_win_SS_array[win_current] += SS_var
                dic_win_NA_genotypes_mean[win_current] += NA_genotypes
                dic_win_NVar[win_current] += 1
            else:
                condensed_dist = squareform(dic_win_SS_array[win_current])
                linkresult = shc.linkage(condensed_dist, method='complete')
                clusters = shc.fcluster(linkresult, max_d, criterion='distance')
                #Z = shc.ward(pdist(dic_win_SS_array[win_current]))
                #clusters = shc.fcluster(Z, max_d, criterion='distance')
                if dic_win_NVar[win_current] > 0:
                    frac_NA_genotypes_mean = dic_win_NA_genotypes_mean[win_current]/dic_win_NVar[win_current]
                else:
                    frac_NA_genotypes_mean = dic_win_NA_genotypes_mean[win_current]
                if not table_hap_groups_pre.empty:
                    max_cluster_pre = table_hap_groups_pre['haplotype_group_pre'].max()
                    table_hap_groups = pd.DataFrame(data={'labs': list(haplotypes), 'haplotype_group': list(clusters+max_cluster_pre), 'Mean_NAN_haplo': list(frac_NA_genotypes_mean)})
                    comparisonDF = table_hap_groups.merge(table_hap_groups_pre, on='labs', how='left').dropna().groupby(['haplotype_group', 'haplotype_group_pre']).agg({'count'})['labs'].reset_index()
                    idx = comparisonDF.groupby('haplotype_group')['count'].transform(max) == comparisonDF['count']
                    idx2 = comparisonDF[idx].groupby('haplotype_group_pre')['count'].transform(max) == comparisonDF[idx]['count']
                    previous_hap_groups = comparisonDF[idx][idx2].groupby('haplotype_group_pre').first().reset_index().groupby('haplotype_group').first().reset_index()
                    table_hap_groups_merge = table_hap_groups.merge( previous_hap_groups[['haplotype_group', 'haplotype_group_pre']], on='haplotype_group', how='left')
                    table_hap_groups_merge['haplotype_group_ed'] = [int(table_hap_groups_merge['haplotype_group'][x]) if math.isnan(table_hap_groups_merge['haplotype_group_pre'][x]) else int(table_hap_groups_merge['haplotype_group_pre'][x]) for x in range(0, table_hap_groups_merge.shape[0])]
                    table_hap_groups_merge['win'] = win_current*step_size
                    table_hap_groups_merge[['win', 'labs', 'haplotype_group_ed', 'Mean_NAN_haplo']].to_csv(f"{output_name}.txt", sep="\t", index=False, header=False, mode='a')
                    table_hap_groups_pre = table_hap_groups_merge[['labs', 'haplotype_group_ed']].rename(columns={"labs": "labs", "haplotype_group_ed": "haplotype_group_pre"})
                    del dic_win_SS_array[win_current]
                    del dic_win_NA_genotypes_mean[win_current]
                    del dic_win_NVar[win_current] 
                else:
                    table_hap_groups = pd.DataFrame(data={'labs': list(haplotypes), 'haplotype_group': list(clusters), 'Mean_NAN_haplo': list(frac_NA_genotypes_mean)})
                    table_hap_groups['win'] = win_current*step_size
                    table_hap_groups[['win', 'labs', 'haplotype_group', 'Mean_NAN_haplo']].to_csv(f"{output_name}.txt", sep="\t", index=False, header=False, mode='a')
                    table_hap_groups_pre = table_hap_groups[['labs', 'haplotype_group']].rename(columns={"labs": "labs", "haplotype_group": "haplotype_group_pre"})
                    del dic_win_SS_array[win_current]
                    del dic_win_NA_genotypes_mean[win_current]
                    del dic_win_NVar[win_current]
        list_win_pre = list_win


```



bash script: ClusteringHaplotypes_overlappingWin.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00
#SBATCH -J clusterHaplotypes

genotype_file=$1
win_size=$2
step_size=$3
maxdif=$4
output_name=$5

# python3 /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesPopPar_Win10000_chr01_missing_genotypeTable.txt 1000000 50000 2000 Cluster_hapIDs_chr01
python3 ClusteringHaplotypes_overlappingWin.py $genotype_file $win_size $step_size $maxdif $output_name 

```



Job:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/03_Nhaplotypes

mamba activate mypython3

#/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr$chN"_"genotypeTable_miltipleAlleles.txt

# using window size of 1Mb and min number of SS 1000. Step size 500000
for chN in $( seq -w 1 12); do echo $chN
sbatch ClusteringHaplotypes_overlappingWin.sh /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr$chN"_"missing_genotypeTable.txt 1000000 500000 1000 haplotypeGroups_min1000SS_1MbWin_chr$chN 
done
# 05.04.2024:2456015..2456026

# using window size of 10kb and min number of SS 10. Step size 5000
for chN in $( seq -w 1 12); do echo $chN
sbatch ClusteringHaplotypes_overlappingWin.sh /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr$chN"_"missing_genotypeTable.txt 10000 5000 10 haplotypeGroups_min10SS_10kbWin_chr$chN 
done
# 05.04.2024:2456027..2456038

# using window size of 50kb and min number of SS 50. Step size 25000
for chN in $( seq -w 1 12); do echo $chN
sbatch ClusteringHaplotypes_overlappingWin.sh /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr$chN"_"missing_genotypeTable.txt 50000 25000 50 haplotypeGroups_min50SS_50kbWin_chr$chN 
done
# 05.04.2024:2456039..2456050

# using window size of 100kb and min number of SS 100. Step size 50000
for chN in $( seq -w 1 12); do echo $chN
sbatch ClusteringHaplotypes_overlappingWin.sh /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr$chN"_"missing_genotypeTable.txt 100000 50000 100 haplotypeGroups_min100SS_100kbWin_chr$chN 
done
# 05.04.2024:2456039..2456050

# using window size of 200kb and min number of SS 200. Step size 100000
for chN in $( seq -w 1 12); do echo $chN
sbatch ClusteringHaplotypes_overlappingWin.sh /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr$chN"_"missing_genotypeTable.txt 200000 100000 200 haplotypeGroups_min200SS_200kbWin_chr$chN 
done
# 05.04.2024:2456051..2456062


# using window size of 500kb and min number of SS 500. Step size 250000
for chN in $( seq -w 1 12); do echo $chN
sbatch ClusteringHaplotypes_overlappingWin.sh /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr$chN"_"missing_genotypeTable.txt 500000 250000 500 haplotypeGroups_min500SS_500kbWin_chr$chN 
done
# 05.04.2024:2456063..2456074


# using window size of 5kb and min number of SS 5. Step size 2500
for chN in $( seq -w 1 12); do echo $chN
sbatch ClusteringHaplotypes_overlappingWin.sh /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr$chN"_"missing_genotypeTable.txt 5000 2500 5 haplotypeGroups_min5SS_5kbWin_chr$chN 
done
# 05.04.2024:2456087..2456098


# using window size of 2kb and min number of SS 2. Step size 1000
for chN in $( seq -w 1 12); do echo $chN
sbatch ClusteringHaplotypes_overlappingWin.sh /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr$chN"_"missing_genotypeTable.txt 2000 1000 2 haplotypeGroups_min2SS_2kbWin_chr$chN 
done
# 05.04.2024:2456099..2456110



# merging all tables to produce number of haplotypes per window:

grep "" haplotypeGroups_min*SS_*Win_chr*.txt | awk '{if($4<0.9) print $1"\t"$3}' | sed 's/haplotypeGroups_min//g' | sed 's/SS_/\t/g' | sed 's/Win_/\t/g' | sed 's/.txt:/\t/g' | sort | uniq | awk '{print $1"\t"$2"\t"$3"\t"$4}' | sort | uniq -c | sed 's/^[ \t]*//'  | awk -F' ' '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$1}' > total_haplotypeGroups_allWinSizes.txt
# done

#for i in {2335499..2335510}; do stop $i ; rm slurm-$i* ; done 

```



### Plots:

```R
# .libPaths("/dss/dsshome1/lxc03/di36guz2/R")
# mamba activate env_others
.libPaths()

setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/03_Nhaplotypes")

library(ggplot2)
library(tidyverse)

centromere_table <- read.csv("centromeres.csv", header = T)

chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

chrID <- "01"

# for (chrID in chr_list) {

table_hapG <- read_table("total_haplotypeGroups_allWinSizes.txt", F)
names(table_hapG) <- c("minSS", "WinSize", "chr", "pos", "NHap")



table_hapG %>% 
  mutate(WinSize=factor(WinSize, levels=c("2kb", "5kb", "10kb", "50kb", "100kb", "200kb", "500kb", "1Mb"))) %>%
  ggplot(aes(pos/1000000, NHap, colour=WinSize, alpha=WinSize)) +
  geom_point()+
  # geom_line(alpha=0.2)+
  # scale_fill_viridis_c(option = "magma", direction = -1, limits = c(1, 40))+
  # scale_colour_viridis_c(option = "magma", direction = -1)+
  scale_x_continuous(breaks=seq(0,100,10), name = "Position (Mb)", expand = c(0,0))+
  # scale_y_continuous(name = "N. Haplotypes", expand = c(0,0))+
  scale_alpha_manual(values=c(0.01, 0.01, 0.01, 0.1, 0.2, 0.3, 0.1))+
  facet_grid(WinSize ~ chr, scales = "free_x", space = "free_x")+
  # ggtitle(paste0("Chr:", chrID))+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"), 
        legend.position = "none")


ggsave(paste0("12_NumberHaplotypes_diffWinSizes.png"), height=10, width=14)
ggsave(paste0("12_NumberHaplotypes_diffWinSizes.pdf"), height=10, width=14)




table_hapG %>% 
  mutate(WinSize=factor(WinSize, levels=c("2kb", "5kb", "10kb", "100kb", "200kb", "500kb", "1Mb"))) %>%
  group_by(WinSize) %>%
  summarise(meanNHap=median(NHap), 
            minNHap=min(NHap),
            maxNHap=max(NHap)) 


WinSizeU=500000

table_hapG %>% 
  mutate(WinSize=factor(WinSize, levels=c("2kb", "5kb", "10kb", "50kb", "100kb", "200kb", "500kb", "1Mb"))) %>%
  mutate(win_ed=floor(pos/WinSizeU)) %>%
  group_by(WinSize, chr, win_ed) %>%
  summarise(meanNHap=mean(NHap), 
            minNHap=min(NHap),
            maxNHap=max(NHap)) %>%
  group_by(WinSize) %>%
  summarise(meanNHapLargeW=mean(meanNHap), 
            minNHapLargeW=min(meanNHap), 
            maxNHapLargeW=max(meanNHap))



table_hapG %>% 
  filter(WinSize=="10kb") %>%
  mutate(win_ed=floor(pos/WinSizeU)) %>%
  group_by(WinSize, chr, win_ed) %>%
  summarise(MeanNHaplotypes=mean(NHap)) %>%
  pull(MeanNHaplotypes) %>%
  mean()

WinSizeU=100000
table_hapG %>% 
  mutate(WinSize=factor(WinSize, levels=c("2kb", "5kb", "10kb", "50kb", "100kb", "200kb", "500kb", "1Mb"))) %>%
  mutate(win_ed=floor(pos/WinSizeU)) %>%
  group_by(WinSize, chr, win_ed) %>%
  summarise(meanNHap=mean(NHap), 
            minNHap=min(NHap),
            maxNHap=max(NHap)) %>%
  ggplot(aes(meanNHap))+
  geom_histogram(binwidth = 1, fill="black") +
  #scale_x_continuous(breaks=seq(0,1,0.02)) +
  #labs(x = "Pairwise divergence SS/SNP",y = "Count") +
  facet_grid(WinSize ~ . , scales = "free")+
  labs(x = "NUmber of Haplotypes",
       y = "Count") +
  theme_classic()


table_hapG %>% 
  mutate(WinSize=factor(WinSize, levels=c("2kb", "5kb", "10kb", "50kb", "100kb", "200kb", "500kb", "1Mb"))) %>%
  ggplot(aes(NHap))+
  geom_histogram(binwidth = 1, fill="#993300") +
  #scale_x_continuous(breaks=seq(0,1,0.02)) +
  #labs(x = "Pairwise divergence SS/SNP",y = "Count") +
  facet_grid(WinSize ~ . , scales = "free")+
  labs(x = "Number of Haplotypes",
       y = "Count") +
  theme_classic()

ggsave("13_distributionNumberHaplotypes_diffWinSizes.png", height=8, width=3)
ggsave("13_distributionNumberHaplotypes_diffWinSizes.pdf", height=8, width=3)






table_hapG %>% 
  filter(WinSize=="10kb") %>%
  pull(NHap) %>%
  quantile(0.95)
head()

mean(table_pop$mean_pi_pw)

table_hapG %>% 
  filter(WinSize=="10kb") %>%
  ggplot(aes(NHap))+
  geom_histogram(binwidth = 1, fill="#AA2E25", colour="#AA2E25") +
  geom_vline(xintercept = table_hapG %>% 
               filter(WinSize=="10kb") %>%
               pull(NHap) %>%
               mean(), color = "red") +
  geom_vline(xintercept = table_hapG %>% 
               filter(WinSize=="10kb") %>%
               pull(NHap) %>%
               quantile(0.95), color = "#084276") +
  scale_x_continuous(breaks=seq(0,45,5)) +
  scale_y_continuous(breaks=seq(0,20000,5000), labels = seq(0,20,5)) +
  labs(x = "Number of Haplotypes",
       y = "Window Count") +
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.title = element_text(size=16), 
        axis.text = element_text(size=14))


ggsave("14_v2_distributionNumberHaplotypes_10KbWin.png", height=2.6, width=3.9)
ggsave("14_v2_distributionNumberHaplotypes_10KbWin.pdf", height=2.6, width=3.9)

```



# Proportion shared regions:





```R

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
  names(table_hapG) <- c("pos", "sample", "hapGroup", "Mean_NAN_haplo")
  
  table_hapG_spread <-  table_hapG %>% 
    select(-Mean_NAN_haplo) %>%
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






chrN="01"
for (chrN in chr_list) {
  lengthHaplBlocks_all <- read.csv(paste0("Table_SharedBlocks_chr",chrN,".csv"), F)
  names(lengthHaplBlocks_all) <- c("chr", "haplotype1", "haplotype2", "startpos", "endpos")
  if (chrN=="01"){
    allChr_lengthHaplBlocks_all <- lengthHaplBlocks_all
  } else{
    allChr_lengthHaplBlocks_all <- allChr_lengthHaplBlocks_all %>%
      rbind(lengthHaplBlocks_all)
  }
}

# raw distribution on block lengths:
allChr_lengthHaplBlocks_all <- allChr_lengthHaplBlocks_all %>%
  mutate(lenHaplBlocks=endpos-startpos+10000)


allChr_lengthHaplBlocks_all %>%
  arrange(-lenHaplBlocks) %>%
  head(30)

allChr_lengthHaplBlocks_all %>%
  filter(lenHaplBlocks>10000) %>%
  # ggplot(aes(log10(lenHaplBlocks)))+
  ggplot(aes(lenHaplBlocks/1000000))+
  # geom_density(fill="black")+
  geom_histogram(fill="black")+
  geom_vline(xintercept = mean(allChr_lengthHaplBlocks_all$lenHaplBlocks/1000000), 
             #linetype="dotted", 
             color = "red") +
  # scale_x_continuous(breaks=seq(0,1000,200))+
  scale_x_continuous(trans="log10") +
  labs(x = "Length Shared Blocks (Mb)",
       y = "Count") +
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor = element_line(color = "grey80"))
# axis.text = element_text(size=24), 
# axis.title = element_text(size=24))

ggsave(paste0("09_a_RawShareHaplLenBlocks_ed_10kbWin_min10SS.png"), height=5, width=7)
ggsave(paste0("09_a_RawShareHaplLenBlocks_ed_10kbWin_min10SS.pdf"), height=5, width=7)


meanChr_allChr_lengthHaplBlocks_all <- allChr_lengthHaplBlocks_all %>%
  group_by(chr) %>%
  summarise(mean_lenHaplBlocks=mean(lenHaplBlocks), 
            max_lenHaplBlocks=max(lenHaplBlocks))


allChr_lengthHaplBlocks_all %>%
  filter(lenHaplBlocks>10000) %>%
  # ggplot(aes(log10(lenHaplBlocks)))+
  ggplot(aes(lenHaplBlocks/1000000))+
  # geom_density(fill="black")+
  geom_histogram(fill="#0033CC")+
  geom_vline(data=meanChr_allChr_lengthHaplBlocks_all, 
             aes(xintercept = mean_lenHaplBlocks/1000000), 
             #linetype="dotted", 
             color = "red") +
  # geom_vline(data=meanChr_allChr_lengthHaplBlocks_all, 
  #            aes(xintercept = max_lenHaplBlocks/1000000), 
  #            #linetype="dotted", 
  #            color = "blue") +
  # scale_x_continuous(breaks=seq(0,1000,200))+
  scale_x_continuous(trans="log10") +
  labs(x = "Length Shared Blocks (Mb)",
       y = "Count") +
  facet_grid(chr~., scale="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey98"), 
        panel.grid.minor = element_line(color = "grey98"))
# axis.text = element_text(size=16), 
# axis.title = element_text(size=16))

ggsave(paste0("09_a2_RawShareHaplLenBlocks_perChr_ed_10kbWin_min10SS.png"), height=10, width=2.5)
ggsave(paste0("09_a2_RawShareHaplLenBlocks_perChr_ed_10kbWin_min10SS.pdf"), height=10, width=2.5)


allChr_lengthHaplBlocks_all %>%
  group_by(chr, haplotype1, haplotype2) %>%
  summarise(max_lenHaplBlocks=max(lenHaplBlocks)) %>%
  ggplot(aes(max_lenHaplBlocks/1000000))+
  geom_histogram(fill="#990000")+
  scale_x_continuous(trans="log10") +
  labs(x = "Max Length Shared Blocks (Mb)",
       y = "Count") +
  facet_grid(chr~., scale="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey98"), 
        panel.grid.minor = element_line(color = "grey98"))


ggsave(paste0("09_a3_MaxRawShareHaplLenBlocks_perChr_ed_10kbWin_min10SS.png"), height=10, width=2.5)
ggsave(paste0("09_a3_MaxRawShareHaplLenBlocks_perChr_ed_10kbWin_min10SS.pdf"), height=10, width=2.5)




# plots with average per pairwise comparison:
chrN="01"
for (chrN in chr_list) {
  lengthHaplBlocks_all <- read.csv(paste0("Table_SharedBlocks_chr",chrN,".csv"), F)
  names(lengthHaplBlocks_all) <- c("chr", "haplotype1", "haplotype2", "startpos", "endpos")
  
  lengthHaplBlocks_all %>%
    mutate(lenHaplBlocks=endpos-startpos) %>%
    group_by(chr, haplotype1, haplotype2) %>%
    summarise(NSharedBlocks=n(), 
              totalLenShared=sum(lenHaplBlocks), 
              meanLenShared=mean(lenHaplBlocks)) %>%
    left_join(DM_ref_chrSizes, by="chr") %>%
    mutate(frac=totalLenShared/totalLen) %>%
    write.table(file=paste0("Table_PropSharedBlocks.csv"), quote = F, row.names = F, col.names = F, sep = ",", append=T)
}




lengthHaplBlocks_all <- read.csv(paste0("Table_PropSharedBlocks.csv"), F)
names(lengthHaplBlocks_all) <- c("chr", "haplotype1", "haplotype2", "NSharedBlocks", "totalLenShared", "meanLen", "totalLen", "frac")
head(lengthHaplBlocks_all)

lengthHaplBlocks_all_meanValues <- lengthHaplBlocks_all %>%
  group_by(chr) %>%
  summarise(mean_frac=mean(frac), 
            min_frac=min(frac), 
            max_frac=max(frac)) 

lengthHaplBlocks_all_meanValues

lengthHaplBlocks_all_meanValues %>%
  arrange(mean_frac) %>%
  pull(mean_frac) %>%
  summary()

# chr   mean_frac min_frac max_frac
# chr01     0.193   0.0449    0.826
# chr02     0.181   0.0258    0.801
# chr03     0.250   0.0721    0.677
# chr04     0.189   0.0423    0.808
# chr05     0.233   0.0254    0.924
# chr06     0.214   0.0357    0.770
# chr07     0.181   0.0182    0.940
# chr08     0.194   0.0377    0.885
# chr09     0.212   0.0216    0.914
# chr10     0.402   0.0586    0.817
# chr11     0.199   0.0383    0.885
# chr12     0.236   0.0191    0.882


lengthHaplBlocks_all %>%
  ggplot(aes(frac))+
  geom_histogram(fill="#CC6633")+
  geom_vline(data=lengthHaplBlocks_all_meanValues, aes(xintercept = mean_frac), 
             #linetype="dotted", 
             color = "red") +
  labs(x = "Total Length Shared Fraction",
       y = "Pair-wise Haplotype Count") +
  facet_grid(chr~., scale="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey98"), 
        panel.grid.minor = element_line(color = "grey98"))

ggsave(paste0("09_MeanShareHaplBlocks_10kbWin_min10SS.png"), height=10, width=2.5)
ggsave(paste0("09_MeanShareHaplBlocks_10kbWin_min10SS.pdf"), height=10, width=2.5)



lengthHaplBlocks_all %>%
  ggplot(aes(frac))+
  geom_histogram(fill="black")+
  geom_vline(data=lengthHaplBlocks_all_meanValues, aes(xintercept = mean_frac), 
             #linetype="dotted", 
             color = "red") +
  scale_x_continuous(breaks=seq(0,1,0.1))+
  labs(x = "Total Length Shared Fraction",
       y = "Pair-wise Haplotype Count") +
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor = element_line(color = "grey80"),
        axis.text = element_text(size=24), 
        axis.title = element_text(size=24))

ggsave(paste0("09_MeanShareHaplBlocks_ed_10kbWin_min10SS.png"), height=5, width=7)
ggsave(paste0("09_MeanShareHaplBlocks_ed_10kbWin_min10SS.pdf"), height=5, width=7)







lengthHaplBlocks_all_meanValues <- lengthHaplBlocks_all %>%
  group_by(chr) %>%
  summarise(mean_Len=mean(meanLen), 
            min_Len=min(meanLen), 
            max_Len=max(meanLen)) 


lengthHaplBlocks_all %>%
  ggplot(aes(log10(meanLen)))+
  geom_histogram(fill="black")+
  geom_vline(data=lengthHaplBlocks_all_meanValues, aes(xintercept = log10(mean_Len)), 
             #linetype="dotted", 
             color = "red") +
  labs(x = "Log10 Mean Length Shared Block",
       y = "Pair-wise Haplotype Count") +
  facet_grid(chr~., scale="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))

lengthHaplBlocks_all %>% 
  ggplot(aes((meanLen/1000)))+
  geom_histogram(fill="black")+
  geom_vline(data=lengthHaplBlocks_all_meanValues, aes(xintercept = (mean_Len/1000)), 
             #linetype="dotted", 
             color = "red") +
  labs(x = "Mean Length Shared Block (kb)",
       y = "Pair-wise Haplotype Count") +
  scale_x_continuous(breaks=seq(0,1000,200))+
  facet_grid(chr~., scale="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))

ggsave(paste0("09_b_MeanShareHaplLenBlocks_10kbWin_min10SS.png"), height=10, width=3)
ggsave(paste0("09_b_MeanShareHaplLenBlocks_10kbWin_min10SS.pdf"), height=10, width=3)



lengthHaplBlocks_all %>%
  ggplot(aes(meanLen/1000))+
  geom_histogram(fill="black")+
  geom_vline(data=lengthHaplBlocks_all_meanValues, aes(xintercept = mean_Len/1000), 
             #linetype="dotted", 
             color = "red") +
  scale_x_continuous(breaks=seq(0,1000,200))+
  labs(x = "Mean Length Shared Block (Kb)",
       y = "Pair-wise Haplotype Count") +
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor = element_line(color = "grey80"),
        axis.text = element_text(size=24), 
        axis.title = element_text(size=24))

ggsave(paste0("09_b_MeanShareHaplLenBlocks_ed_10kbWin_min10SS.png"), height=5, width=7)
ggsave(paste0("09_b_MeanShareHaplLenBlocks_ed_10kbWin_min10SS.pdf"), height=5, width=7)


```



## Proportion shared within individuals



```R

# .libPaths("/dss/dsshome1/lxc03/di36guz2/R")
# mamba activate env_others
.libPaths()


setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/05_diversity_withinInd/ROH")
#setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/03_Nhaplotypes")

library(ggplot2)
library(tidyverse)

chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

DM_ref_chrSizes <- data.frame(chr = c("chr01", "chr02", "chr03", "chr04", "chr05", "chr06", "chr07", "chr08", "chr09", "chr10", "chr11", "chr12"), 
                              totalLen = c(88591686, 46102915, 60707570, 69236331, 55599697, 59091578, 57639317, 59226000, 67600300, 61044151, 46777387, 59670755))

chrN <- "01"


for (chrN in chr_list) {
  #print(chrN) }
  #lengthHaplBlocks_all = data.frame(chr="tem", haplotype1="tem", haplotype2="tem", startpos=0, endpos=0)
  table_hapG <- read_table(paste0("../../03_Nhaplotypes/haplotypeGroups_min10SS_chr",chrN,".txt"), F)
  names(table_hapG) <- c("pos", "sample", "hapGroup", "Mean_NAN_haplo")
  
  table_hapG_spread <-  table_hapG %>% 
    select(-Mean_NAN_haplo) %>%
    spread(key=sample, value=hapGroup) %>%
    arrange(pos)
  
  list_samples <- unique(table_hapG$sample)
  
  for (hapi in seq(1,length(list_samples), 4)){
    for (i in seq(hapi,(hapi+3))){
      for (j in seq(hapi,(hapi+3))){
        if (j > i){
          table_hapG_2Hapl <- table_hapG_spread %>%
            select(pos, list_samples[i], list_samples[j]) 
          
          sameHapl <- as.vector(((table_hapG_2Hapl[,2] == table_hapG_2Hapl[,3])*1))
          
          lengthHaplBlocks <- table_hapG_2Hapl %>%
            mutate(sameHapl) %>%
            arrange(pos) %>%
            #mutate(cont_haplo=(sameHapl==lag(sameHapl) & sameHapl==1)*1) 
            mutate(cont_haplo=(sameHapl + lag(sameHapl) > 0)*1) %>%
            mutate(sequence = data.table::rleid(cont_haplo == 1)) %>%
            filter(sameHapl==1) %>%
            filter(cont_haplo == 1) %>%
            group_by(sequence) %>%
            summarise(startpos = min(pos), 
                      endpos=max(pos)) %>%
            mutate(chr=paste0("chr", chrN), 
                   sample=str_split(list_samples[i], "_")[[1]][1],
                   haplotype1=list_samples[i], 
                   haplotype2=list_samples[j]) %>%
            select(chr, sample, haplotype1, haplotype2, startpos, endpos)
          write.table(lengthHaplBlocks, file=paste0("Table_SharedBlocksWithinInd_chr", chrN, ".csv"), 
                      quote = F, row.names = F, col.names = F, sep = ",", append=T)
        }
      }
    }
  }
}





chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
chrN <- "01"
for (hapi in seq(1,500)){
  print(hapi)
  haplotypesUsed <- sample(1:40, 4, replace=F)
  for (chrN in chr_list) {
    #print(chrN) }
    #lengthHaplBlocks_all = data.frame(chr="tem", haplotype1="tem", haplotype2="tem", startpos=0, endpos=0)
    table_hapG <- read_table(paste0("../../03_Nhaplotypes/haplotypeGroups_min10SS_chr",chrN,".txt"), F)
    names(table_hapG) <- c("pos", "sample", "hapGroup", "Mean_NAN_haplo")
    
    list_samples <- unique(table_hapG$sample)
    
    table_hapG_spread <-  table_hapG %>% 
      select(-Mean_NAN_haplo) %>%
      spread(key=sample, value=hapGroup) %>%
      arrange(pos)
    
    for (i in seq(1,4)){
      for (j in seq(1,4)){
        if (j > i){
          table_hapG_2Hapl <- table_hapG_spread %>%
            select(pos, list_samples[haplotypesUsed[i]], list_samples[haplotypesUsed[j]]) 
          
          sameHapl <- as.vector(((table_hapG_2Hapl[,2] == table_hapG_2Hapl[,3])*1))
          
          lengthHaplBlocks <- table_hapG_2Hapl %>%
            mutate(sameHapl) %>%
            arrange(pos) %>%
            #mutate(cont_haplo=(sameHapl==lag(sameHapl) & sameHapl==1)*1) 
            mutate(cont_haplo=(sameHapl + lag(sameHapl) > 0)*1) %>%
            mutate(sequence = data.table::rleid(cont_haplo == 1)) %>%
            filter(sameHapl==1) %>%
            filter(cont_haplo == 1) %>%
            group_by(sequence) %>%
            summarise(startpos = min(pos), 
                      endpos=max(pos)) %>%
            mutate(chr=paste0("chr", chrN), 
                   sample=hapi,
                   haplotype1=list_samples[haplotypesUsed[i]], 
                   haplotype2=list_samples[haplotypesUsed[j]]) %>%
            select(chr, sample, haplotype1, haplotype2, startpos, endpos)
          write.table(lengthHaplBlocks, file=paste0("RandomHap_Table_SharedBlocksWithinInd_chr", chrN, ".csv"), 
                      quote = F, row.names = F, col.names = F, sep = ",", append=T)
        }
      }
    }
  }
}



chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
chrN="01"
for (chrN in chr_list) {
  lengthHaplBlocks_all <- read.csv(paste0("Table_SharedBlocksWithinInd_chr",chrN,".csv"), F)
  names(lengthHaplBlocks_all) <- c("chr", "sample", "haplotype1", "haplotype2", "startpos", "endpos")
  
  lengthHaplBlocks_all %>%
    mutate(lenHaplBlocks=endpos-startpos) %>%
    group_by(chr, sample, haplotype1, haplotype2) %>%
    summarise(NSharedBlocks=n(), 
              totalLenShared=sum(lenHaplBlocks), 
              mean_blockLength=mean(lenHaplBlocks)) %>%
    left_join(DM_ref_chrSizes, by="chr") %>%
    mutate(frac=totalLenShared/totalLen) %>%
    write.table(file=paste0("AllHap_Table_PropSharedBlocks.csv"), 
                quote = F, row.names = F, col.names = F, sep = ",", append=T)
  
  lengthHaplBlocks_all %>%
    mutate(lenHaplBlocks=endpos-startpos) %>%
    group_by(chr, sample, haplotype1, haplotype2) %>%
    summarise(NSharedBlocks=n(), 
              totalLenShared=sum(lenHaplBlocks), 
              mean_blockLength=mean(lenHaplBlocks)) %>%
    group_by(chr, sample) %>%
    summarise(mean_NSharedBlocks=mean(NSharedBlocks), 
              mean_totalLenShared=mean(totalLenShared), 
              sample_mean_blockLength=mean(mean_blockLength)) %>%
    left_join(DM_ref_chrSizes, by="chr") %>%
    mutate(mean_frac=mean_totalLenShared/totalLen) %>%
    write.table(file=paste0("AllperSample_Table_PropSharedBlocks.csv"), 
                quote = F, row.names = F, col.names = F, sep = ",", append=T)
}




chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
chrN="01"
for (chrN in chr_list) {
  lengthHaplBlocks_all <- read.csv(paste0("RandomHap_Table_SharedBlocksWithinInd_chr",chrN,".csv"), F)
  names(lengthHaplBlocks_all) <- c("chr", "sample", "haplotype1", "haplotype2", "startpos", "endpos")
  
  lengthHaplBlocks_all %>%
    mutate(lenHaplBlocks=endpos-startpos) %>%
    group_by(chr, sample, haplotype1, haplotype2) %>%
    summarise(NSharedBlocks=n(), 
              totalLenShared=sum(lenHaplBlocks), 
              mean_blockLength=mean(lenHaplBlocks)) %>%
    left_join(DM_ref_chrSizes, by="chr") %>%
    mutate(frac=totalLenShared/totalLen) %>%
    write.table(file=paste0("AllRandomHap_AllHap_Table_PropSharedBlocks.csv"), 
                quote = F, row.names = F, col.names = F, sep = ",", append=T)
  
  lengthHaplBlocks_all %>%
    mutate(lenHaplBlocks=endpos-startpos) %>%
    group_by(chr, sample, haplotype1, haplotype2) %>%
    summarise(NSharedBlocks=n(), 
              totalLenShared=sum(lenHaplBlocks), 
              mean_blockLength=mean(lenHaplBlocks)) %>%
    group_by(chr, sample) %>%
    summarise(mean_NSharedBlocks=mean(NSharedBlocks), 
              mean_totalLenShared=mean(totalLenShared), 
              sample_mean_blockLength=mean(mean_blockLength)) %>%
    left_join(DM_ref_chrSizes, by="chr") %>%
    mutate(mean_frac=mean_totalLenShared/totalLen) %>%
    write.table(file=paste0("AllRandomHap_AllperSample_Table_PropSharedBlocks.csv"), 
                quote = F, row.names = F, col.names = F, sep = ",", append=T)
}






lengthHaplBlocks_all <- read.csv(paste0("AllperSample_Table_PropSharedBlocks.csv"), F)
#names(lengthHaplBlocks_all) <- c("chr", "sample", "haplotype1", "haplotype2", "NSharedBlocks", "totalLenShared", "totalLen", "frac")
names(lengthHaplBlocks_all) <- c("chr", "sample", "meanNSharedBlocks", "meantotalLenShared", "meanblockLength", "chrSize", "meanfrac")



random_lengthHaplBlocks_all <- read.csv(paste0("AllRandomHap_AllperSample_Table_PropSharedBlocks.csv"), F)
#names(lengthHaplBlocks_all) <- c("chr", "sample", "haplotype1", "haplotype2", "NSharedBlocks", "totalLenShared", "totalLen", "frac")
names(random_lengthHaplBlocks_all) <- c("chr", "sample", "meanNSharedBlocks", "meantotalLenShared", "meanblockLength", "chrSize", "meanfrac")




random_lengthHaplBlocks_all %>%
  ggplot(aes(meanfrac))+
  geom_histogram(fill="black")+
  geom_vline(data=lengthHaplBlocks_all, aes(xintercept = meanfrac), 
             #linetype="dotted", 
             color = "red") +
  labs(x = "Total Length Shared Fraction",
       y = "") +
  facet_grid(chr~., scale="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))

ggsave(paste0("01_MeanShareHaplBlocks_perInd_10kbWin_min10SS.png"), height=10, width=5)
ggsave(paste0("01_MeanShareHaplBlocks_perInd_10kbWin_min10SS.pdf"), height=10, width=5)



random_lengthHaplBlocks_all %>%
  ggplot(aes(meantotalLenShared))+
  geom_histogram(fill="black")+
  geom_vline(data=lengthHaplBlocks_all, aes(xintercept = meantotalLenShared), 
             #linetype="dotted", 
             color = "red") +
  labs(x = "Mean Total Length Shared",
       y = "") +
  facet_grid(chr~., scale="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))





random_lengthHaplBlocks_all %>%
  filter(meanblockLength<300000) %>%
  ggplot(aes(meanblockLength))+
  geom_histogram(fill="black")+
  geom_vline(data=lengthHaplBlocks_all, aes(xintercept = meanblockLength), 
             #linetype="dotted", 
             color = "red") +
  labs(x = "Mean Shared Block Length",
       y = "") +
  facet_grid(chr~., scale="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))

ggsave(paste0("02_MeanShareBlockLenght_perInd_10kbWin_min10SS.png"), height=10, width=5)
ggsave(paste0("02_MeanShareBlockLenght_perInd_10kbWin_min10SS.pdf"), height=10, width=5)




random_lengthHaplBlocks_all %>%
  ggplot(aes(meanNSharedBlocks))+
  geom_histogram(fill="black")+
  geom_vline(data=lengthHaplBlocks_all, aes(xintercept = meanNSharedBlocks), 
             #linetype="dotted", 
             color = "red") +
  labs(x = "Mean Number Shared Blocks",
       y = "") +
  facet_grid(chr~., scale="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))

ggsave(paste0("03_MeanShareBlockLenght_perInd_10kbWin_min10SS.png"), height=10, width=5)
ggsave(paste0("03_MeanShareBlockLenght_perInd_10kbWin_min10SS.pdf"), height=10, width=5)






random_lengthHaplBlocks_all %>%
  group_by(sample) %>%
  summarise(meanfrac=mean(meanfrac)) %>%
  ggplot(aes(meanfrac))+
  geom_histogram(fill="black")+
  geom_vline(data=lengthHaplBlocks_all %>%
               group_by(sample) %>%
               summarise(meanfrac=mean(meanfrac)), aes(xintercept = meanfrac), 
             #linetype="dotted", 
             color = "red") +
  labs(x = "Total Length Shared Fraction",
       y = "") +
  # facet_grid(chr~., scale="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))

ggsave(paste0("04_totalMeanShareHaplBlocks_perInd_10kbWin_min10SS.png"), height=5, width=5)
ggsave(paste0("04_totalMeanShareHaplBlocks_perInd_10kbWin_min10SS.pdf"), height=5, width=5)

```



## Proportion Shared within Individuals - removing GAPs



```R

# .libPaths("/dss/dsshome1/lxc03/di36guz2/R")
# mamba activate env_others
.libPaths()


setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/05_diversity_withinInd/ROH")
#setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/03_Nhaplotypes")

library(ggplot2)
library(tidyverse)

chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

DM_ref_chrSizes <- data.frame(chr = c("chr01", "chr02", "chr03", "chr04", "chr05", "chr06", "chr07", "chr08", "chr09", "chr10", "chr11", "chr12"), 
                              totalLen = c(88591686, 46102915, 60707570, 69236331, 55599697, 59091578, 57639317, 59226000, 67600300, 61044151, 46777387, 59670755))

chrN <- "01"

for (chrN in chr_list) {
  #print(chrN) }
  #lengthHaplBlocks_all = data.frame(chr="tem", haplotype1="tem", haplotype2="tem", startpos=0, endpos=0)
  table_hapG <- read_table(paste0("../../03_Nhaplotypes/haplotypeGroups_min10SS_chr",chrN,".txt"), F)
  names(table_hapG) <- c("pos", "sample", "hapGroup", "Mean_NAN_haplo")
  
  table_hapG_spread <-  table_hapG %>% 
    filter(Mean_NAN_haplo<0.8) %>%
    select(-Mean_NAN_haplo) %>%
    spread(key=sample, value=hapGroup, fill = NA) %>%
    arrange(pos)
  
  list_samples <- unique(table_hapG$sample)
  
  for (hapi in seq(1,length(list_samples), 4)){
    for (i in seq(hapi,(hapi+3))){
      for (j in seq(hapi,(hapi+3))){
        if (j > i){
          table_hapG_2Hapl <- table_hapG_spread %>%
            select(pos, list_samples[i], list_samples[j]) 
          
          sameHapl <- as.vector((((table_hapG_2Hapl[,2] == table_hapG_2Hapl[,3]) & !(is.na(table_hapG_2Hapl[,2])) & !(is.na(table_hapG_2Hapl[,3])))*1))
          
          lengthHaplBlocks <- table_hapG_2Hapl %>%
            mutate(sameHapl) %>%
            #mutate(cont_haplo=(sameHapl==lag(sameHapl) & sameHapl==1)*1) 
            mutate(cont_haplo=(sameHapl + lag(sameHapl) > 0)*1) %>%
            mutate(sequence = data.table::rleid(cont_haplo == 1)) %>% 
            filter(sameHapl==1) %>%
            filter(cont_haplo == 1) %>% 
            group_by(sequence) %>%
            summarise(startpos = min(pos), 
                      endpos=max(pos)) %>%
            mutate(chr=paste0("chr", chrN), 
                   sample=str_split(list_samples[i], "_")[[1]][1],
                   haplotype1=list_samples[i], 
                   haplotype2=list_samples[j]) %>%
            select(chr, sample, haplotype1, haplotype2, startpos, endpos)
          write.table(lengthHaplBlocks, file=paste0("RemGapsTable_SharedBlocksWithinInd_chr", chrN, ".csv"), 
                      quote = F, row.names = F, col.names = F, sep = ",", append=T)
        }
      }
    }
  }
}





chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
chrN <- "01"
for (hapi in seq(1,500)){
  print(hapi)
  haplotypesUsed <- sample(1:40, 4, replace=F)
  for (chrN in chr_list) {
    #print(chrN) }
    #lengthHaplBlocks_all = data.frame(chr="tem", haplotype1="tem", haplotype2="tem", startpos=0, endpos=0)
    table_hapG <- read_table(paste0("../../03_Nhaplotypes/haplotypeGroups_min10SS_chr",chrN,".txt"), F)
    names(table_hapG) <- c("pos", "sample", "hapGroup", "Mean_NAN_haplo")
    
    list_samples <- unique(table_hapG$sample)
    
    table_hapG_spread <-  table_hapG %>% 
      filter(Mean_NAN_haplo<0.8) %>%
      select(-Mean_NAN_haplo) %>%
      spread(key=sample, value=hapGroup, fill = NA) %>%
      arrange(pos)
    
    for (i in seq(1,4)){
      for (j in seq(1,4)){
        if (j > i){
          table_hapG_2Hapl <- table_hapG_spread %>%
            select(pos, list_samples[haplotypesUsed[i]], list_samples[haplotypesUsed[j]]) 
          
          sameHapl <- as.vector((((table_hapG_2Hapl[,2] == table_hapG_2Hapl[,3]) & !(is.na(table_hapG_2Hapl[,2])) & !(is.na(table_hapG_2Hapl[,3])))*1))
          
          lengthHaplBlocks <- table_hapG_2Hapl %>%
            mutate(sameHapl) %>%
            arrange(pos) %>%
            #mutate(cont_haplo=(sameHapl==lag(sameHapl) & sameHapl==1)*1) 
            mutate(cont_haplo=(sameHapl + lag(sameHapl) > 0)*1) %>%
            mutate(sequence = data.table::rleid(cont_haplo == 1)) %>%
            filter(sameHapl==1) %>%
            filter(cont_haplo == 1) %>%
            group_by(sequence) %>%
            summarise(startpos = min(pos), 
                      endpos=max(pos)) %>%
            mutate(chr=paste0("chr", chrN), 
                   sample=hapi,
                   haplotype1=list_samples[haplotypesUsed[i]], 
                   haplotype2=list_samples[haplotypesUsed[j]]) %>%
            select(chr, sample, haplotype1, haplotype2, startpos, endpos)
          write.table(lengthHaplBlocks, file=paste0("RemGapsRandomHap_Table_SharedBlocksWithinInd_chr", chrN, ".csv"), 
                      quote = F, row.names = F, col.names = F, sep = ",", append=T)
        }
      }
    }
  }
}







chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
chrN="01"
for (chrN in chr_list) {
  lengthHaplBlocks_all <- read.csv(paste0("RemGapsTable_SharedBlocksWithinInd_chr",chrN,".csv"), F)
  names(lengthHaplBlocks_all) <- c("chr", "sample", "haplotype1", "haplotype2", "startpos", "endpos")
  
  lengthHaplBlocks_all %>%
    mutate(lenHaplBlocks=endpos-startpos) %>%
    group_by(chr, sample, haplotype1, haplotype2) %>%
    summarise(NSharedBlocks=n(), 
              totalLenShared=sum(lenHaplBlocks), 
              mean_blockLength=mean(lenHaplBlocks)) %>%
    left_join(DM_ref_chrSizes, by="chr") %>%
    mutate(frac=totalLenShared/totalLen) %>%
    write.table(file=paste0("AllHap_RemGaps_Table_PropSharedBlocks.csv"), 
                quote = F, row.names = F, col.names = F, sep = ",", append=T)
  
  lengthHaplBlocks_all %>%
    mutate(lenHaplBlocks=endpos-startpos) %>%
    group_by(chr, sample, haplotype1, haplotype2) %>%
    summarise(NSharedBlocks=n(), 
              totalLenShared=sum(lenHaplBlocks), 
              mean_blockLength=mean(lenHaplBlocks)) %>%
    group_by(chr, sample) %>%
    summarise(mean_NSharedBlocks=mean(NSharedBlocks), 
              mean_totalLenShared=mean(totalLenShared), 
              sample_mean_blockLength=mean(mean_blockLength)) %>%
    left_join(DM_ref_chrSizes, by="chr") %>%
    mutate(mean_frac=mean_totalLenShared/totalLen) %>%
    write.table(file=paste0("AllperSample_RemGaps_Table_PropSharedBlocks.csv"), 
                quote = F, row.names = F, col.names = F, sep = ",", append=T)
}




chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
chrN="01"
for (chrN in chr_list) {
  lengthHaplBlocks_all <- read.csv(paste0("RemGapsRandomHap_Table_SharedBlocksWithinInd_chr",chrN,".csv"), F)
  names(lengthHaplBlocks_all) <- c("chr", "sample", "haplotype1", "haplotype2", "startpos", "endpos")
  
  lengthHaplBlocks_all %>%
    mutate(lenHaplBlocks=endpos-startpos) %>%
    group_by(chr, sample, haplotype1, haplotype2) %>%
    summarise(NSharedBlocks=n(), 
              totalLenShared=sum(lenHaplBlocks), 
              mean_blockLength=mean(lenHaplBlocks)) %>%
    left_join(DM_ref_chrSizes, by="chr") %>%
    mutate(frac=totalLenShared/totalLen) %>%
    write.table(file=paste0("AllRandomHap_RemGaps_AllHap_Table_PropSharedBlocks.csv"), 
                quote = F, row.names = F, col.names = F, sep = ",", append=T)
  
  lengthHaplBlocks_all %>%
    mutate(lenHaplBlocks=endpos-startpos) %>%
    group_by(chr, sample, haplotype1, haplotype2) %>%
    summarise(NSharedBlocks=n(), 
              totalLenShared=sum(lenHaplBlocks), 
              mean_blockLength=mean(lenHaplBlocks)) %>%
    group_by(chr, sample) %>%
    summarise(mean_NSharedBlocks=mean(NSharedBlocks), 
              mean_totalLenShared=mean(totalLenShared), 
              sample_mean_blockLength=mean(mean_blockLength)) %>%
    left_join(DM_ref_chrSizes, by="chr") %>%
    mutate(mean_frac=mean_totalLenShared/totalLen) %>%
    write.table(file=paste0("AllRandomHap_RemGaps_AllperSample_Table_PropSharedBlocks.csv"), 
                quote = F, row.names = F, col.names = F, sep = ",", append=T)
}






lengthHaplBlocks_all <- read.csv(paste0("AllperSample_RemGaps_Table_PropSharedBlocks.csv"), F)
#names(lengthHaplBlocks_all) <- c("chr", "sample", "haplotype1", "haplotype2", "NSharedBlocks", "totalLenShared", "totalLen", "frac")
names(lengthHaplBlocks_all) <- c("chr", "sample", "meanNSharedBlocks", "meantotalLenShared", "meanblockLength", "chrSize", "meanfrac")



random_lengthHaplBlocks_all <- read.csv(paste0("AllRandomHap_RemGaps_AllperSample_Table_PropSharedBlocks.csv"), F)
#names(lengthHaplBlocks_all) <- c("chr", "sample", "haplotype1", "haplotype2", "NSharedBlocks", "totalLenShared", "totalLen", "frac")
names(random_lengthHaplBlocks_all) <- c("chr", "sample", "meanNSharedBlocks", "meantotalLenShared", "meanblockLength", "chrSize", "meanfrac")




random_lengthHaplBlocks_all %>%
  ggplot(aes(meanfrac))+
  geom_histogram(fill="black")+
  geom_vline(data=lengthHaplBlocks_all, aes(xintercept = meanfrac), 
             #linetype="dotted", 
             color = "red") +
  labs(x = "Total Length Shared Fraction",
       y = "") +
  facet_grid(chr~., scale="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))

ggsave(paste0("05_RemGaps_MeanShareHaplBlocks_perInd_10kbWin_min10SS.png"), height=10, width=5)
ggsave(paste0("05_RemGaps_MeanShareHaplBlocks_perInd_10kbWin_min10SS.pdf"), height=10, width=5)



random_lengthHaplBlocks_all %>%
  ggplot(aes(meantotalLenShared))+
  geom_histogram(fill="black")+
  geom_vline(data=lengthHaplBlocks_all, aes(xintercept = meantotalLenShared), 
             #linetype="dotted", 
             color = "red") +
  labs(x = "Mean Total Length Shared",
       y = "") +
  facet_grid(chr~., scale="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))





random_lengthHaplBlocks_all %>%
  filter(meanblockLength<300000) %>%
  ggplot(aes(meanblockLength))+
  geom_histogram(fill="black")+
  geom_vline(data=lengthHaplBlocks_all, aes(xintercept = meanblockLength), 
             #linetype="dotted", 
             color = "red") +
  labs(x = "Mean Shared Block Length",
       y = "") +
  facet_grid(chr~., scale="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))

ggsave(paste0("06_RemGaps_MeanShareBlockLenght_perInd_10kbWin_min10SS.png"), height=10, width=5)
ggsave(paste0("06_RemGaps_MeanShareBlockLenght_perInd_10kbWin_min10SS.pdf"), height=10, width=5)




random_lengthHaplBlocks_all %>%
  ggplot(aes(meanNSharedBlocks))+
  geom_histogram(fill="black")+
  geom_vline(data=lengthHaplBlocks_all, aes(xintercept = meanNSharedBlocks), 
             #linetype="dotted", 
             color = "red") +
  labs(x = "Mean Number Shared Blocks",
       y = "") +
  facet_grid(chr~., scale="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))

ggsave(paste0("07_RemGaps_MeanShareBlockLenght_perInd_10kbWin_min10SS.png"), height=10, width=5)
ggsave(paste0("07_RemGaps_MeanShareBlockLenght_perInd_10kbWin_min10SS.pdf"), height=10, width=5)






random_lengthHaplBlocks_all %>%
  group_by(sample) %>%
  summarise(meanfrac=mean(meanfrac)) %>%
  ggplot(aes(meanfrac))+
  geom_histogram(fill="black")+
  geom_vline(data=lengthHaplBlocks_all %>%
               group_by(sample) %>%
               summarise(meanfrac=mean(meanfrac)), aes(xintercept = meanfrac), 
             #linetype="dotted", 
             color = "red", alpha=0.5) +
  labs(x = "Total Length Shared Fraction",
       y = "") +
  # facet_grid(chr~., scale="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))

ggsave(paste0("08_RemGaps_totalMeanShareHaplBlocks_perInd_10kbWin_min10SS.png"), height=5, width=5)
ggsave(paste0("08_RemGaps_totalMeanShareHaplBlocks_perInd_10kbWin_min10SS.pdf"), height=5, width=5)

```



## Proportion Shared within Individuals - removing GAPs - removing pericentromeric regions



```R

# .libPaths("/dss/dsshome1/lxc03/di36guz2/R")
# mamba activate env_others
.libPaths()


setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/05_diversity_withinInd/ROH")
#setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/03_Nhaplotypes")

library(ggplot2)
library(tidyverse)

chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

DM_ref_chrSizes <- data.frame(chr = c("chr01", "chr02", "chr03", "chr04", "chr05", "chr06", "chr07", "chr08", "chr09", "chr10", "chr11", "chr12"), 
                              totalLen = c(88591686, 46102915, 60707570, 69236331, 55599697, 59091578, 57639317, 59226000, 67600300, 61044151, 46777387, 59670755))

centromere_table <- read.csv("centromeres.csv", header = T)

chrN <- "01"

for (chrN in chr_list) {
  #print(chrN) }
  #lengthHaplBlocks_all = data.frame(chr="tem", haplotype1="tem", haplotype2="tem", startpos=0, endpos=0)
  table_hapG <- read_table(paste0("../../03_Nhaplotypes/haplotypeGroups_min10SS_chr",chrN,".txt"), F)
  names(table_hapG) <- c("pos", "sample", "hapGroup", "Mean_NAN_haplo")
  
  # head(table_hapG)
  chrValues <- centromere_table %>%
    filter(chr==paste0("chr", chrN))
  
  table_hapG_spread <-  table_hapG %>%
    filter(Mean_NAN_haplo<0.8) %>%
    filter((pos<chrValues$start) | (pos>chrValues$end)) %>%
    select(-Mean_NAN_haplo) %>%
    spread(key=sample, value=hapGroup, fill = NA) %>%
    arrange(pos)
  
  list_samples <- unique(table_hapG$sample)
  
  for (hapi in seq(1,length(list_samples), 4)){
    for (i in seq(hapi,(hapi+3))){
      for (j in seq(hapi,(hapi+3))){
        if (j > i){
          table_hapG_2Hapl <- table_hapG_spread %>%
            select(pos, list_samples[i], list_samples[j]) 
          
          sameHapl <- as.vector((((table_hapG_2Hapl[,2] == table_hapG_2Hapl[,3]) & !(is.na(table_hapG_2Hapl[,2])) & !(is.na(table_hapG_2Hapl[,3])))*1))
          
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
                   sample=str_split(list_samples[i], "_")[[1]][1],
                   haplotype1=list_samples[i], 
                   haplotype2=list_samples[j]) %>%
            select(chr, sample, haplotype1, haplotype2, startpos, endpos)
          write.table(lengthHaplBlocks, file=paste0("RemCentromeres_RemGapsTable_SharedBlocksWithinInd_chr", chrN, ".csv"), 
                      quote = F, row.names = F, col.names = F, sep = ",", append=T)
        }
      }
    }
  }
}





chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
chrN <- "01"
for (hapi in seq(1,500)){
  print(hapi)
  haplotypesUsed <- sample(1:40, 4, replace=F)
  for (chrN in chr_list) {
    #print(chrN) }
    #lengthHaplBlocks_all = data.frame(chr="tem", haplotype1="tem", haplotype2="tem", startpos=0, endpos=0)
    table_hapG <- read_table(paste0("../../03_Nhaplotypes/haplotypeGroups_min10SS_chr",chrN,".txt"), F)
    names(table_hapG) <- c("pos", "sample", "hapGroup", "Mean_NAN_haplo")
    
    list_samples <- unique(table_hapG$sample)
    
    chrValues <- centromere_table %>%
      filter(chr==paste0("chr", chrN))
    
    table_hapG_spread <-  table_hapG %>% 
      filter(Mean_NAN_haplo<0.8) %>%
      filter((pos<chrValues$start) | (pos>chrValues$end)) %>%
      select(-Mean_NAN_haplo) %>%
      spread(key=sample, value=hapGroup, fill = NA) %>%
      arrange(pos)
    
    for (i in seq(1,4)){
      for (j in seq(1,4)){
        if (j > i){
          table_hapG_2Hapl <- table_hapG_spread %>%
            select(pos, list_samples[haplotypesUsed[i]], list_samples[haplotypesUsed[j]]) 
          
          sameHapl <- as.vector((((table_hapG_2Hapl[,2] == table_hapG_2Hapl[,3]) & !(is.na(table_hapG_2Hapl[,2])) & !(is.na(table_hapG_2Hapl[,3])))*1))
          
          lengthHaplBlocks <- table_hapG_2Hapl %>%
            mutate(sameHapl) %>%
            arrange(pos) %>%
            #mutate(cont_haplo=(sameHapl==lag(sameHapl) & sameHapl==1)*1) 
            mutate(cont_haplo=(sameHapl + lag(sameHapl) > 0)*1) %>%
            mutate(sequence = data.table::rleid(cont_haplo == 1)) %>%
            filter(cont_haplo == 1) %>%
            group_by(sequence) %>%
            summarise(startpos = min(pos), 
                      endpos=max(pos)) %>%
            mutate(chr=paste0("chr", chrN), 
                   sample=hapi,
                   haplotype1=list_samples[haplotypesUsed[i]], 
                   haplotype2=list_samples[haplotypesUsed[j]]) %>%
            select(chr, sample, haplotype1, haplotype2, startpos, endpos)
          write.table(lengthHaplBlocks, file=paste0("RemCentromeres_RemGapsRandomHap_Table_SharedBlocksWithinInd_chr", chrN, ".csv"), 
                      quote = F, row.names = F, col.names = F, sep = ",", append=T)
        }
      }
    }
  }
}







chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
chrN="01"
for (chrN in chr_list) {
  lengthHaplBlocks_all <- read.csv(paste0("RemCentromeres_RemGapsTable_SharedBlocksWithinInd_chr",chrN,".csv"), F)
  names(lengthHaplBlocks_all) <- c("chr", "sample", "haplotype1", "haplotype2", "startpos", "endpos")
  
  lengthHaplBlocks_all %>%
    mutate(lenHaplBlocks=endpos-startpos) %>%
    group_by(chr, sample, haplotype1, haplotype2) %>%
    summarise(NSharedBlocks=n(), 
              totalLenShared=sum(lenHaplBlocks), 
              mean_blockLength=mean(lenHaplBlocks)) %>%
    left_join(DM_ref_chrSizes, by="chr") %>%
    mutate(frac=totalLenShared/totalLen) %>%
    write.table(file=paste0("AllHap_RemCentromeres_RemGaps_Table_PropSharedBlocks.csv"), 
                quote = F, row.names = F, col.names = F, sep = ",", append=T)
  
  lengthHaplBlocks_all %>%
    mutate(lenHaplBlocks=endpos-startpos) %>%
    group_by(chr, sample, haplotype1, haplotype2) %>%
    summarise(NSharedBlocks=n(), 
              totalLenShared=sum(lenHaplBlocks), 
              mean_blockLength=mean(lenHaplBlocks)) %>%
    group_by(chr, sample) %>%
    summarise(mean_NSharedBlocks=mean(NSharedBlocks), 
              mean_totalLenShared=mean(totalLenShared), 
              sample_mean_blockLength=mean(mean_blockLength)) %>%
    left_join(DM_ref_chrSizes, by="chr") %>%
    mutate(mean_frac=mean_totalLenShared/totalLen) %>%
    write.table(file=paste0("AllperSample_RemCentromeres_RemGaps_Table_PropSharedBlocks.csv"), 
                quote = F, row.names = F, col.names = F, sep = ",", append=T)
}




chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
chrN="01"
for (chrN in chr_list) {
  lengthHaplBlocks_all <- read.csv(paste0("RemCentromeres_RemGapsRandomHap_Table_SharedBlocksWithinInd_chr",chrN,".csv"), F)
  names(lengthHaplBlocks_all) <- c("chr", "sample", "haplotype1", "haplotype2", "startpos", "endpos")
  
  lengthHaplBlocks_all %>%
    mutate(lenHaplBlocks=endpos-startpos) %>%
    group_by(chr, sample, haplotype1, haplotype2) %>%
    summarise(NSharedBlocks=n(), 
              totalLenShared=sum(lenHaplBlocks), 
              mean_blockLength=mean(lenHaplBlocks)) %>%
    left_join(DM_ref_chrSizes, by="chr") %>%
    mutate(frac=totalLenShared/totalLen) %>%
    write.table(file=paste0("AllRandomHap_RemCentromeres_RemGaps_AllHap_Table_PropSharedBlocks.csv"), 
                quote = F, row.names = F, col.names = F, sep = ",", append=T)
  
  lengthHaplBlocks_all %>%
    mutate(lenHaplBlocks=endpos-startpos) %>%
    group_by(chr, sample, haplotype1, haplotype2) %>%
    summarise(NSharedBlocks=n(), 
              totalLenShared=sum(lenHaplBlocks), 
              mean_blockLength=mean(lenHaplBlocks)) %>%
    group_by(chr, sample) %>%
    summarise(mean_NSharedBlocks=mean(NSharedBlocks), 
              mean_totalLenShared=mean(totalLenShared), 
              sample_mean_blockLength=mean(mean_blockLength)) %>%
    left_join(DM_ref_chrSizes, by="chr") %>%
    mutate(mean_frac=mean_totalLenShared/totalLen) %>%
    write.table(file=paste0("AllRandomHap_RemCentromeres_RemGaps_AllperSample_Table_PropSharedBlocks.csv"), 
                quote = F, row.names = F, col.names = F, sep = ",", append=T)
}






lengthHaplBlocks_all <- read.csv(paste0("AllperSample_RemCentromeres_RemGaps_Table_PropSharedBlocks.csv"), F)
#names(lengthHaplBlocks_all) <- c("chr", "sample", "haplotype1", "haplotype2", "NSharedBlocks", "totalLenShared", "totalLen", "frac")
names(lengthHaplBlocks_all) <- c("chr", "sample", "meanNSharedBlocks", "meantotalLenShared", "meanblockLength", "chrSize", "meanfrac")



random_lengthHaplBlocks_all <- read.csv(paste0("AllRandomHap_RemCentromeres_RemGaps_AllperSample_Table_PropSharedBlocks.csv"), F)
#names(lengthHaplBlocks_all) <- c("chr", "sample", "haplotype1", "haplotype2", "NSharedBlocks", "totalLenShared", "totalLen", "frac")
names(random_lengthHaplBlocks_all) <- c("chr", "sample", "meanNSharedBlocks", "meantotalLenShared", "meanblockLength", "chrSize", "meanfrac")




random_lengthHaplBlocks_all %>%
  ggplot(aes(meanfrac))+
  geom_histogram(fill="black")+
  geom_vline(data=lengthHaplBlocks_all, aes(xintercept = meanfrac), 
             #linetype="dotted", 
             color = "red") +
  labs(x = "Total Length Shared Fraction",
       y = "") +
  facet_grid(chr~., scale="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))

ggsave(paste0("12_RemCentromeres_RemGaps_MeanShareHaplBlocks_perInd_10kbWin_min10SS.png"), height=10, width=5)
ggsave(paste0("12_RemCentromeres_RemGaps_MeanShareHaplBlocks_perInd_10kbWin_min10SS.pdf"), height=10, width=5)



random_lengthHaplBlocks_all %>%
  ggplot(aes(meantotalLenShared))+
  geom_histogram(fill="black")+
  geom_vline(data=lengthHaplBlocks_all, aes(xintercept = meantotalLenShared), 
             #linetype="dotted", 
             color = "red") +
  labs(x = "Mean Total Length Shared",
       y = "") +
  facet_grid(chr~., scale="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))





random_lengthHaplBlocks_all %>%
  filter(meanblockLength<300000) %>%
  ggplot(aes(meanblockLength))+
  geom_histogram(fill="black")+
  geom_vline(data=lengthHaplBlocks_all, aes(xintercept = meanblockLength), 
             #linetype="dotted", 
             color = "red") +
  labs(x = "Mean Shared Block Length",
       y = "") +
  facet_grid(chr~., scale="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))

ggsave(paste0("13_RemCentromeres_RemGaps_MeanShareBlockLenght_perInd_10kbWin_min10SS.png"), height=10, width=5)
ggsave(paste0("13_RemCentromeres_RemGaps_MeanShareBlockLenght_perInd_10kbWin_min10SS.pdf"), height=10, width=5)




random_lengthHaplBlocks_all %>%
  ggplot(aes(meanNSharedBlocks))+
  geom_histogram(fill="black")+
  geom_vline(data=lengthHaplBlocks_all, aes(xintercept = meanNSharedBlocks), 
             #linetype="dotted", 
             color = "red") +
  labs(x = "Mean Number Shared Blocks",
       y = "") +
  facet_grid(chr~., scale="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))

ggsave(paste0("14_RemCentromeres_RemGaps_MeanShareBlockLenght_perInd_10kbWin_min10SS.png"), height=10, width=5)
ggsave(paste0("14_RemCentromeres_RemGaps_MeanShareBlockLenght_perInd_10kbWin_min10SS.pdf"), height=10, width=5)






random_lengthHaplBlocks_all %>%
  group_by(sample) %>%
  summarise(meanfrac=mean(meanfrac)) %>%
  ggplot(aes(meanfrac))+
  geom_histogram(fill="black")+
  geom_vline(data=lengthHaplBlocks_all %>%
               group_by(sample) %>%
               summarise(meanfrac=mean(meanfrac)), aes(xintercept = meanfrac), 
             #linetype="dotted", 
             color = "red", alpha=0.5) +
  labs(x = "Total Length Shared Fraction",
       y = "") +
  # facet_grid(chr~., scale="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))

ggsave(paste0("15_RemCentromeres_RemGaps_totalMeanShareHaplBlocks_perInd_10kbWin_min10SS.png"), height=5, width=5)
ggsave(paste0("15_RemCentromeres_RemGaps_totalMeanShareHaplBlocks_perInd_10kbWin_min10SS.pdf"), height=5, width=5)

```









## Proportion Shared Gaps within Individuals 



```R

# .libPaths("/dss/dsshome1/lxc03/di36guz2/R")
# mamba activate env_others
.libPaths()


setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/05_diversity_withinInd/ROH")
#setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/03_Nhaplotypes")

library(ggplot2)
library(tidyverse)

chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

DM_ref_chrSizes <- data.frame(chr = c("chr01", "chr02", "chr03", "chr04", "chr05", "chr06", "chr07", "chr08", "chr09", "chr10", "chr11", "chr12"), 
                              totalLen = c(88591686, 46102915, 60707570, 69236331, 55599697, 59091578, 57639317, 59226000, 67600300, 61044151, 46777387, 59670755))

chrN <- "01"

for (chrN in chr_list) {
  #print(chrN) }
  #lengthHaplBlocks_all = data.frame(chr="tem", haplotype1="tem", haplotype2="tem", startpos=0, endpos=0)
  table_hapG <- read_table(paste0("../../03_Nhaplotypes/haplotypeGroups_min10SS_chr",chrN,".txt"), F)
  names(table_hapG) <- c("pos", "sample", "hapGroup", "Mean_NAN_haplo")
  
  table_hapG_spread <-  table_hapG %>% 
    filter(Mean_NAN_haplo>0.8) %>%
    select(-Mean_NAN_haplo) %>%
    mutate(hapGroup=1) %>%
    spread(key=sample, value=hapGroup, fill = 0) %>%
    arrange(pos)
  
  list_samples <- unique(table_hapG$sample)
  
  for (hapi in seq(1,length(list_samples), 4)){
    table_hapG_4Hapl <- table_hapG_spread %>%
      select(pos, list_samples[hapi], list_samples[hapi+1], list_samples[hapi+2], list_samples[hapi+3])
    
    gapHapl <- as.vector(((table_hapG_4Hapl[,2] + table_hapG_4Hapl[,3] + table_hapG_4Hapl[,4] + table_hapG_4Hapl[,5])>0)*1)
    
    lengthHaplBlocks <- table_hapG_4Hapl %>%
      mutate(gapHapl) %>%
      #mutate(cont_haplo=(sameHapl==lag(sameHapl) & sameHapl==1)*1) 
      mutate(cont_haplo=(gapHapl + lag(gapHapl) > 0)*1) %>%
      mutate(sequence = data.table::rleid(cont_haplo == 1)) %>% 
      filter(cont_haplo == 1) %>% 
      group_by(sequence) %>%
      summarise(startpos = min(pos), 
                endpos=max(pos)) %>%
      mutate(chr=paste0("chr", chrN), 
             sample=str_split(list_samples[hapi], "_")[[1]][1]) %>%
      select(chr, sample, startpos, endpos)
    write.table(lengthHaplBlocks, file=paste0("GapsTable_SharedBlocksWithinInd_chr", chrN, ".csv"), 
                quote = F, row.names = F, col.names = F, sep = ",", append=T)
  }
}



chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
chrN <- "01"
for (hapi in seq(1,500)){
  print(hapi)
  haplotypesUsed <- sample(1:40, 4, replace=F)
  for (chrN in chr_list) {
    #print(chrN) }
    #lengthHaplBlocks_all = data.frame(chr="tem", haplotype1="tem", haplotype2="tem", startpos=0, endpos=0)
    table_hapG <- read_table(paste0("../../03_Nhaplotypes/haplotypeGroups_min10SS_chr",chrN,".txt"), F)
    names(table_hapG) <- c("pos", "sample", "hapGroup", "Mean_NAN_haplo")
    
    list_samples <- unique(table_hapG$sample)
    
    table_hapG_spread <-  table_hapG %>% 
      filter(Mean_NAN_haplo>0.8) %>%
      select(-Mean_NAN_haplo) %>%
      mutate(hapGroup=1) %>%
      spread(key=sample, value=hapGroup, fill = 0) %>%
      arrange(pos)
    
    table_hapG_4Hapl <- table_hapG_spread %>%
      select(pos, list_samples[haplotypesUsed[1]], list_samples[haplotypesUsed[2]], list_samples[haplotypesUsed[3]], list_samples[haplotypesUsed[4]])
    
    gapHapl <- as.vector(((table_hapG_4Hapl[,2] + table_hapG_4Hapl[,3] + table_hapG_4Hapl[,4] + table_hapG_4Hapl[,5])>0)*1)
    
    lengthHaplBlocks <- table_hapG_4Hapl %>%
      mutate(gapHapl) %>%
      #mutate(cont_haplo=(sameHapl==lag(sameHapl) & sameHapl==1)*1) 
      mutate(cont_haplo=(gapHapl + lag(gapHapl) > 0)*1) %>%
      mutate(sequence = data.table::rleid(cont_haplo == 1)) %>% 
      filter(cont_haplo == 1) %>% 
      group_by(sequence) %>%
      summarise(startpos = min(pos), 
                endpos=max(pos)) %>%
      mutate(chr=paste0("chr", chrN), 
             sample=hapi) %>%
      select(chr, sample, startpos, endpos)
    write.table(lengthHaplBlocks, file=paste0("GapsRandomHap_Table_SharedBlocksWithinInd_chr", chrN, ".csv"), 
                quote = F, row.names = F, col.names = F, sep = ",", append=T)
  }
}



chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
chrN="01"
for (chrN in chr_list) {
  lengthHaplBlocks_all <- read.csv(paste0("GapsTable_SharedBlocksWithinInd_chr",chrN,".csv"), F)
  names(lengthHaplBlocks_all) <- c("chr", "sample", "startpos", "endpos")
  
  lengthHaplBlocks_all %>%
    mutate(lenHaplBlocks=endpos-startpos) %>%
    group_by(chr, sample) %>%
    summarise(NSharedBlocks=n(), 
              totalLenShared=sum(lenHaplBlocks), 
              mean_blockLength=mean(lenHaplBlocks)) %>%
    left_join(DM_ref_chrSizes, by="chr") %>%
    mutate(frac=totalLenShared/totalLen) %>%
    write.table(file=paste0("AllperSample_Gaps_Table_PropSharedBlocks.csv"), 
                quote = F, row.names = F, col.names = F, sep = ",", append=T)
}




chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
chrN="01"
for (chrN in chr_list) {
  lengthHaplBlocks_all <- read.csv(paste0("GapsRandomHap_Table_SharedBlocksWithinInd_chr",chrN,".csv"), F)
  names(lengthHaplBlocks_all) <- c("chr", "sample", "startpos", "endpos")
  
  lengthHaplBlocks_all %>%
    mutate(lenHaplBlocks=endpos-startpos) %>%
    group_by(chr, sample) %>%
    summarise(NSharedBlocks=n(), 
              totalLenShared=sum(lenHaplBlocks), 
              mean_blockLength=mean(lenHaplBlocks)) %>%
    left_join(DM_ref_chrSizes, by="chr") %>%
    mutate(frac=totalLenShared/totalLen) %>%
    write.table(file=paste0("AllRandomHap_Gaps_AllperSample_Table_PropSharedBlocks.csv"), 
                quote = F, row.names = F, col.names = F, sep = ",", append=T)
}


lengthHaplBlocks_all <- read.csv(paste0("AllperSample_Gaps_Table_PropSharedBlocks.csv"), F)
names(lengthHaplBlocks_all) <- c("chr", "sample", "NSharedBlocks", "totalLenShared", "meanblockLength", "chrSize", "frac")
head(lengthHaplBlocks_all)


random_lengthHaplBlocks_all <- read.csv(paste0("AllRandomHap_Gaps_AllperSample_Table_PropSharedBlocks.csv"), F)
names(random_lengthHaplBlocks_all) <- c("chr", "sample", "NSharedBlocks", "totalLenShared", "meanblockLength", "chrSize", "frac")




random_lengthHaplBlocks_all %>%
  ggplot(aes(frac))+
  geom_histogram(fill="black")+
  geom_vline(data=lengthHaplBlocks_all, aes(xintercept = frac), 
             #linetype="dotted", 
             color = "red") +
  labs(x = "Total Length Shared Gaps Fraction",
       y = "") +
  facet_grid(chr~., scale="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))

ggsave(paste0("09_Gaps_MeanShareGapBlocks_perInd_10kbWin_min10SS.png"), height=10, width=5)
ggsave(paste0("09_Gaps_MeanShareGapBlocks_perInd_10kbWin_min10SS.pdf"), height=10, width=5)



random_lengthHaplBlocks_all %>%
  ggplot(aes(totalLenShared))+
  geom_histogram(fill="black")+
  geom_vline(data=lengthHaplBlocks_all, aes(xintercept = totalLenShared), 
             #linetype="dotted", 
             color = "red") +
  labs(x = "Total Length Shared Gap",
       y = "") +
  facet_grid(chr~., scale="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))





random_lengthHaplBlocks_all %>%
  # filter(meanblockLength<300000) %>%
  ggplot(aes(meanblockLength))+
  geom_histogram(fill="black")+
  geom_vline(data=lengthHaplBlocks_all, aes(xintercept = meanblockLength), 
             #linetype="dotted", 
             color = "red") +
  labs(x = "Mean Shared Gap Block Length",
       y = "") +
  facet_grid(chr~., scale="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))

ggsave(paste0("10_Gaps_MeanShareBlockLenght_perInd_10kbWin_min10SS.png"), height=10, width=5)
ggsave(paste0("10_Gaps_MeanShareBlockLenght_perInd_10kbWin_min10SS.pdf"), height=10, width=5)




random_lengthHaplBlocks_all %>%
  ggplot(aes(NSharedBlocks))+
  geom_histogram(fill="black")+
  geom_vline(data=lengthHaplBlocks_all, aes(xintercept = NSharedBlocks), 
             #linetype="dotted", 
             color = "red") +
  labs(x = "Mean Number Shared Gap Blocks",
       y = "") +
  facet_grid(chr~., scale="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))

ggsave(paste0("11_Gaps_MeanShareBlockLenght_perInd_10kbWin_min10SS.png"), height=10, width=5)
ggsave(paste0("11_Gaps_MeanShareBlockLenght_perInd_10kbWin_min10SS.pdf"), height=10, width=5)



```



# LD



## Excluding singletons



### Remove singletons and subsampling genotype table:



First the genotype table was edited to: 1. alleles were changed to "1" the most common allele and "0" other alleles. Missing alleles were maintained as "nan". 2. Singleton variants were removed. 3. only one variant site was included every X distance (subsampling of variants). Distance was defined by the parameter "min_dis". In this case 500 bp was used. 4. At least there should be 4 genotypes (4 haplotypes without nan) in a variant site in order to consider it. 



python script: editGenotypeTable_subsampling_removesingletons.py

```python
import numpy as np
import sys
import os

file="test_AllSamplesMultiAlleles_chr01_missing_genotypeTable.txt"
output="LDtest2"
min_dis = int(500)

file = str(sys.argv[1])
output=str(sys.argv[2])
min_dis = int(sys.argv[3])


def most_frequent(List):
    unique, counts = np.unique(List, return_counts=True)
    index = np.argmax(counts)
    return unique[index]

def chage_genotype(genotype_list, allele):
	return ['1' if b==allele else '0' if b!='nan' else 'nan' for b in genotype_list]


#output file 1: LD output:
output_table = open(f'{output}.txt', "w")

start_current = 0

for line in open(file, "r"):
	if line[0:3] == "CHR":
		line_values  = line.strip().split("\t")
		line_join = '\t'.join(map(str, line_values[0:3]+line_values[5:]))
		output_table.write(f'{line_join}\n')
	else:
		line_values  = line.strip().split("\t")
		chr, start, end, ref, alt = line_values[0:5]
		if int(start) - start_current > min_dis:
			genotypes = [x for x in line_values[5:]]
			Non_nan_genotypes = [x != 'nan' for x in genotypes] 
			index_Non_nan_genotypes = [x for x in range(len(Non_nan_genotypes)) if Non_nan_genotypes[x]] 
			filtered_genotypes = [genotypes[i] for i in index_Non_nan_genotypes]
			if len(filtered_genotypes) > 5:
				allele_1_genotypes = most_frequent(filtered_genotypes)
				alleles_genotypes = [(x == allele_1_genotypes)*1 for x in filtered_genotypes]
				total = len(alleles_genotypes)
				Np1=sum(alleles_genotypes)
				if (Np1 < total-1) and (total > 10):
					new_genotypes = chage_genotype(genotypes, allele_1_genotypes)
					line_join = '\t'.join(map(str, line_values[0:3]+new_genotypes))
					output_table.write(f'{line_join}\n')
					start_current = int(start)


output_table.close()

```



bash script: editGenotypeTable_subsampling_removesingletons.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=1-00:00:00
#SBATCH -J subSamplingGT

file=$1
output=$2
min_dis=$3

python3 editGenotypeTable_subsampling_removesingletons.py $file $output $min_dis 

```

job:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/04_LD

mamba activate mypython3

# /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr01_missing_genotypeTable.txt

minDis=500
for i in $(seq -w 1 12)
do chr="chr"$i ; echo $chr
sbatch editGenotypeTable_subsampling_removesingletons.sh /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr$i"_"missing_genotypeTable.txt AllSamplesMultiAlleles_chr$i"_"missing_genotypeTable_subSampling$minDis"bp"_NoSingletons $minDis
done
# 05.04.2024: 2456111..2456122


```



### Calculate LD:

Python script: calculate_LD_withoutSingletons.py

```python
import numpy as np
import sys
import os

file="test_AllSamplesMultiAlleles_chr01_missing_genotypeTable.txt"
file="AllSamplesMultiAlleles_chr01_missing_genotypeTable_subSampling500bp_NoSingletons.txt"
output="LDtest2"
start_current = int(0)
end_pos = int(2000)
max_distance = 20000000

file = str(sys.argv[1])
output=str(sys.argv[2])
start_current = int(sys.argv[3])
end_pos = int(sys.argv[4])
max_distance = int(sys.argv[5])


#output file 1: LD output:
output_popParameters = open(f'{output}.txt', "w")

for line in open(file, "r"):
	if line[0:3] == "CHR":
		#output_popParameters.write(f'chr\tstart\tstart2\ttotal\tp1\tq1\tp1q1\tD\tDp\tr2\n')
		output_popParameters.write(f'chr\tstart\tstart2\ttotal\tp1\tq1\tp1q1\tD\tr2\n')
	else:
		line_values  = line.strip().split("\t")
		chr, start, end = line_values[0:3]
		if int(start) > end_pos:
			break
		if int(start) > start_current:
			genotypes = [x for x in line_values[3:]]
			Non_nan_genotypes = [x != 'nan' for x in genotypes] 
			index_Non_nan_genotypes = [x for x in range(len(genotypes)) if Non_nan_genotypes[x]] 
			start_current = int(start)
			start2_current = 0
			for line2 in open(file, "r"):
				if line2[0:3] == "CHR":
					continue
				else:
					line_values  = line2.strip().split("\t")
					chr2, start2, end2 = line_values[0:3]
					if (int(start2) > int(start)) and (int(start2) > start2_current):
						genotypes2 = [x for x in line_values[3:]]
						Non_nan_genotypes2 = [x != 'nan' for x in genotypes2] 
						index_Non_nan_genotypes2 = [x for x in range(len(genotypes2)) if Non_nan_genotypes2[x]] 
						total_non_na = list(set(index_Non_nan_genotypes).intersection(set(index_Non_nan_genotypes2)))
						filterd_genotypes = [int(genotypes[i]) for i in total_non_na] 
						filterd_genotypes2 = [int(genotypes2[i]) for i in total_non_na]
						total = len(filterd_genotypes)
						if total > 0:
							Np1=sum(filterd_genotypes)
							Nq1=sum(filterd_genotypes2)
							p1=Np1/total
							q1=Nq1/total
							if (0 < p1 < 1) and (0 < q1 < 1):
								start2_current = int(start2)
								p1q1 = sum([(sum(i) == 2)*1 for i in zip(filterd_genotypes, filterd_genotypes2)])/total
								D=p1q1-(p1*q1)
								#if D<0:
								#	Dp = D/(max([(-p1*q1),(-(1-p1)*(1-q1))]))
								#else:
								#	Dp =  D/(min([(p1*(1-q1)),((1-p1)*(q1))])) 
								r2=(D**2)/(p1*(1-p1)*q1*(1-q1))
								output_popParameters.write(f'{chr}\t{start}\t{start2}\t{total}\t{p1}\t{q1}\t{p1q1}\t{D}\t{r2}\n')
								if (start2_current - start_current > max_distance):
									break


output_popParameters.close()


```



Bash script job: calculate_LD_withoutSingletons.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=8-00:00:00
#SBATCH -J LD_cal

file=$1
output=$2
start_current=$3
end_pos=$4
max_distance=$5

python3 calculate_LD_withoutSingletons.py $file $output $start_current $end_pos $max_distance

```



Job: 

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/04_LD

mamba activate mypython3

# /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr01_missing_genotypeTable.txt

grep -v scaffold_ ../../../02_reference/DM_1-3_516_R44_potato_genome_assembly.v6.1/scatter.interval_list | grep "@SQ" | awk '{print $2"\t"$3}' | sed 's@SN:@@g' | sed 's@LN:@@g' > chr_size.txt

# each chromosome was split into 10 to speed up the job:
max_distance=30000000
max_distance_name=30Mb
for i in $(seq -w 1 12)
do chr="chr"$i ; echo $chr
mkdir LD_NoSingletons_$chr"_maxDis"$max_distance_name
chrSize=$(grep $chr chr_size.txt | awk '{print $2}') ; echo $chrSize
splitsize=$(echo $chrSize | awk '{printf "%.0f\n", $1/10}') ; echo $splitsize
for startpos in $(seq -w 0 $splitsize $chrSize)
do finalpos=$(echo $startpos | awk -v Ssize=$splitsize '{print $1+Ssize}')
echo $startpos $finalpos
sbatch calculate_LD_withoutSingletons.sh AllSamplesMultiAlleles_chr$i"_"missing_genotypeTable_subSampling500bp_NoSingletons.txt ./LD_NoSingletons_$chr"_maxDis"$max_distance_name"/LD_chr"$i"_Start"$startpos"_End"$finalpos"_maxDis"$max_distance_name $startpos $finalpos $max_distance
done
done
# old without non-aligned regions: 2156404..2156528






#  Version using the whole chromosome distance:
max_distance=100000000
max_distance_name=100Mb
for i in $(seq -w 1 12)
do chr="chr"$i ; echo $chr
mkdir LD_NoSingletons_$chr"_maxDis"$max_distance_name
chrSize=$(grep $chr chr_size.txt | awk '{print $2}') ; echo $chrSize
splitsize=$(echo $chrSize | awk '{printf "%.0f\n", $1/10}') ; echo $splitsize
for startpos in $(seq -w 0 $splitsize $chrSize)
do finalpos=$(echo $startpos | awk -v Ssize=$splitsize '{print $1+Ssize}')
echo $startpos $finalpos
sbatch calculate_LD_withoutSingletons.sh AllSamplesMultiAlleles_chr$i"_"missing_genotypeTable_subSampling500bp_NoSingletons.txt ./LD_NoSingletons_$chr"_maxDis"$max_distance_name"/LD_chr"$i"_Start"$startpos"_End"$finalpos"_maxDis"$max_distance_name $startpos $finalpos $max_distance
done
done
# 05.04.2024: 2456123..2456247

# for i in {2183956..2184011}; do stop $i ; rm slurm-$i* ; done

```



sampleLines.py

```python
#!/bin/python3
import os
import sys
import random
import argparse

parser = argparse.ArgumentParser(description="Returns random lines from text files, as fast as possible.")
parser.add_argument('--file', required=True, help="the file to parse")
parser.add_argument('-n', required=True, help="amount of lines you want to return")
parser.add_argument('--min-len', dest='minLen', type=int, default=0, help="minimum line length, defaults to 0")
parser.add_argument('--max-len', dest='maxLen', type=int, default=float('inf'), help="maximum line length, defaults to ∞")
args = parser.parse_args()

# len mistake catch
if args.minLen > args.maxLen:
    sys.exit("I can't think of any strings that are longer than " + str(args.minLen) + " while being shorter than " + str(args.maxLen) + ".")

filename = args.file
filesizeBytes = os.path.getsize(filename)
#print(filesizeBytes)

# you can find this with: sudo blockdev --getbsz $(df --output=source $(pwd) | grep '/')
blockSize = 4096               # Bytes

import random
f = open(filename, "rt")
count = 0
bytesRead = 0
while (count < 10000):
  buffer = f.read(blockSize)
  if not buffer: break
  count += buffer.count('\n')
  bytesRead += blockSize
f.close()

# less then 1 blocksize file fix
if bytesRead > filesizeBytes:
    bytesRead = filesizeBytes

bytesPerLine = bytesRead/count
#print(bytesPerLine)

totalLinesEst = filesizeBytes / bytesPerLine
#print("Estimated lines in file: " + str(totalLinesEst))

# exact output count
resultlines = 0

fh = open(filename)
while resultlines < int(args.n):
    linestart = random.randint(0, int(totalLinesEst))
    readstart = (linestart * bytesPerLine) - bytesPerLine
    if readstart < bytesPerLine:
        readstart = 0
    else:
        fh.seek(readstart)
    scratch = fh.readline()
    line = fh.readline().replace('\n','')
    while not (args.minLen <= len(line) <= args.maxLen):
        line = fh.readline().replace('\n','')
    print(line)
    resultlines += 1
fh.close()

```



**SubSampling LD table:**

bash script: subsampling_LDtable.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=8-00:00:00
#SBATCH -J SubSample_LD


max_distance=$1
max_distance_name=$2
sampling_lines=$3
chr="chr"$4

#for i in $(seq -w 1 12)
#i=01
#chr="chr"$i
#chrSize=$(grep $chr chr_size.txt | awk '{print $2}') ; echo $chrSize
#ProSampling=$(echo $chrSize | awk -v SSize=$sampling_lines -v maxD=/$max_distance '{print SSize/(((($1/500)-1)*($1/500))/2)}') ; echo $ProSampling
# cat ./LD_NoSingletons_$chr"_maxDis"$max_distance_name"/LD_"$chr"_"* | awk -v PROB=$ProSampling -v CHR=$chr -v distanceN=$max_distance_name 'BEGIN  {srand()} !/^$/ { if (rand() <= PROB || FNR==1) print > "zLD_"CHR"_NoSingletons_"distanceN".txt"}' 
cat ./LD_NoSingletons_$chr"_maxDis"$max_distance_name"/LD_"$chr"_"* > total_LD_"$chr"_"maxDis"$max_distance_name.subSample$sampling_lines.txt

python3 sampleLines.py -n $sampling_lines --file total_LD_"$chr"_"maxDis"$max_distance_name.subSample$sampling_lines.txt > all_LD_"$chr"_"maxDis"$max_distance_name.subSample$sampling_lines.txt

awk '{print $1"\t"$3-$2"\t"$9}' all_LD_"$chr"_"maxDis"$max_distance_name.subSample$sampling_lines.txt > dis_all_LD_"$chr"_"maxDis"$max_distance_name.subSample$sampling_lines.txt

rm total_LD_"$chr"_"maxDis"$max_distance_name.subSample$sampling_lines.txt
#done



```

Job:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/04_LD

mamba activate mypython3

sbatch subsampling_LDtable.sh 30000000 30Mb 1000000
# old: 2157149


#max_distance=$1
#max_distance_name=$2
#sampling_lines=$3
#chr="chr"$4

for i in $(seq -w 1 12)
do echo $i
if [ $i -lt 6 ]
then
  echo $i
  sbatch tem_subsampling_LDtable.sh 100000000 100Mb 3000000 $i
fi
done
# 30.11.2023: 2253639..2253643 (1 to 5)

for i in $(seq -w 1 12)
do echo $i
sbatch subsampling_LDtable.sh 100000000 100Mb 3000000 $i
done
# 07.04.2024: 2457530..2457541


# for i in {2253634..2253638}; do stop $i ; rm slu*""$i.out	; done

```



### Produce mean values per 50kb window:

python scrip: meanLDvalues_subsampling.py

```python
import os
import sys
import glob
import random
import math
import numpy as np



chrID = str("chr01")
win_size = int(10000)
chr_size = int(90000000)
output_name = str("test")

chrID = str(sys.argv[1])
win_size = int(sys.argv[2])
chr_size = int(sys.argv[3])
output_name = str(sys.argv[4])

FilenamesList = glob.glob('LD_NoSingletons_'+chrID+'_maxDis100Mb/LD_'+chrID+'_*.txt')

NWin = math.floor(int(chr_size) / win_size)+5
sumR2 = np.zeros((NWin, NWin))
NLines = np.zeros((NWin, NWin))

output_subsampling = open(f'Subsampling_python_{chrID}_maxDis100Mb_WinSize{win_size}.txt', "w")
output_meanVal = open(f'{output_name}_{chrID}_maxDis100Mb_WinSize{win_size}.txt', "w")

for file in FilenamesList:
    for line in open(file, "r"):
        if line[0:5] == "chr\ts":
            continue
        else:
            if random.random() < 0.0005:
                output_subsampling.write(f'{line}')
            line_values = line.strip().split("\t")
            win1 = math.floor(_ / win_size)
            win2 = math.floor(int(line_values[2]) / win_size)
            sumR2[win1, win2] += float(line_values[8])
            NLines[win1, win2] += 1


for win1 in range(NWin):
    for win2 in range(NWin):
        if win2 >= win1:
            if NLines[win1, win2] > 0:
                output_meanVal.write(f'{win1*win_size}\t{win2*win_size}\t{sumR2[win1, win2]/NLines[win1, win2]}\t{NLines[win1, win2]}\n')



output_subsampling.close()
output_meanVal.close()


```



bash script: meanLDvalues_subsampling.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00
#SBATCH -J meanLD_cal

chrID=$1
win_size=$2
chr_size=$3
output_name=$4

python3 meanLDvalues_subsampling.py $chrID $win_size $chr_size $output_name


```

Job:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/04_LD

mamba activate mypython3

winSize=10000
for i in $(seq -w 1 12)
do chr="chr"$i ; echo $chr
chrSize=$(grep $chr chr_size.txt | awk '{print $2}') ; echo $chrSize
sbatch meanLDvalues_subsampling.sh $chr $winSize $chrSize meanLDValue
done
# old: 98..2157509


winSize=30000
for i in $(seq -w 1 12)
do chr="chr"$i ; echo $chr
chrSize=$(grep $chr chr_size.txt | awk '{print $2}') ; echo $chrSize
sbatch meanLDvalues_subsampling.sh $chr $winSize $chrSize MeanLDValue
done
# old: 2157872..2157883

winSize=50000
for i in $(seq -w 1 12)
do chr="chr"$i ; echo $chr
chrSize=$(grep $chr chr_size.txt | awk '{print $2}') ; echo $chrSize
sbatch meanLDvalues_subsampling.sh $chr $winSize $chrSize MeanLDValue
done
# 07.04.2024: 2457542..2457553


winSize=100000
for i in $(seq -w 1 12)
do chr="chr"$i ; echo $chr
chrSize=$(grep $chr chr_size.txt | awk '{print $2}') ; echo $chrSize
sbatch meanLDvalues_subsampling.sh $chr $winSize $chrSize MeanLDValue
done
#  07.04.2024: 2457554..2457565

winSize=1000000
for i in $(seq -w 1 12)
do chr="chr"$i ; echo $chr
chrSize=$(grep $chr chr_size.txt | awk '{print $2}') ; echo $chrSize
sbatch ./meanLDvalues_subsampling.sh $chr $winSize $chrSize MeanLDValue
done
#  20.06.2024: 2594112..2594123

# for i in {2159042..2159053}; do stop $i ; rm slurm-$i* ; done

for i in $(seq -w 1 12)
do chr="chr"$i ; echo $chr
awk '{print $1"\t"$3-$2"\t"$9}' Subsampling_python_$chr"_"maxDis100Mb_WinSize50000.txt > dis_Subsampling_python_$chr"_"maxDis100Mb_WinSize50000.txt
done
#done

```



#### Plot:

```R

# .libPaths("/dss/dsshome1/lxc03/di36guz2/R")
# mamba activate env_others
.libPaths()

setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/04_LD")

library(ggplot2)
library(tidyverse)

# centromere_table <- read.csv("centromeres.csv", header = T)

chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

chrID <- "01"

WinSize <- 50000
for (chrID in chr_list) {
  
  WinSize <- 50000
  
  table_LD <- read_table(paste0("MeanLDValue_chr",chrID,"_maxDis100Mb_WinSize",WinSize,".txt"), F)
  names(table_LD) <- c("pos1", "pos2", "meanLD", "NVar")
  
  table_LD %>%
    filter(pos1 != pos2) %>%
    ggplot(aes(pos1/1000000, pos2/1000000, fill=meanLD, colour=meanLD)) +
    geom_tile()+
    #scico::scale_fill_scico(palette = "lajolla")+
    scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, 1))+
    scale_colour_viridis_c(option = "magma", direction = -1, limits = c(0, 1), guide = 'none')+
    scale_x_continuous(breaks=seq(0,100,5), name = "Position (Mb)", expand = c(0,0))+
    scale_y_continuous(breaks=seq(0,100,5), name = "Position (Mb)", expand = c(0,0))+
    # scale_y_discrete(name = "Haplotype", expand = c(0,0))+
    ggtitle(paste0("Chr:", chrID))+
    labs(fill="LD (R^2)")+
    theme_classic()+
    theme(panel.grid.major = element_line(color = "grey80"),
          panel.grid.minor = element_line(color = "grey80"))
  
  ggsave(paste0("02_heatmapLD_chr", chrID, "_WinSize", WinSize, ".png"), height=10, width=12)
  
  WinSize <- 200000
  table_LD %>%
    filter(pos1 != pos2) %>%
    mutate(pos1= floor(pos1/WinSize)*WinSize,
           pos2= floor(pos2/WinSize)*WinSize) %>%
    rename(rawmeanLD=meanLD) %>%
    group_by(pos1, pos2) %>%
    # summarise(meanLD=(mean(rawmeanLD)*(0.6-0.1))+0.2) %>%
    # summarise(meanLD=if_else(meanLD>0.6, 0.6,
    #                          if_else(meanLD<0.1, 0.1, meanLD))) %>%
    summarise(meanLD=if_else(mean(rawmeanLD)>0.5, 0.5,
                             if_else(mean(rawmeanLD)<0.1, 0.1, mean(rawmeanLD)))) %>%
    ggplot(aes(pos2/1000000, pos1/1000000, fill=meanLD, colour=meanLD)) +
    geom_tile()+
    #scico::scale_fill_scico(palette = "lajolla")+
    scale_fill_gradient(low="white", high="darkgreen") +
    scale_colour_gradient(low="white", high="darkgreen", guide='none') +
    # scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, 1))+
    # scale_colour_viridis_c(option = "magma", direction = -1, limits = c(0, 1), guide = 'none')+
    scale_x_continuous(breaks=seq(0,100,5), name = "Position (Mb)", expand = c(0,0))+
    scale_y_continuous(breaks=seq(0,100,5), name = "Position (Mb)", expand = c(0,0), position = "right")+
    # scale_y_discrete(name = "Haplotype", expand = c(0,0))+
    labs(fill="LD (R^2)")+
    theme_classic()+
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "none")
  
  ggsave(paste0("02_heatmapLD_Simple_chr", chrID, "_WinSize", WinSize, ".png"), height=12, width=12)
}

# just to produce the scale:
WinSize <- 200000
table_LD %>%
  filter(pos1 != pos2) %>%
  mutate(pos1= floor(pos1/WinSize)*WinSize,
         pos2= floor(pos2/WinSize)*WinSize) %>%
  rename(rawmeanLD=meanLD) %>%
  group_by(pos1, pos2) %>%
  summarise(meanLD=if_else(rawmeanLD>0.6, 0.6,
                           if_else(rawmeanLD<0.1, 0.1, rawmeanLD))) %>%
  # summarise(meanLD=mean(rawmeanLD)) %>%
  ggplot(aes(pos2/1000000, pos1/1000000, fill=meanLD, colour=meanLD)) +
  geom_tile()+
  #scico::scale_fill_scico(palette = "lajolla")+
  scale_fill_gradient(low="white", high="darkgreen", limits = c(0, 1)) +
  scale_colour_gradient(low="white", high="darkgreen", limits = c(0, 1), guide='none') +
  # scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, 1))+
  # scale_colour_viridis_c(option = "magma", direction = -1, limits = c(0, 1), guide = 'none')+
  scale_x_continuous(breaks=seq(0,100,5), name = "Position (Mb)", expand = c(0,0))+
  scale_y_continuous(breaks=seq(0,100,5), name = "Position (Mb)", expand = c(0,0), 
                     position = "right")+
  # scale_y_discrete(name = "Haplotype", expand = c(0,0))+
  labs(fill="LD (R^2)")+
  theme_classic()+
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "bottom")

ggsave(paste0("02_heatmapLD_Simple_chr", chrID, "_WinSize", WinSize, "_withLegend.png"), height=12, width=12)




```



### Linkage within 500kb distance



```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/04_LD

for i in $(seq -w 1 12)
do chr="chr"$i ; echo $chr
awk -v CHR=$chr '{if($2-$1<1500000) print CHR"\t"$0}' MeanLDValue_$chr"_"maxDis100Mb_WinSize50000.txt > within1.5Mb_MeanLDValue_$chr"_"maxDis100Mb_WinSize50000.txt
sort within1.5Mb_MeanLDValue_$chr"_"maxDis100Mb_WinSize50000.txt | uniq > ed_within1.5Mb_MeanLDValue_$chr"_"maxDis100Mb_WinSize50000.txt
done


cat ed_within1.5Mb_MeanLDValue_ch*maxDis100Mb_WinSize50000.txt > ed_within1.5Mb_MeanLDValue_all_maxDis100Mb_WinSize50000.txt

# done

```



#### Plot:



```R


# .libPaths("/dss/dsshome1/lxc03/di36guz2/R")
# mamba activate env_others
.libPaths()

setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/04_LD")

library(ggplot2)
library(tidyverse)

# centromere_table <- read.csv("centromeres.csv", header = T)

chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

chrID <- "01"

# plot within 500kb
LD_table <- read_table(paste0("ed_within1.5Mb_MeanLDValue_all_maxDis100Mb_WinSize50000.txt"), F)
names(LD_table) <- c("chr", "pos1", "pos2", "r2", "Nvar")

winSize=500000
meanLD_table <- LD_table %>%
  mutate(pos=floor(pos1/winSize)*winSize) %>%
  mutate(pos2w=floor(pos2/winSize)*winSize) %>%
  filter(pos==pos2w) %>%
  group_by(chr, pos) %>%
  dplyr::summarise(mean_r2=mean(r2)) 

quantile_95 <- meanLD_table %>% 
  pull(mean_r2) %>%
  quantile(0.93)

quantile_95

meanLD_table %>%
  ggplot(aes(pos/1000000, mean_r2))+
  geom_ribbon(aes(ymax=ifelse(mean_r2>mean(meanLD_table$mean_r2), mean_r2, mean(meanLD_table$mean_r2)), 
                  ymin=mean(meanLD_table$mean_r2)), fill="#91b873") +
  geom_line(alpha=0.8)+
  geom_point(data=meanLD_table %>%
               filter(mean_r2>quantile_95) %>%
               mutate(posy=(0.1)), aes(pos/1000000, posy), colour="#203f08") +
  # geom_hline(yintercept=0.1, linetype="dashed", color = "red")+
  geom_line(data=centromere_table %>%
              gather(key="posB", value = "pos", -chr) %>%
              mutate(posy = 0), 
            aes(pos/1000000,posy), colour="darkgreen", linewidth=1.5, alpha=0.5) +
  #geom_point(alpha=0.04)+
  geom_hline(yintercept=mean(meanLD_table$mean_r2), linetype="dashed", color = "red")+
  scale_y_continuous(breaks=seq(0,1,0.1), limits = c(0,1)) +
  scale_x_continuous(breaks=seq(0,100,15)) +
  facet_grid(. ~ chr, scale="free", space="free")+
  labs(x = "Pos. (Mb)",
       y = "mean LD (R^2)") +
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12), 
        axis.title = element_text(size=14))


ggsave("03_MeanLD_within500kbWin.png", height=2.5, width=15)
ggsave("03_MeanLD_within500kbWin.pdf", height=2.5, width=15)



```



### Linkage decay plot:



```R
# .libPaths("/dss/dsshome1/lxc03/di36guz2/R")
# mamba activate env_others
.libPaths()

setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/04_LD")

library(ggplot2)
library(tidyverse)
library(ggrepel)
# library(geomtextpath)


out <- read_table(paste0("dis_Subsampling_python_chr03_maxDis100Mb_WinSize50000.txt"), F)
names(out) <- c("chr", "Distance", "Pearson_r2")


out %>%
  ggplot(aes(Distance)) +
  geom_histogram()

out %>%
  ggplot(aes(Distance/1000000, Pearson_r2)) +
  geom_point(alpha=0.01)

# The decay of r2 with distance was fitted using Hill and Weir expectation of r2 between adjacent sites (Hill & Weir 1988). In accordance with previous work (Remington et al. 2001), we used the equation
# 	LD.data~((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance)))
# 	where n is the sample size and C, the parameter to be estimated, represents the product of the population recombination parameter (rho=4Ne*r) and the distance in base pairs.


ch<-"01"
chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

# sampling <- 3000000

for (ch in chr_list){
  # out <- read_table(paste0("dis_all_LD_chr", ch, "_maxDis100Mb.subSample.txt"), F)
  out <- read_table(paste0("dis_Subsampling_python_chr", ch, "_maxDis100Mb_WinSize50000.txt"), F)
  names(out) <- c("chr", "Distance", "Pearson_r2")
  
  out1 <- out[order(out$Distance),]
  distance <- out1$Distance
  LD.data <- out1$Pearson_r2
  n <- 40
  HW.st <- c(C=0.1)
  HW.nonlinear <- nls(LD.data~((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),start=HW.st,control=nls.control(maxiter=500))
  tt <- summary(HW.nonlinear)
  new.rho <- tt$parameters[1]			# 
  # out$fpoints <-((10+new.rho*distance)/((2+new.rho*distance)*(11+new.rho*distance)))*(1+((3+new.rho*distance)*(12+12*new.rho*distance+(new.rho*distance)^2))/(n*(2+new.rho*distance)*(11+new.rho*distance)))
  # HalfDecay <- max(LD.data)/2
  # down <- out1[max(which(fpoints>=HalfDecay)), ]	#
  # up <- out1[min(which(fpoints<HalfDecay)), ]		#
  
  tem_predicted_regression <- data.frame(chr=paste0("chr", ch), 
                                         distance_seq=seq(floor(min(out$Distance)/1000)*1000, floor(max(out$Distance)/1000)*1000, 1000)) %>%
    mutate(fpoints=((10+new.rho*distance_seq)/((2+new.rho*distance_seq)*(11+new.rho*distance_seq)))*(1+((3+new.rho*distance_seq)*(12+12*new.rho*distance_seq+(new.rho*distance_seq)^2))/(n*(2+new.rho*distance_seq)*(11+new.rho*distance_seq))))
  
  dis_LD50=tem_predicted_regression$distance_seq[which(abs(0.5-tem_predicted_regression$fpoints) == min(abs(0.5-tem_predicted_regression$fpoints)))]
  dis_LD25=tem_predicted_regression$distance_seq[which(abs(0.25-tem_predicted_regression$fpoints) == min(abs(0.25-tem_predicted_regression$fpoints)))]
  dis_LD1=tem_predicted_regression$distance_seq[which(abs(0.1-tem_predicted_regression$fpoints) == min(abs(0.1-tem_predicted_regression$fpoints)))]
  
  if (ch=="01"){
    
    predicted_regression <- tem_predicted_regression
    
    data.frame(chr=paste0("chr", ch), 
               rho=new.rho, 
               disLD50=dis_LD50, 
               disLD25=dis_LD25, 
               disLD1=dis_LD1) %>%
      write.table(file=paste0("./summary_values_perChr_subSampling.csv"), sep = ",", row.names = F, col.names = T)
    
  } else {
    predicted_regression <- rbind(predicted_regression, 
                                  tem_predicted_regression)
    
    data.frame(chr=paste0("chr", ch), 
               rho=new.rho, 
               disLD50=dis_LD50, 
               disLD25=dis_LD25, 
               disLD1=dis_LD1) %>%
      write.table(file=paste0("./summary_values_perChr_subSampling.csv"), sep = ",", row.names = F, col.names = F, append = T)
    
  }
}

# write.table(predicted_regression, file=paste0("./regressionValues_subSampling.csv"), sep = ",", row.names = F, col.names = T)

predicted_regression <- read_csv(paste0("./regressionValues_subSampling.csv"), T)


predicted_regression_labels <- predicted_regression %>%
  filter(distance_seq==1000000) 


ggplot(predicted_regression, aes(log10(distance_seq), fpoints, group=chr)) +
  geom_line()+
  # geom_textline(hjust = .7)+
  geom_text_repel(aes(label = chr), data = predicted_regression_labels,
                  nudge_x = -1.5,
                  nudge_y = -0.04,
                  color = "grey30", 
                  size = 5)+
  labs(x = "Distance (log10)",
       y = "R^2") +
  # facet_grid(chr ~ ., scale="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor = element_line(color = "grey80"),
        axis.text = element_text(size=24), 
        axis.title = element_text(size=24))

ggsave(paste0("01_decay_LD_subSampling.png"), height=8, width=10)
ggsave(paste0("01_decay_LD_subSampling.pdf"), height=8, width=10)


summary_values_perChr <- read_csv(paste0("./summary_values_perChr_subSampling.csv"), T)
summary_values_perChr

arrange(summary_values_perChr, disLD25) %>%
  select(-disLD50) 
# 


```



Summary table with a subsampling of 3 million lines per chromosome:



| Chr   | rho         | disLD<0.25 (bp) | disLD<0.1 (bp) |
| ----- | ----------- | --------------- | -------------- |
| chr02 | 0.00000259  | 794000          | 4062000        |
| chr03 | 0.00000182  | 1126000         | 5761000        |
| chr12 | 0.00000177  | 1163000         | 5949000        |
| chr11 | 0.00000166  | 1238000         | 6337000        |
| chr07 | 0.00000158  | 1299000         | 6648000        |
| chr04 | 0.00000146  | 1403000         | 7177000        |
| chr06 | 0.00000136  | 1510000         | 7728000        |
| chr01 | 0.00000120  | 1712000         | 8758000        |
| chr08 | 0.00000109  | 1876000         | 9599000        |
| chr09 | 0.000000960 | 2141000         | 10953000       |
| chr05 | 0.000000885 | 2320000         | 11872000       |
| chr10 | 0.000000543 | 3782000         | 19353000       |

New version 08.04.2024



```
1 chr02 0.00000184  1115000  5708000
 2 chr03 0.00000139  1479000  7567000
 3 chr06 0.00000112  1830000  9366000
 4 chr11 0.00000107  1919000  9822000
 5 chr01 0.000000882 2328000 11913000
 6 chr04 0.000000860 2389000 12222000
 7 chr08 0.000000728 2821000 14437000
 8 chr07 0.000000559 3673000 18795000
 9 chr09 0.000000521 3939000 20158000
10 chr05 0.000000427 4815000 24641000
11 chr12 0.000000391 5250000 26862000
12 chr10 0.000000384 5350000 27375000
```



# Recombination rate

Linkage map data were taken from Guihermeda Silva Perira 2021 

Script: recombinationRate.py

```python
import pandas as pd
import numpy as np
import sys

file = str(sys.argv[1])
output = str(sys.argv[2])
winSize = int(sys.argv[3])

# Read the data from the text file
df1 = pd.read_csv(file, sep='\t')

# Define window size
window_size = winSize

# Initialize an empty list for recombination rates per window
window_data = []

# Iterate over each chromosome
for chrom, chrom_group in df1.groupby('Chrom'):
    # Calculate recombination rate (centimorgan/Mbp) for each interval within the chromosome
    chrom_group['RecombinationRate'] = chrom_group['MapPosition'].diff() / (chrom_group['Position'].diff() / 1e6)
    # Remove the first row of each chromosome because it will have NaN for recombination rate
    chrom_group = chrom_group.dropna(subset=['RecombinationRate'])
    # Create a function to map position to window
    def position_to_window(position, window_size):
        return int((position // window_size) + 1)
    # Iterate over each row in the chromosome group to fill recombination rates per window
    for i in range(1, len(chrom_group)):
        start_pos = chrom_group.iloc[i-1]['Position']
        end_pos = chrom_group.iloc[i]['Position']
        rate = chrom_group.iloc[i]['RecombinationRate']
        start_window = position_to_window(start_pos, window_size)
        end_window = position_to_window(end_pos, window_size)
        for window in range(start_window, end_window + 1):
            window_start_pos = (window - 1) * window_size
            window_end_pos = window * window_size
            overlap_start = max(start_pos, window_start_pos)
            overlap_end = min(end_pos, window_end_pos)
            # Calculate the proportion of the interval that falls within the current window
            proportion = (overlap_end - overlap_start) / (end_pos - start_pos)
            #window_rate = rate * proportion
            size_bp = (end_pos - start_pos) * proportion
            window_data.append([int(chrom), window, size_bp, rate])

# Convert the list to a dataframe
recombination_per_window = pd.DataFrame(window_data, columns=['Chrom', 'Window', 'size_bp', 'RecombinationRate'])

recombination_per_window['GenDistance'] = (recombination_per_window['size_bp']*recombination_per_window['RecombinationRate'])/winSize

# Group by Chrom and Window and calculate the mean recombination rate per window
recombination_per_window = recombination_per_window[['Chrom', 'Window', 'GenDistance']].groupby(['Chrom', 'Window'], as_index=False).sum()

# Ensure the 'Chrom' column is of integer type
recombination_per_window['Chrom'] = recombination_per_window['Chrom'].astype(int)

# Save the output to a new text file
recombination_per_window.to_csv(output, sep='\t', index=False)

```



Job:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/08_recombinationRateCorrelations

mamba activate env_others

python3 recombinationRate.py linkageMap_data.txt rec_per10KbWin.txt 10000
python3 recombinationRate.py linkageMap_data.txt rec_per100KbWin.txt 100000
python3 recombinationRate.py linkageMap_data.txt rec_per500KbWin.txt 500000
python3 recombinationRate.py linkageMap_data.txt rec_per1MbWin.txt 1000000


```



# Distribution of SV along reference genome



script: extractSeq.py

```python
import argparse

def parse_fasta(filename):
    sequences = {}
    with open(filename, 'r') as file:
        seq_name = ''
        seq_data = ''
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if seq_name:
                    sequences[seq_name] = seq_data
                seq_name = line[1:]
                seq_data = ''
            else:
                seq_data += line
        if seq_name:
            sequences[seq_name] = seq_data
    return sequences

def write_fasta(filename, seq_name, seq_data):
    with open(filename, 'w') as file:
        file.write(f'>{seq_name}\n')
        for i in range(0, len(seq_data), 80):
            file.write(seq_data[i:i+80] + '\n')

def main(input_file, output_file, old_name, new_name):
    sequences = parse_fasta(input_file)
    if old_name not in sequences:
        print(f'Sequence with name {old_name} not found in {input_file}')
        return
    write_fasta(output_file, new_name, sequences[old_name])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract and rename a sequence from a FASTA file.')
    parser.add_argument('input_file', help='Input FASTA file')
    parser.add_argument('output_file', help='Output FASTA file')
    parser.add_argument('old_name', help='Name of the sequence to extract from input file')
    parser.add_argument('new_name', help='Name of the sequence in the output file')
    
    args = parser.parse_args()
    
    main(args.input_file, args.output_file, args.old_name, args.new_name)

```



Produce alignment file:

script: minimap2.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=1:00:00
#SBATCH -J minimap2


ref=$1
sample=$2
output=$3

minimap2 -ax asm20 --eqx -t 3 $ref $sample > $output
bam_out=$(echo $output | sed 's@.sam@.bam@g')
samtools view -@ 3 -S -b $output > $bam_out

```

Job:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/01_nucmer_syri_vs_DM/reference/cultivar_seqSyri

chr=10

for chr in $(seq -w 1 12); do echo $chr

ls -1 /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/01_nucmer_syri_vs_DM/nucmer_syri_*_vs_DM/refDM_chr$chr""_vs_qry_*/out.delta | grep -v "_O_" | grep -v "_Y_" | sed 's@out.delta@@g' > list_samples_chr$chr.txt

sed 's@/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/01_nucmer_syri_vs_DM/nucmer_syri_@@g' list_samples_chr$chr.txt | sed 's@_vs_DM/refDM_chr@\t@g' | awk '{print $1}' | sort | uniq > list_sampleName_chr$chr.txt

#samplename=A
for samplename in $(cat list_sampleName_chr$chr.txt) ; do echo $samplename

#hapOrder=1
for hapOrder in $(seq 1 4); do echo $hapOrder

line=$(grep "_"$samplename"_" list_samples_chr$chr.txt | head -n $hapOrder | tail -n 1)
#line=/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/01_nucmer_syri_vs_DM/nucmer_syri_A_vs_DM/refDM_chr01_vs_qry_hap1/

sample=$(echo $line | sed 's@/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/01_nucmer_syri_vs_DM/nucmer_syri_@@g' | sed 's@_vs_DM/refDM_@\t@g' | sed 's@_vs_qry_@\t@g' | sed 's@/syri.out@@g' | awk '{print $1}')
chr=$(echo $line | sed 's@/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/01_nucmer_syri_vs_DM/nucmer_syri_@@g' | sed 's@_vs_DM/refDM_@\t@g' | sed 's@_vs_qry_@\t@g' | sed 's@/syri.out@@g' | awk '{print $2}' | sed 's@chr@@g')
hap=$(echo $line | sed 's@/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/01_nucmer_syri_vs_DM/nucmer_syri_@@g' | sed 's@_vs_DM/refDM_@\t@g' | sed 's@_vs_qry_@\t@g' | sed 's@/@@g' | awk '{print $3}')
path_out=$(echo $line)
hap_syri=$(head -100 $path_out"syri.out" | awk '{print $6}' | sort | uniq | grep -v "-")

culFasta=$(ls -1 ../cultivars/chr$chr"/"$samplename"_"$hap"/"*fasta)
sequenceUsed=$(grep  ">" $culFasta | grep -v PGA | sed 's@>@@g')

python extractSeq.py $culFasta chr$chr"_"$sample"_"$hap.fasta $sequenceUsed $hap_syri

sbatch minimap2.sh ../DM/DM_chr$chr.fa chr$chr"_"$sample"_"$hap.fasta $path_out"out.sam"

done
done
done
# 27.06.2024: 2609119..2609606

```

Run Syri

script: syri.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=1:00:00
#SBATCH -J syri


ref=$1
cul=$2
aligment=$3
sample=$4

mkdir -p $sample

syri --nc 3 -F B --cigar --dir $sample -c $aligment -r $ref -q $cul -k --lf $sample.syri.log


```



job:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/01_nucmer_syri_vs_DM/reference/cultivar_seqSyri

chr=10

mamba activate mymsyd

for chr in $(seq -w 1 12); do echo $chr
#samplename=A
ls -1 /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/01_nucmer_syri_vs_DM/nucmer_syri_*_vs_DM/refDM_chr$chr""_vs_qry_*/out.delta | grep -v "_O_" | grep -v "_Y_" | sed 's@out.delta@@g' > list_samples_chr$chr.txt

sed 's@/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/01_nucmer_syri_vs_DM/nucmer_syri_@@g' list_samples_chr$chr.txt | sed 's@_vs_DM/refDM_chr@\t@g' | awk '{print $1}' | sort | uniq > list_sampleName_chr$chr.txt

for samplename in $(cat list_sampleName_chr$chr.txt) ; do echo $samplename
#hapOrder=1
for hapOrder in $(seq 1 4); do echo $hapOrder

line=$(grep "_"$samplename"_" list_samples_chr$chr.txt | head -n $hapOrder | tail -n 1)
#line=/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/01_nucmer_syri_vs_DM/nucmer_syri_A_vs_DM/refDM_chr01_vs_qry_hap1/
sample=$(echo $line | sed 's@/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/01_nucmer_syri_vs_DM/nucmer_syri_@@g' | sed 's@_vs_DM/refDM_@\t@g' | sed 's@_vs_qry_@\t@g' | sed 's@/syri.out@@g' | awk '{print $1}')
chr=$(echo $line | sed 's@/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/01_nucmer_syri_vs_DM/nucmer_syri_@@g' | sed 's@_vs_DM/refDM_@\t@g' | sed 's@_vs_qry_@\t@g' | sed 's@/syri.out@@g' | awk '{print $2}' | sed 's@chr@@g')
hap=$(echo $line | sed 's@/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/01_nucmer_syri_vs_DM/nucmer_syri_@@g' | sed 's@_vs_DM/refDM_@\t@g' | sed 's@_vs_qry_@\t@g' | sed 's@/@@g' | awk '{print $3}')
path_out=$(echo $line)

sbatch syri.sh ../01_nucmer_syri_vs_DM/reference/DM/DM_chr$chr.fa ../01_nucmer_syri_vs_DM/reference/cultivar_seqSyri/chr$chr"_"$sample"_"$hap.fasta $path_out"out.bam" Syri_chr$chr"_"$sample"_"$hap 

done
done
done
# 27.06.2024: 2609713..2610192

```



Running Msyd

script: msyd.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2:00:00
#SBATCH -J msyd


samples=$1
output=$2

msyd call -i $samples -o $output

```



job:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/10_SVComparison

ls -1 /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/01_nucmer_syri_vs_DM/


mamba activate mymsyd

#chr=01
for chr in $(seq -w 1 12); do echo $chr
echo -e "#name\taln\tsyri\tvcf" > genomes_chr$chr.tsv

ls -1 /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/01_nucmer_syri_vs_DM/nucmer_syri_*_vs_DM/refDM_chr$chr""_vs_qry_*/out.delta | grep -v "_O_" | grep -v "_Y_" | sed 's@out.delta@@g' > list_samples_chr$chr.txt

for line in $(cat list_samples_chr$chr.txt); do echo $line
name=$(echo $line | sed 's@/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/01_nucmer_syri_vs_DM/nucmer_syri_@@g' | sed 's@_vs_DM/refDM_chr.*_vs_qry_@_@g' | sed 's@/@@g')
fasta=$(echo /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/01_nucmer_syri_vs_DM/reference/cultivar_seqSyri/chr$chr"_"$name".fasta")
echo -e $name"\t"$line"out.bam\tSyri_chr"$chr"_"$name"/syri.out\tSyri_chr"$chr"_"$name"/syri.vcf\t"$fasta >> genomes_chr$chr.tsv
done

sbatch msyd.sh genomes_chr$chr.tsv msyd_genomes_chr$chr.pff
done
# 27306.2024: 2610193..2610204

```

Variation of syntonic regions along the genome:

transfor file:

script: process_msyd.py

```python
import sys

def process_sv_table(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        header = infile.readline()
        outfile.write(header)

        for line in infile:
            columns = line.strip().split('\t')
            new_columns = []

            # Process the columns
            for i, col in enumerate(columns):
                if i < 4:  # Assuming the first four columns are chromosome, start, end, and annotation
                    new_columns.append(col)
                else:
                    if col == "-":
                        new_columns.append('0')
                    else:
                        try:
                            # If it can be converted to a number, keep it as is
                            float(col)
                            new_columns.append(col)
                        except ValueError:
                            # Otherwise, replace with 1
                            new_columns.append('1')
            outfile.write('\t'.join(new_columns) + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python process_sv.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    process_sv_table(input_file, output_file)

```

job:

```bash

python process_msyd.py msyd_genomes_chr10.pff msyd_genomes_chr10_sym.txt

for chr in $(seq -w 1 12); do echo $chr
python process_msyd.py msyd_genomes_chr$chr.pff msyd_genomes_chr$chr"_"sym.txt
done

```

Produce table per window:

script: process_sv_windows.py

```python
import sys

def process_sv_table(input_file, output_file):
    window_size = 1000000  # 1 Mb

    def calculate_pi(sample_data):
        # Calculate genetic diversity (pi)
        n = len(sample_data)
        if n < 2:
            return 0.0
        
        pairwise_differences = 0
        for i in range(n):
            for j in range(i + 1, n):
                if sample_data[i] != sample_data[j]:
                    pairwise_differences += 1
        
        num_comparisons = n * (n - 1) / 2
        pi = pairwise_differences / num_comparisons if num_comparisons > 0 else 0.0
        return pi

    with open(input_file, 'r') as infile:
        header = infile.readline().strip().split('\t')
        samples = header[4:]  # Get sample names
        rows = [line.strip().split('\t') for line in infile]

    # Convert columns to appropriate types
    processed_rows = []
    for row in rows:
        if len(row) < 4:
            print(f"Skipping line due to insufficient columns: {row}")
            continue
        try:
            chrom = row[0]
            start = int(row[1])
            end = int(row[2])
            ann = row[3]
            sample_data = row[4:]
            processed_rows.append([chrom, start, end, ann] + sample_data)
        except ValueError as e:
            print(f"Skipping line due to conversion error: {row}, error: {e}")
            continue

    # Calculate windows and overlap
    results = []
    for row in processed_rows:
        chrom, start, end, ann, *sample_data = row
        window_start = (start // window_size) * window_size
        window_end = window_start + window_size

        while window_start <= end:
            overlap_start = max(start, window_start)
            overlap_end = min(end, window_end)
            overlap_bp = max(0, overlap_end - overlap_start + 1)

            if overlap_bp > 0:
                num_samples_with_1 = sum(1 for x in sample_data if x == '1')
                num_samples_with_0 = sum(1 for x in sample_data if x == '0')
                pi = calculate_pi(sample_data)
                results.append([chrom, window_start, window_end, start, end, overlap_bp, num_samples_with_1, num_samples_with_0, pi])

            window_start += window_size
            window_end += window_size

    # Write results to output file
    with open(output_file, 'w') as outfile:
        outfile.write('CHR\tWINDOW_START\tWINDOW_END\tSTART\tEND\tOVERLAP_BP\tNUM_SAMPLES_WITH_1\tNUM_SAMPLES_WITH_0\tPI\n')
        for result in results:
            outfile.write('\t'.join(map(str, result)) + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python process_sv_windows.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    process_sv_table(input_file, output_file)

```

job:

```bash
python process_sv_windows.py msyd_genomes_chr10_sym.txt win_msyd_genomes_chr10_sym.txt

for chr in $(seq -w 1 12); do echo $chr
python process_sv_windows.py msyd_genomes_chr$chr"_"sym.txt win_msyd_genomes_chr$chr"_"sym.txt
done

```

diversity per window:

script: calculate_weighted_pi.py

```python
import sys

def calculate_weighted_pi(input_file, output_file, window_size=1000000):
    windows = {}
    chromosomes = set()
    
    with open(input_file, 'r') as infile:
        header = infile.readline().strip().split('\t')
        for line in infile:
            columns = line.strip().split('\t')
            chrom = columns[0]
            window_start = int(columns[1])
            overlap_bp = int(columns[5])
            pi = float(columns[8])
            
            window_key = (chrom, window_start)
            if window_key not in windows:
                windows[window_key] = {'total_weighted_pi': 0, 'total_bp': 0}
            
            weight = overlap_bp / window_size
            weighted_pi = pi * weight
            
            windows[window_key]['total_weighted_pi'] += weighted_pi
            windows[window_key]['total_bp'] += overlap_bp
            chromosomes.add(chrom)

    # Ensure all windows are present
    for chrom in chromosomes:
        max_end = max(key[1] for key in windows if key[0] == chrom)
        for start in range(0, max_end, window_size):
            window_key = (chrom, start)
            if window_key not in windows:
                windows[window_key] = {'total_weighted_pi': 0, 'total_bp': 0}

    with open(output_file, 'w') as outfile:
        outfile.write('chr\tpos\tpi\n')
        for window_key in sorted(windows.keys()):
            chrom, window_start = window_key
            total_weighted_pi = windows[window_key]['total_weighted_pi']
            total_bp = windows[window_key]['total_bp']
            average_pi = total_weighted_pi if total_bp > 0 else 0
            outfile.write(f'{chrom}\t{window_start}\t{average_pi}\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python calculate_weighted_pi.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    calculate_weighted_pi(input_file, output_file)

```

job:

```bash

python calculate_weighted_pi.py win_msyd_genomes_chr01_sym.txt mean_win_msyd_genomes_chr01_sym.txt

for chr in $(seq -w 1 12); do echo $chr
python calculate_weighted_pi.py win_msyd_genomes_chr$chr"_"sym.txt mean_win_msyd_genomes_chr$chr"_"sym.txt
done

```



# Correlation between genetic diversity (Pi) and number of haplotypes/LD/Recombination Rate



```R
## Model including diversity in haplotype ancestry:

setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes")

library(ggplot2)
library(tidyverse)
library(car)
library(lme4)



# loading Pi table
table_pop_pi <- read_table("AllSamplesMultiAlleles_Win10000_ConsideringmissingData_AllCHR_popParameters.txt", F)
names(table_pop_pi) <- c("chr", "win", "sum_pi", "SS", "NVar", "mean_pi_pw", "W_theta", "Taj_D", "mean_Nsamples")

head(table_pop_pi)

winSize = 500000
table_pop_pi_winSize <- table_pop_pi %>%
  mutate(pos=floor(win/winSize)*winSize) %>%
  group_by(chr, pos) %>%
  dplyr::summarise(meanPi=mean(mean_pi_pw)) 
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
  left_join(meanLD_table, by=c("chr", "pos")) %>% 
  filter(MeanNHaplotypes<20)

head(table_pop_merged)


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


table_pop_merged %>%
  # filter(NAN_haplo == 0) %>%
  ggplot(aes(MeanNHaplotypes, meanPi))+
  geom_point(alpha=0.4)+
  geom_smooth(method='lm', formula= y~x, alpha=0.5, colour="black")


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
  mutate(meanPi=0.035, 
         MeanNHaplotypes=10)

stats_pi_NHap %>% 
  filter(chr!="All")


table_pop_merged %>% 
  # filter(MeanNHaplotypes<20) %>%
  ggplot(aes(MeanNHaplotypes, meanPi))+
  geom_point(alpha=0.4)+
  geom_smooth(method='lm', formula= y~x, alpha=0.5, colour="black")+
  geom_text(data=stats_pi_NHap %>% filter(chr!="All"), 
            aes(label=label), size=3)+
  # facet_grid(chr ~ .)+
  facet_grid(.~chr)+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor.x = element_line(color = "grey80"),
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12)) 

ggsave("21_CorNHap_Pi_mean500kb_perChr.png", height=4, width=18)
ggsave("21_CorNHap_Pi_mean500kb_perChr.pdf", height=4, width=18)


stats_pi_NHap <- read.csv("summary_lm_Pi_MeanNHaplotypes.csv", T) %>%
  mutate(label = paste0("AdjR2=", format(round(adj_r2, 2), nsmall = 2), 
                        "\nP-Val=", format(signif(p.value, digits=3), scientific = T, nsmall = 2))) %>%
  mutate(meanPi=0.035, 
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

ggsave("22_CorNHap_Pi_mean500kb.png", height=4, width=5)
ggsave("22_CorNHap_Pi_mean500kb.pdf", height=4, width=5)







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
  mutate(meanPi=0.035, 
         MeanMaxNShared=15)

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

ggsave("23_CorMeanMaxNShared_Pi_mean500kb_perChr.png", height=4, width=18)
ggsave("23_CorMeanMaxNShared_Pi_mean500kb_perChr.pdf", height=4, width=18)


stats_pi_MeanMaxNShared <- read.csv("./summary_lm_Pi_MeanMaxNShared.csv", T) %>%
  mutate(label = paste0("AdjR2=", format(round(adj_r2, 2), nsmall = 2), 
                        "\nP-Val=", format(signif(p.value, digits=3), scientific = T, nsmall = 2))) %>%
  mutate(meanPi=0.03, 
         MeanMaxNShared=20)

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

ggsave("24_CorMeanMaxNShared_Pi_mean500kb.png", height=4, width=5)
ggsave("24_CorMeanMaxNShared_Pi_mean500kb.pdf", height=4, width=5)


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
  mutate(meanPi=0.035, 
         mean_r2=0.5)



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

ggsave("24_CorR2_Pi_mean500kb_perChr.png", height=4, width=18)
ggsave("24_CorR2_Pi_mean500kb_perChr.pdf", height=4, width=18)



stats_pi_NHap <- read.csv("summary_lm_Pi_R2.csv", T) %>%
  mutate(label = paste0("AdjR2=", format(round(adj_r2, 2), nsmall = 2), 
                        "\nP-Val=", format(signif(p.value, digits=3), scientific = T, nsmall = 2))) %>%
  mutate(meanPi=0.03, 
         mean_r2=0.5)

table_pop_merged %>%
  ggplot(aes(mean_r2, meanPi))+
  geom_point(alpha=0.4)+
  geom_smooth(method='lm', formula= y~x, alpha=0.5, colour="black")+
  geom_text(data=stats_pi_NHap %>% filter(chr=="All"), 
            aes(label=label), size=3)+
  # ylim(0, NA)+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor.x = element_line(color = "grey80"),
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12)) 

ggsave("25_CorR2_Pi_mean500kb.png", height=4, width=5)
ggsave("25_CorR2_Pi_mean500kb.pdf", height=4, width=5)












# loading recombination values

setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/08_recombinationRateCorrelations/plots")
Rec_table <- read_table(paste0("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/08_recombinationRateCorrelations/rec_per1MbWin.txt"), T)

winSize=1000000
meanRec_table <- Rec_table %>%
  mutate(chr=ifelse(Chrom<10, paste0("chr0", Chrom), paste0("chr", Chrom))) %>%
  mutate(pos=(Window-1)*winSize) %>%
  mutate(RecombinationRate=GenDistance) %>%
  select(chr, pos, RecombinationRate)

head(meanRec_table)

meanRec_table %>%
  # filter(RecombinationRate<100) %>%
  ggplot(aes(pos/1000000, (RecombinationRate)))+
  geom_point(colour="#b01e29") +
  # geom_area(aes(fill=factor(low_pi))) +
  geom_line(alpha=0.8)+
  # scale_y_continuous(breaks=seq(0,1,0.01), limits = c(0,0.035)) +
  scale_x_continuous(breaks=seq(0,100,15)) +
  facet_grid(. ~ chr, scale="free", space="free")+
  labs(x = "Pos. (Mb)",
       y = "Recombination (cM/Mb)") +
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12), 
        axis.title = element_text(size=14))

ggsave("01_ed_recombination_mean1Mb.png", height=2.5, width=15)
ggsave("01_ed_recombination_mean1Mb.pdf", height=2.5, width=15)





head(table_pop_merged)

table_pop_merged2 <- table_pop_merged %>%
  mutate(pos=floor(pos/winSize)*winSize) %>%
  group_by(chr, pos) %>%
  dplyr::summarise(meanPi=mean(meanPi), 
                   mean_r2=mean(mean_r2)) %>%
  left_join(meanRec_table, by=c("chr", "pos")) 



## Recombination Vs Pi

head(table_pop_merged2)
m4_Pi_Rec <-lm(meanPi ~ RecombinationRate, data=table_pop_merged2)
Anova(m4_Pi_Rec, Type="III") ; summary(m4_Pi_Rec)

chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

data.frame(chr="chr", 
           r2="r2", 
           adj_r2="adj_r2",
           p.value="p.value") %>%
  write.table(file=paste0("./summary_lm_Pi_Rec.csv"), sep = ",", row.names = F, col.names = F, append=F)

for (chrID in chr_list) {
  chrUsed=paste0("chr",chrID)
  m4_Pi_Rec<-lm(meanPi ~ RecombinationRate, data=table_pop_merged2 %>% 
                  filter(chr==chrUsed) ) 
  #Anova(m4_Pi_LD, Type="III") ; 
  res <- summary(m4_Pi_Rec)
  p.value<-res$coefficients[2,4]
  r2 <- res$r.squared
  adj_r2 <- res$adj.r.squared
  data.frame(chr=chrUsed, 
             r2=r2, 
             adj_r2=adj_r2,
             p.value=p.value) %>%
    write.table(file=paste0("./summary_lm_Pi_Rec.csv"), sep = ",", row.names = F, col.names = F, append=T)
}


m4_Pi_Rec<-lm(meanPi ~ RecombinationRate, data=table_pop_merged2)
res <- summary(m4_Pi_Rec)
p.value<-res$coefficients[2,4]
r2 <- res$r.squared
adj_r2 <- res$adj.r.squared
data.frame(chr="All", 
           r2=r2, 
           adj_r2=adj_r2,
           p.value=p.value) %>%
  write.table(file=paste0("./summary_lm_Pi_Rec.csv"), sep = ",", row.names = F, col.names = F, append=T)



stats_pi_Rec <- read.csv("summary_lm_Pi_Rec.csv", T) %>%
  mutate(label = paste0("AdjR2=", format(round(adj_r2, 2), nsmall = 2), 
                        "\nP-Val=", format(signif(p.value, digits=3), scientific = T, nsmall = 2))) %>%
  mutate(meanPi=0.035, 
         RecombinationRate=5)



table_pop_merged2 %>% 
  ggplot(aes(RecombinationRate, meanPi))+
  geom_point(alpha=0.4)+
  geom_smooth(method='lm', formula= y~x, alpha=0.5, colour="black")+
  geom_text(data=stats_pi_Rec %>% filter(chr!="All"),
            aes(label=label), size=3)+
  ylim(0, NA)+
  facet_grid(. ~ chr)+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor.x = element_line(color = "grey80"),
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12)) 

ggsave("02_CorRec_Pi_mean1Mb_perChr.png", height=4, width=18)
ggsave("02_CorRec_Pi_mean1Mb_perChr.pdf", height=4, width=18)



stats_pi_NHap <- read.csv("summary_lm_Pi_Rec.csv", T) %>%
  mutate(label = paste0("AdjR2=", format(round(adj_r2, 2), nsmall = 2), 
                        "\nP-Val=", format(signif(p.value, digits=3), scientific = T, nsmall = 2))) %>%
  mutate(meanPi=0.035, 
         RecombinationRate=5)

table_pop_merged2 %>%
  ggplot(aes(RecombinationRate, meanPi))+
  geom_point(alpha=0.4)+
  geom_smooth(method='lm', formula= y~x, alpha=0.5, colour="black")+
  geom_text(data=stats_pi_NHap %>% filter(chr=="All"), 
            aes(label=label), size=3)+
  # ylim(0, NA)+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor.x = element_line(color = "grey80"),
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12)) 

ggsave("03_CorRec_Pi_mean1Mb.png", height=4, width=5)
ggsave("03_CorRec_Pi_mean1Mb.pdf", height=4, width=5)









## Recombination Vs LD

head(combined_ld_rec_data2)
m4_LD_Rec <-lm(r2 ~ CumulativeRecombinationRate, data=combined_ld_rec_data2)
Anova(m4_LD_Rec, Type="III") ; summary(m4_LD_Rec)

chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

data.frame(chr="chr", 
           r2="r2", 
           adj_r2="adj_r2",
           p.value="p.value") %>%
  write.table(file=paste0("./summary_lm_LD_Rec.csv"), sep = ",", row.names = F, col.names = F, append=F)

for (chrID in chr_list) {
  chrUsed=paste0("chr",chrID)
  m4_LD_Rec<-lm(r2 ~ CumulativeRecombinationRate, data=combined_ld_rec_data2 %>% 
                  filter(chr==chrUsed) ) 
  #Anova(m4_Pi_LD, Type="III") ; 
  res <- summary(m4_LD_Rec)
  p.value<-res$coefficients[2,4]
  r2 <- res$r.squared
  adj_r2 <- res$adj.r.squared
  data.frame(chr=chrUsed, 
             r2=r2, 
             adj_r2=adj_r2,
             p.value=p.value) %>%
    write.table(file=paste0("./summary_lm_LD_Rec.csv"), sep = ",", row.names = F, col.names = F, append=T)
}


m4_LD_Rec<-lm(r2 ~ CumulativeRecombinationRate, data=combined_ld_rec_data2)
res <- summary(m4_LD_Rec)
p.value<-res$coefficients[2,4]
r2 <- res$r.squared
adj_r2 <- res$adj.r.squared
data.frame(chr="All", 
           r2=r2, 
           adj_r2=adj_r2,
           p.value=p.value) %>%
  write.table(file=paste0("./summary_lm_LD_Rec.csv"), sep = ",", row.names = F, col.names = F, append=T)



stats_LD_Rec <- read.csv("summary_lm_LD_Rec.csv", T) %>%
  mutate(label = paste0("AdjR2=", format(round(adj_r2, 2), nsmall = 2), 
                        "\nP-Val=", format(signif(p.value, digits=3), scientific = T, nsmall = 2))) %>%
  mutate(r2=0.75, 
         CumulativeRecombinationRate=10)



combined_ld_rec_data2 %>%
  filter(distance<2000000) %>%
  ggplot(aes(CumulativeRecombinationRate, r2))+
  geom_point(alpha=0.2)+
  geom_smooth(method='lm', formula= y~x, alpha=0.5, colour="black")+
  geom_text(data=stats_LD_Rec %>% filter(chr!="All"),
            aes(label=label), size=3)+
  ylim(0, NA)+
  facet_grid(. ~ chr)+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor.x = element_line(color = "grey80"),
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12)) 

ggsave("02_CorRec_LD_mean1Mb_perChr.png", height=4, width=18)
ggsave("02_CorRec_LD_mean1Mb_perChr.pdf", height=4, width=18)
getwd()


stats_LD_Rec <- read.csv("summary_lm_LD_Rec.csv", T) %>%
  mutate(label = paste0("AdjR2=", format(round(adj_r2, 2), nsmall = 2), 
                        "\nP-Val=", format(signif(p.value, digits=3), scientific = T, nsmall = 2))) %>%
  mutate(r2=0.75, 
         CumulativeRecombinationRate=10)

combined_ld_rec_data2 %>%
  filter(distance<2000000) %>%
  ggplot(aes(CumulativeRecombinationRate, r2))+
  geom_point(alpha=0.2)+
  # geom_smooth(method='loess', formula= y~x, alpha=0.5, colour="black")+
  geom_smooth(method='gam', formula= y ~ s(x, bs = "cs"), alpha=0.5, colour="black")+
  geom_text(data=stats_LD_Rec %>% filter(chr=="All"), 
            aes(label=label), size=3)+
  # ylim(0, NA)+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor.x = element_line(color = "grey80"),
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12)) 

ggsave("03_CorRec_LD_mean1Mb.png", height=4, width=5)
ggsave("03_CorRec_LD_mean1Mb.pdf", height=4, width=5)













# correlation LD Vs recombination 
# residual analyses:

### 
setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/08_recombinationRateCorrelations/plots")
library(dplyr)
library(zoo)

recomb_data <- table_pop_merged2 %>%
  select(chr, pos, RecombinationRate)


head(recomb_data)


# Function to calculate cumulative recombination rate between two positions
get_cumulative_recomb_rate <- function(chrom, pos1, pos2, recomb_data) {
  sub_data <- recomb_data %>%
    filter(chr == chrom & pos >= pos1 & pos <= pos2)
  cum_recomb_rate <- sum(sub_data$RecombinationRate, na.rm = T)
  return(cum_recomb_rate)
}


chrID="12"
# Define the list of chromosomes
chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

# Loop through each chromosome and process the data
combined_ld_rec_data <- data.frame()  ; for (chrID in chr_list) {
  # Load LD data
  ld_data <- read.table(paste0("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/04_LD/MeanLDValue_chr", chrID, "_maxDis100Mb_WinSize1000000.txt"), 
                        header=FALSE, stringsAsFactors=FALSE) 
  
  names(ld_data) <- c("pos1", "pos2", "r2", "NVar")
  
  # Apply to LD data
  merged_ld_data <- ld_data %>%
    mutate(chr=paste0("chr", chrID)) %>%
    rowwise() %>%
    mutate(CumulativeRecombinationRate = get_cumulative_recomb_rate(chr, pos1, pos2, recomb_data)) %>%
    ungroup()
  
  combined_ld_rec_data <- bind_rows(combined_ld_rec_data, merged_ld_data)
  
}


combined_ld_rec_data2 <- combined_ld_rec_data %>%
  # filter(pos1!=pos2) %>%
  mutate(distance=ifelse(pos1<pos2, pos2-pos1, pos1-pos2)) 



# Plot LD Decay Curves

combined_ld_rec_data2 %>%
  ggplot(aes(x=distance/1000000, y=r2, color=CumulativeRecombinationRate)) +
  geom_point(alpha=0.05) +
  geom_smooth(method="loess", colour="black") +
  scale_color_gradient(low="blue", high="red") +
  labs(title="LD Decay with Cumulative\nRecombination Rate", 
       x="Distance (bp)", y="LD (r2)", 
       colour='Cumulative Rec. R.')+
  theme_classic()+
  theme(axis.text.y = element_text(size=18, colour="black"),
        axis.text.x = element_text(size=18, colour="black"),
        axis.title = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=20),
        legend.position = "right")

ggsave(filename=paste0("./04_Distance_Vs_LD_Recombination.png"),
       device = png, units="cm",
       width=25, height=15)
ggsave(filename=paste0("./04_Distance_Vs_LD_Recombination.pdf"),
       device = pdf, units="cm",
       width=25, height=15)




combined_ld_rec_data2 %>%
  ggplot(aes(x=distance/1000000, y=r2, color=CumulativeRecombinationRate)) +
  geom_point(alpha=0.05) +
  geom_smooth(method="loess", colour="black") +
  scale_color_gradient(low="blue", high="red") +
  labs(title="LD Decay with Cumulative\nRecombination Rate", 
       x="Distance (bp)", y="LD (r2)", 
       colour='Cumulative Rec. R.')+
  facet_grid(chr ~ .)+
  theme_classic()+
  theme(axis.text.y = element_text(size=18, colour="black"),
        axis.text.x = element_text(size=18, colour="black"),
        axis.title = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=20),
        legend.position = "right")

ggsave(filename=paste0("./05_Distance_Vs_LD_Recombination_perChr.png"),
       device = png, units="cm",
       width=20, height=33)
ggsave(filename=paste0("./05_Distance_Vs_LD_Recombination_perChr.pdf"),
       device = pdf, units="cm",
       width=20, height=33)




# Simple linear regression including cumulative recombination rate
head(combined_ld_rec_data2)
model <- lm(r2 ~ distance + CumulativeRecombinationRate, data=combined_ld_rec_data2)
summary(model)


combined_ld_rec_data2$residuals <- residuals(model)


ggplot(combined_ld_rec_data2, aes(x=distance, y=residuals, color=CumulativeRecombinationRate)) +
  geom_point() +
  geom_smooth(method="loess") +
  scale_color_gradient(low="blue", high="red") +
  labs(title="Residuals of LD Adjusted for Recombination Rate", x="Distance (bp)", y="Residuals")+
  facet_grid(chr ~ .)+
  theme_classic()+
  theme(axis.text.y = element_text(size=18, colour="black"),
        axis.text.x = element_text(size=18, colour="black"), 
        panel.grid.major.y = element_line(),
        axis.title = element_text(size=20), 
        legend.text = element_text(size=18),
        legend.title = element_text(size=20),
        legend.position = "right")


combined_ld_rec_data2 %>%
  # filter(pos1 != pos2) %>%
  # filter(chr=="chr10") %>%
  # filter(distance<4000000) %>%
  mutate(residuals2=ifelse(residuals<0, 0, residuals)) %>%
  ggplot(aes(pos1/1000000, pos2/1000000, fill=residuals2, colour=residuals2)) +
  geom_tile()+
  #scico::scale_fill_scico(palette = "lajolla")+
  scale_fill_viridis_c(option = "magma", direction = -1)+
  scale_colour_viridis_c(option = "magma", direction = -1, guide = 'none')+
  # scale_fill_viridis_c(option = "magma", direction = -1, limits = c(-0.25, 1))+
  # scale_colour_viridis_c(option = "magma", direction = -1, limits = c(-0.25, 1), guide = 'none')+
  scale_x_continuous(breaks=seq(0,100,15), name = "Position (Mb)", expand = c(0,0))+
  scale_y_continuous(breaks=seq(0,100,15), name = "Position (Mb)", expand = c(0,0))+
  # scale_y_discrete(name = "Haplotype", expand = c(0,0))+
  # ggtitle(paste0("Chr:", chrID))+
  # labs(fill="LD (R^2)")+
  facet_grid(. ~ chr, scale="free", space="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"), 
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12), 
        axis.title = element_text(size=14))


ggsave("06_Residual_MeanLD_heatmap.png", height=3, width=20)
ggsave("06_Residual_MeanLD_heatmap.pdf", height=3, width=20)



mean_combined_ld_rec_data2 <- combined_ld_rec_data2 %>%
  filter(distance<4000000) %>%
  mutate(pos=pos1) %>%
  group_by(chr, pos) %>%
  dplyr::summarise(mean_residuals=mean(residuals)) 


quantile_95 <- mean_combined_ld_rec_data2 %>% 
  pull(mean_residuals) %>%
  quantile(0.95)

quantile_95


centromere_table <- read.csv("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/04_LD/centromeres.csv", header = T)
head(centromere_table)


mean_combined_ld_rec_data2 %>%
  ggplot(aes(pos/1000000, mean_residuals))+
  geom_ribbon(aes(ymax=ifelse(mean_residuals>quantile_95, mean_residuals, quantile_95),
                  ymin=quantile_95), fill="#91b873") +
  # geom_ribbon(aes(ymax=ifelse(mean_residuals>0, mean_residuals, 0),
  #                 ymin=0), fill="#91b873") +
  geom_line(alpha=0.8)+
  geom_point(data=mean_combined_ld_rec_data2 %>%
               filter(mean_residuals>quantile_95) %>%
               mutate(posy=(-0.2)), aes(pos/1000000, posy), colour="#203f08", alpha=0.4) +
  geom_line(data=centromere_table %>%
              gather(key="posB", value = "pos", -chr) %>%
              mutate(posy = -0.3),
            aes(pos/1000000,posy), colour="darkgreen", linewidth=1.5, alpha=0.5) +
  geom_hline(yintercept=quantile_95, linetype="dashed", color = "red")+
  # geom_hline(yintercept=0, linetype="dashed", color = "red")+
  scale_y_continuous(breaks=seq(0,0.5,0.25), labels = c("0.0", "", "0.5")) +
  scale_x_continuous(breaks=seq(0,100,15)) +
  facet_grid(. ~ chr, scale="free", space="free")+
  labs(x = "Pos. (Mb)",
       y = "Residual LD") +
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12), 
        axis.title = element_text(size=14))


ggsave("07_Residual_MeanLD_within4Mb.png", height=2, width=15)
ggsave("07_Residual_MeanLD_within4Mb.pdf", height=2, width=15)





# Model including diversity in haplotype ancestry and SVs:
setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/08_recombinationRateCorrelations/plots")
diversity_hap_ancestry <- read.csv("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/05_ABBA_BABA/DSuite/diversity_hap_ancestry_1MbWin.csv", T) 
head(diversity_hap_ancestry)


combined_ld_rec_data3 <- combined_ld_rec_data %>%
  filter(pos1==pos2) %>%
  mutate(pos=pos1) %>%
  select(chr, pos, r2, CumulativeRecombinationRate) %>%
  left_join(diversity_hap_ancestry, by=c("chr", "pos"))

head(combined_ld_rec_data3)


# loading diversity in SV: 

combined_df_sv <- data.frame()

# Loop through chromosome numbers 1 to 12
for (i in 1:12) {
  # Construct the file name
  file_name <- sprintf("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/10_SVComparison/mean_win_msyd_genomes_chr%02d_sym.txt", i)
  
  # Read the file into a data frame
  chr_df <- read.table(file_name, header = TRUE, sep = "\t")
  
  # Combine the data frame with the previously read data
  combined_df_sv <- bind_rows(combined_df_sv, chr_df)
}

# Print the combined data frame
head(combined_df_sv)
names(combined_df_sv) <- c("chr", "pos", "divSV")


combined_ld_rec_data4 <- combined_ld_rec_data3 %>%
  left_join(combined_df_sv, by=c("chr", "pos"))





combined_ld_rec_data4 %>%
  head()

model2 <- lm(r2 ~ CumulativeRecombinationRate + DivHapAnc , data=combined_ld_rec_data4)
summary(model2)

model3 <- lm(r2 ~ CumulativeRecombinationRate + DivHapAnc + chr, data=combined_ld_rec_data4)
summary(model3)

model4 <- lm(r2 ~ DivHapAnc, data=combined_ld_rec_data4)
summary(model4)



combined_ld_rec_data4$residuals <- residuals(model2)

quantile_95 <- combined_ld_rec_data4 %>% 
  pull(residuals) %>%
  quantile(0.95)

quantile_95


combined_ld_rec_data4 %>%
  ggplot(aes(pos/1000000, residuals))+
  geom_ribbon(aes(ymax=ifelse(residuals>quantile_95, residuals, quantile_95),
                  ymin=quantile_95), fill="#91b873") +
  geom_line(alpha=0.8)+
  geom_point(data=combined_ld_rec_data4 %>%
               filter(residuals>quantile_95) %>%
               mutate(posy=(-0.2)), aes(pos/1000000, posy), colour="#203f08", alpha=0.4) +
  geom_line(data=centromere_table %>%
              ungroup() %>%
              mutate(numeros_str = sprintf("%02d", as.numeric(as.character(chr))),
                     chr = factor(paste0("chr", numeros_str))) %>%
              select(-numeros_str) %>%
              gather(key="posB", value = "pos", -chr) %>%
              mutate(posy = -0.3),
            aes(pos/1000000,posy), colour="darkgreen", linewidth=1.5, alpha=0.5) +
  geom_hline(yintercept=quantile_95, linetype="dashed", color = "red")+
  # geom_hline(yintercept=0, linetype="dashed", color = "red")+
  scale_y_continuous(breaks=seq(0,0.5,0.25), labels = c("0.0", "", "0.5")) +
  scale_x_continuous(breaks=seq(0,100,15)) +
  facet_grid(. ~ chr, scale="free", space="free")+
  labs(x = "Pos. (Mb)",
       y = "Residual LD") +
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12), 
        axis.title = element_text(size=14))


ggsave("08_Residualr2_recombination_divHapAnc_within4Mb.png", height=2, width=15)
ggsave("08_Residualr2_recombination_divHapAnc_within4Mb.pdf", height=2, width=15)



combined_ld_rec_data4 %>%
  ggplot(aes(pos/1000000, divSV))+
  geom_line(alpha=0.8)+
  geom_point(alpha=0.4) +
  # geom_line(data=centromere_table %>%
  #             gather(key="posB", value = "pos", -chr) %>%
  #             mutate(posy = -0.3),
  #           aes(pos/1000000,posy), colour="darkgreen", linewidth=1.5, alpha=0.5) +
  # geom_hline(yintercept=quantile_95, linetype="dashed", color = "red")+
  # scale_y_continuous(breaks=seq(0,0.5,0.25), labels = c("0.0", "", "0.5")) +
  scale_x_continuous(breaks=seq(0,100,15)) +
  facet_grid(. ~ chr, scale="free", space="free")+
  labs(x = "Pos. (Mb)",
       y = "Diversity SV") +
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12), 
        axis.title = element_text(size=14))


ggsave("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/10_SVComparison/plots/01_DivSV_alongGenome_1MbWin.png", height=2, width=15)
ggsave("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/10_SVComparison/plots/01_DivSV_alongGenome_1MbWin.pdf", height=2, width=15)





## SV Vs LD

head(combined_ld_rec_data4)
m4_LD_sv <-lm(r2 ~ divSV, data=combined_ld_rec_data4)
Anova(m4_LD_sv, Type="III") ; summary(m4_LD_sv)

chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

data.frame(chr="chr", 
           r2="r2", 
           adj_r2="adj_r2",
           p.value="p.value") %>%
  write.table(file=paste0("./summary_lm_LD_sv.csv"), sep = ",", row.names = F, col.names = F, append=F)

for (chrID in chr_list) {
  chrUsed=paste0("chr",chrID)
  m4_LD_sv<-lm(r2 ~ divSV, data=combined_ld_rec_data4 %>% 
                  filter(chr==chrUsed) ) 
  #Anova(m4_Pi_LD, Type="III") ; 
  res <- summary(m4_LD_sv)
  p.value<-res$coefficients[2,4]
  r2 <- res$r.squared
  adj_r2 <- res$adj.r.squared
  data.frame(chr=chrUsed, 
             r2=r2, 
             adj_r2=adj_r2,
             p.value=p.value) %>%
    write.table(file=paste0("./summary_lm_LD_sv.csv"), sep = ",", row.names = F, col.names = F, append=T)
}


m4_LD_sv<-lm(r2 ~ divSV, data=combined_ld_rec_data4)
res <- summary(m4_LD_sv)
p.value<-res$coefficients[2,4]
r2 <- res$r.squared
adj_r2 <- res$adj.r.squared
data.frame(chr="All", 
           r2=r2, 
           adj_r2=adj_r2,
           p.value=p.value) %>%
  write.table(file=paste0("./summary_lm_LD_sv.csv"), sep = ",", row.names = F, col.names = F, append=T)



stats_LD_sv <- read.csv("summary_lm_LD_sv.csv", T) %>%
  mutate(label = paste0("AdjR2=", format(round(adj_r2, 2), nsmall = 2), 
                        "\nP-Val=", format(signif(p.value, digits=3), scientific = T, nsmall = 2))) %>%
  mutate(divSV=0.1, 
         r2=0.75)

head(stats_LD_sv)

combined_ld_rec_data4 %>% 
  ggplot(aes(divSV, r2))+
  geom_point(alpha=0.4)+
  geom_smooth(method='lm', formula= y~x, alpha=0.5, colour="black")+
  geom_text(data=stats_LD_sv %>% filter(chr!="All"),
            aes(label=label), size=3)+
  ylim(0, NA)+
  facet_grid(. ~ chr)+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor.x = element_line(color = "grey80"),
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12)) 

getwd()

ggsave("02_CorLD_sv_mean1Mb_perChr.png", height=4, width=18)
ggsave("02_CorLD_sv_mean1Mb_perChr.pdf", height=4, width=18)



stats_LD_sv <- read.csv("summary_lm_LD_sv.csv", T) %>%
  mutate(label = paste0("AdjR2=", format(round(adj_r2, 2), nsmall = 2), 
                        "\nP-Val=", format(signif(p.value, digits=3), scientific = T, nsmall = 2))) %>%
  mutate(divSV=0.1, 
         r2=0.8)


combined_ld_rec_data4 %>% 
  ggplot(aes(divSV, r2))+
  geom_point(alpha=0.4)+
  geom_smooth(method='lm', formula= y~x, alpha=0.5, colour="black")+
  geom_text(data=stats_LD_sv %>% filter(chr=="All"), 
            aes(label=label), size=3)+
  # ylim(0, NA)+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor.x = element_line(color = "grey80"),
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12)) 

ggsave("03_CorLD_sv_mean1Mb.png", height=4, width=5)
ggsave("03_CorLD_sv_mean1Mb.pdf", height=4, width=5)







#using mix models:


# install.packages("lme4")
library(lme4)


head(combined_ld_rec_data4)

model2 <- lm(r2 ~ CumulativeRecombinationRate + DivHapAnc , data=combined_ld_rec_data4)
summary(model2)

model2 <- lm(r2 ~ divSV , data=combined_ld_rec_data4)
summary(model2)

model <- lmer(r2 ~ CumulativeRecombinationRate + divSV + DivHapAnc + (1 | chr), data = combined_ld_rec_data4)
summary(model)


model3 <- lm(r2 ~ CumulativeRecombinationRate + DivHapAnc + chr, data=combined_ld_rec_data4)
summary(model3)

model4 <- lm(r2 ~ DivHapAnc, data=combined_ld_rec_data4)
summary(model4)


model <- lmer(r2 ~ CumulativeRecombinationRate + (1 | chr), data = combined_ld_rec_data4)
summary(model)

model <- lmer(r2 ~ CumulativeRecombinationRate + DivHapAnc + (1 | chr), data = combined_ld_rec_data4)
summary(model)

model <- lmer(r2 ~ DivHapAnc + (1 | chr), data = combined_ld_rec_data4)
summary(model)


model <- lmer(r2 ~ CumulativeRecombinationRate + divSV + DivHapAnc + (1 | chr), data = combined_ld_rec_data4)
summary(model)
library(MuMIn)
r.squaredGLMM(model)


combined_ld_rec_data4$residuals <- residuals(model)

quantile_95 <- combined_ld_rec_data4 %>% 
  pull(residuals) %>%
  quantile(0.95)

quantile_95


combined_ld_rec_data4 %>%
  ggplot(aes(pos/1000000, residuals))+
  geom_ribbon(aes(ymax=ifelse(residuals>quantile_95, residuals, quantile_95),
                  ymin=quantile_95), fill="#91b873") +
  geom_line(alpha=0.8)+
  geom_point(data=combined_ld_rec_data4 %>%
               filter(residuals>quantile_95) %>%
               mutate(posy=(-0.2)), aes(pos/1000000, posy), colour="#203f08", alpha=0.4) +
  geom_line(data=centromere_table %>%
              gather(key="posB", value = "pos", -chr) %>%
              mutate(posy = -0.3),
            aes(pos/1000000,posy), colour="darkgreen", linewidth=1.5, alpha=0.5) +
  geom_hline(yintercept=quantile_95, linetype="dashed", color = "red")+
  # geom_hline(yintercept=0, linetype="dashed", color = "red")+
  scale_y_continuous(breaks=seq(0,0.5,0.25), labels = c("0.0", "", "0.5")) +
  scale_x_continuous(breaks=seq(0,100,15)) +
  facet_grid(. ~ chr, scale="free", space="free")+
  labs(x = "Pos. (Mb)",
       y = "Residual LD") +
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12), 
        axis.title = element_text(size=14))


getwd()
ggsave("09_Residualr2_recombination_divHapAnc_SV_within4Mb.png", height=2, width=15)
ggsave("09_Residualr2_recombination_divHapAnc_SV_within4Mb.pdf", height=2, width=15)








winSize=1000000
head(combined_ld_rec_data4)

combined_ld_rec_data5 <- combined_ld_rec_data4 %>%
  left_join(table_pop_merged %>% 
  mutate(pos=floor(pos/winSize)*winSize) %>%
  group_by(chr, pos) %>%
  dplyr::summarise(meanPi=mean(meanPi), 
                   mean_NHaplotypes=mean(MeanNHaplotypes), 
                   mean_MaxNShared=mean(MeanMaxNShared)), by=c("chr", "pos") )

head(combined_ld_rec_data5)

model2 <- lm(r2 ~ divSV , data=combined_ld_rec_data5)
summary(model2)

model2 <- lm(meanPi ~ divSV , data=combined_ld_rec_data5)
summary(model2)

model2 <- lm(meanPi ~ DivHapAnc, data=combined_ld_rec_data5)
summary(model2)



# Fit the full model
full_model <- lmer(r2 ~ CumulativeRecombinationRate + divSV + DivHapAnc + (1 | chr), data = combined_ld_rec_data5)

# Calculate R^2 for the full model
full_r2 <- r.squaredGLMM(full_model)

# Fit nested models
model_no_cumrecomb <- lmer(r2 ~ divSV + DivHapAnc + (1 | chr), data = combined_ld_rec_data5)
model_no_divsv <- lmer(r2 ~ CumulativeRecombinationRate + DivHapAnc + (1 | chr), data = combined_ld_rec_data5)
model_no_divhapanc <- lmer(r2 ~ CumulativeRecombinationRate + divSV + (1 | chr), data = combined_ld_rec_data5)

# Calculate R^2 for nested models
r2_no_cumrecomb <- r.squaredGLMM(model_no_cumrecomb)
r2_no_divsv <- r.squaredGLMM(model_no_divsv)
r2_no_divhapanc <- r.squaredGLMM(model_no_divhapanc)

# Extract marginal R^2 (variance explained by fixed effects)
full_r2_marginal <- full_r2[1]
r2_no_cumrecomb_marginal <- r2_no_cumrecomb[1]
r2_no_divsv_marginal <- r2_no_divsv[1]
r2_no_divhapanc_marginal <- r2_no_divhapanc[1]

# Calculate contributions
contribution_cumrecomb <- full_r2_marginal - r2_no_cumrecomb_marginal
contribution_divsv <- full_r2_marginal - r2_no_divsv_marginal
contribution_divhapanc <- full_r2_marginal - r2_no_divhapanc_marginal

# Print contributions
cat("Contribution of CumulativeRecombinationRate:", contribution_cumrecomb, "\n")
cat("Contribution of divSV:", contribution_divsv, "\n")
cat("Contribution of DivHapAnc:", contribution_divhapanc, "\n")



```





# Pan-genome collapsing shared haplotypes



## Clustering using overlapping Windows and remove missing data:

These scripts are the same as in the section"Clustering using overlapping Windows and remove missing data" but using more relaxed limits of the minimum number of allowed SS to differentiate haplotypes. In this case, for the pan genome the limit was increased to 55 SS. 

Version using overlaping windows:

script: ClusteringHaplotypes_overlappingWin.py

```python
import os
import sys
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as shc
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
import math

# genotype_file = str("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr01_missing_genotypeTable.txt")
# win_size = int(10000)
# output_name = str("test")
# max_d = int(50)

# this script takes a genotype table, a window size and and output name and: produce a distance matrix and cluster haplotypes with a maximum of around 10SNPs (max_d parameter).
# the output is a table with window start position, sample and haplotype ID after clustering.
# The scripts uses the haplotype groups produced from the previous window to allocate haplotypes ID in the current window. 

# genotype_file = str("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr01_missing_genotypeTable.txt")
# win_size = int(1000000)
# step_size= int(50000)
# max_d = int(2000)
# output_name = str("test")


genotype_file = str(sys.argv[1])
win_size = int(sys.argv[2])
step_size=int(sys.argv[3])
max_d = int(sys.argv[4])
output_name = str(sys.argv[5])

def str_to_int(x):
    if x!='nan':
        return int(x)
    else:
        return int(-1)

dic_win_SS_array = {}
dic_win_NA_genotypes_mean = {}
dic_win_NVar = {}
break_var = 0

for line in open(genotype_file, "r"):
    if line[0:3] == "CHR":
        haplotypes = line.strip().split("\t")[5:]
        list_win_pre = [0] 
        dic_win_SS_array[0] = np.zeros((len(haplotypes), len(haplotypes)))
        dic_win_NA_genotypes_mean[0] = np.zeros((len(haplotypes)))
        dic_win_NVar[0] = 0
        table_hap_groups_pre = pd.DataFrame()
        if os.path.exists(f"{output_name}.txt"):
            os.remove(f"{output_name}.txt")
    else:
        line_values = line.strip().split("\t")
        pos = int(line_values[1])
        list_win = [i for i in range(int((pos-win_size)/step_size),int(pos/step_size)+1) if i>=0]
        list_win_withPre = set(list_win_pre + list_win)
        genotypes = np.asarray([str_to_int(x) for x in line_values[5:]])
        SS_var = genotypes[:,None] != genotypes
        NA_genotypes = (genotypes==(-1))*1
        NA_genotypes_boolean = (genotypes==(-1))
        SS_var[NA_genotypes_boolean, :] = False
        SS_var[:, NA_genotypes_boolean] = False
        for win_current in list_win_withPre:
            if not win_current in dic_win_SS_array:
                dic_win_SS_array[win_current] = np.zeros((len(haplotypes), len(haplotypes)))
                dic_win_NA_genotypes_mean[win_current] = np.zeros((len(haplotypes)))
                dic_win_NVar[win_current] = 0
            if win_current in list_win:
                dic_win_SS_array[win_current] += SS_var
                dic_win_NA_genotypes_mean[win_current] += NA_genotypes
                dic_win_NVar[win_current] += 1
            else:
                #Z = shc.ward(pdist(dic_win_SS_array[win_current]))
                condensed_dist = squareform(dic_win_SS_array[win_current])
                linkresult = shc.linkage(condensed_dist, method='complete')
                clusters = shc.fcluster(linkresult, max_d, criterion='distance')
                #Z = shc.ward(dic_win_SS_array[win_current])
                #clusters = shc.fcluster(condensed_dist, max_d, criterion='distance')
                if dic_win_NVar[win_current] > 0:
                    frac_NA_genotypes_mean = dic_win_NA_genotypes_mean[win_current]/dic_win_NVar[win_current]
                else:
                    frac_NA_genotypes_mean = dic_win_NA_genotypes_mean[win_current]
                if not table_hap_groups_pre.empty:
                    max_cluster_pre = table_hap_groups_pre['haplotype_group_pre'].max()
                    table_hap_groups = pd.DataFrame(data={'labs': list(haplotypes), 'haplotype_group': list(clusters+max_cluster_pre), 'Mean_NAN_haplo': list(frac_NA_genotypes_mean)})
                    comparisonDF = table_hap_groups.merge(table_hap_groups_pre, on='labs', how='left').dropna().groupby(['haplotype_group', 'haplotype_group_pre']).agg({'count'})['labs'].reset_index()
                    idx = comparisonDF.groupby('haplotype_group')['count'].transform(max) == comparisonDF['count']
                    idx2 = comparisonDF[idx].groupby('haplotype_group_pre')['count'].transform(max) == comparisonDF[idx]['count']
                    previous_hap_groups = comparisonDF[idx][idx2].groupby('haplotype_group_pre').first().reset_index().groupby('haplotype_group').first().reset_index()
                    table_hap_groups_merge = table_hap_groups.merge( previous_hap_groups[['haplotype_group', 'haplotype_group_pre']], on='haplotype_group', how='left')
                    table_hap_groups_merge['haplotype_group_ed'] = [int(table_hap_groups_merge['haplotype_group'][x]) if math.isnan(table_hap_groups_merge['haplotype_group_pre'][x]) else int(table_hap_groups_merge['haplotype_group_pre'][x]) for x in range(0, table_hap_groups_merge.shape[0])]
                    table_hap_groups_merge['win'] = win_current*step_size
                    table_hap_groups_merge[['win', 'labs', 'haplotype_group_ed', 'Mean_NAN_haplo']].to_csv(f"{output_name}.txt", sep="\t", index=False, header=False, mode='a')
                    table_hap_groups_pre = table_hap_groups_merge[['labs', 'haplotype_group_ed']].rename(columns={"labs": "labs", "haplotype_group_ed": "haplotype_group_pre"})
                    del dic_win_SS_array[win_current]
                    del dic_win_NA_genotypes_mean[win_current]
                    del dic_win_NVar[win_current]
                else:
                    table_hap_groups = pd.DataFrame(data={'labs': list(haplotypes), 'haplotype_group': list(clusters), 'Mean_NAN_haplo': list(frac_NA_genotypes_mean)})
                    table_hap_groups['win'] = win_current*step_size
                    table_hap_groups[['win', 'labs', 'haplotype_group', 'Mean_NAN_haplo']].to_csv(f"{output_name}.txt", sep="\t", index=False, header=False, mode='a')
                    table_hap_groups_pre = table_hap_groups[['labs', 'haplotype_group']].rename(columns={"labs": "labs", "haplotype_group": "haplotype_group_pre"})
                    del dic_win_SS_array[win_current]
                    del dic_win_NA_genotypes_mean[win_current]
                    del dic_win_NVar[win_current]
        list_win_pre = list_win


```



bash script: ClusteringHaplotypes_overlappingWin.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00
#SBATCH -J clusterHaplotypes

genotype_file=$1
win_size=$2
step_size=$3
maxdif=$4
output_name=$5

# python3 /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesPopPar_Win10000_chr01_missing_genotypeTable.txt 1000000 50000 2000 Cluster_hapIDs_chr01
python3 ClusteringHaplotypes_overlappingWin.py $genotype_file $win_size $step_size $maxdif $output_name 

```



Job:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/06_mergedHaplotypes_PanGenome

mamba activate mypython3

#/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr$chN"_"genotypeTable_miltipleAlleles.txt

# using window size of 10 kb and min number of SS 50. Step size 5000.
for chN in $( seq -w 1 12); do echo $chN
sbatch ClusteringHaplotypes_overlappingWin.sh /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr$chN"_"missing_genotypeTable.txt 10000 5000 50 haplotypeGroups_min50SS_10KbWin_chr$chN 
done
# 06.04.2024: 2456317..2456328


# using window size of 10 kb and min number of SS 20. Step size 5000.
for chN in $( seq -w 1 12); do echo $chN
sbatch ClusteringHaplotypes_overlappingWin.sh /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr$chN"_"missing_genotypeTable.txt 10000 5000 20 haplotypeGroups_min20SS_10KbWin_chr$chN 
done
# 06.04.2024: 2456329..2456340


#for i in {2199995..2200006}; do stop $i ; rm slurm-$i* ; done 

```



## Calculate number of difference per window:

difference_perLargeWin.py

```python
import os
import sys
import numpy as np
import pandas as pd

import scipy.cluster.hierarchy as shc
from scipy.spatial.distance import pdist, squareform
import math

input = 'haplotypeGroups_min10SS_chr06.txt'
winSize=int(1000000)
step_size = int(50000)
output="test"
chr="chr06"

input = str(sys.argv[1])
winSize = int(sys.argv[2])
step_size = int(sys.argv[3])
output = str(sys.argv[4])
chr = str(sys.argv[5])

table = pd.read_table(input, sep="\t", names=['pos', 'sample', 'haplotype', 'missingP'])  
tableSpread = table[['pos', 'sample', 'haplotype']].pivot(index = 'pos', columns = "sample", values = "haplotype").reset_index().rename_axis("", axis = 1)

output_file = open(f'{output}.txt', "w")
output_file.write(f'chr\tpos\thap1\thap2\tN10kbWin\tSS\n')

end = max(tableSpread['pos'])
list_haplotypes = list(tableSpread.columns)[1:]
NHap = len(list_haplotypes)

for start in range(0,end,step_size):
	tableSpreadWin = tableSpread[(tableSpread['pos'] >= start) & (tableSpread['pos'] < start+winSize)]
	SS_array = np.zeros((NHap, NHap))
	pairDistance = [row[:,None] != row for row in tableSpreadWin[list(tableSpreadWin.columns)[1:]].to_numpy()]
	for i in pairDistance:
		SS_array += i
	for i in range(NHap):
		for j in range(NHap):
			if j > i:
				output_file.write(f'{chr}\t{start}\t{list_haplotypes[i]}\t{list_haplotypes[j]}\t{len(pairDistance)}\t{SS_array[i,j]}\n')


```



python script: difference_perLargeWin_withGaps.py

```python
import os
import sys
import numpy as np
import pandas as pd

import scipy.cluster.hierarchy as shc
#from scipy.spatial.distance import pdist
from scipy.spatial.distance import cdist

import math

input = 'haplotypeGroups_min50SS_10KbWin_chr01.txt'
winSize=int(1000000)
step_size = int(50000)
output="test"
chr="chr01"

input = str(sys.argv[1])
winSize = int(sys.argv[2])
step_size = int(sys.argv[3])
output = str(sys.argv[4])
chr = str(sys.argv[5])

table = pd.read_table(input, sep="\t", names=['pos', 'sample', 'haplotype', 'missingP'])  
tableSpread = table[['pos', 'sample', 'haplotype']].pivot(index = 'pos', columns = "sample", values = "haplotype").reset_index().rename_axis("", axis = 1)
tableSpread_GAPS = table[['pos', 'sample', 'missingP']].pivot(index = 'pos', columns = "sample", values = "missingP").reset_index().rename_axis("", axis = 1)

output_file = open(f'{output}.txt', "w")
output_file.write(f'chr\tpos\thap1\thap2\tN10kbWin\tSS\n')

end = max(tableSpread['pos'])
list_haplotypes = list(tableSpread.columns)[1:]
NHap = len(list_haplotypes)
#list_gapHap = list(table[(table['missingP'] >= 0.9)]['haplotype'].unique())


def clean_dis(df1, df2):
	df1na = np.isnan(df1)
	df1clean = df1[~df1na]
	df2clean = df2[~df1na]
	df2na = np.isnan(df2clean)
	df1clean = df1clean[~df2na]
	df2clean = df2clean[~df2na]
	distance = (df1clean != df2clean).sum()
	return distance

for start in range(0,end,step_size):
	tableSpreadWin = tableSpread[(tableSpread['pos'] >= start) & (tableSpread['pos'] < start+winSize)]
	tableSpreadWin_GAPS = tableSpread_GAPS[(tableSpread_GAPS['pos'] >= start) & (tableSpread_GAPS['pos'] < start+winSize)]
	#gapsTable = np.isin(tableSpreadWin[list(tableSpreadWin.columns)[1:]].to_numpy(), list_gapHap)
	haplotypeTable = tableSpreadWin[list(tableSpreadWin.columns)[1:]].to_numpy().astype(float)
	haplotypeTable_GAPS = (tableSpreadWin_GAPS[list(tableSpreadWin_GAPS.columns)[1:]].to_numpy()) > 0.5
	haplotypeTable[haplotypeTable_GAPS] = np.nan
	#distance_haplotypeTable_tem = (cdist(haplotypeTable.T, haplotypeTable.T, 'hamming') - cdist(gapsTable.T, gapsTable.T, 'hamming')) * haplotypeTable.shape[0]
	#distance_haplotypeTable = (cdist(haplotypeTable.T, haplotypeTable.T, 'hamming') - cdist(haplotypeTable_GAPS.T, haplotypeTable_GAPS.T, 'hamming')) * haplotypeTable.shape[0]
	distance_haplotypeTable = (cdist(haplotypeTable.T, haplotypeTable.T, lambda u, v: clean_dis(u, v)))
	for i in range(NHap):
		for j in range(NHap):
			if j > i:
				output_file.write(f'{chr}\t{start}\t{list_haplotypes[i]}\t{list_haplotypes[j]}\t{haplotypeTable.shape[0]}\t{distance_haplotypeTable[i,j]}\n')



```



Bash script: difference_perLargeWin_withGaps.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=12:00:00
#SBATCH -J diffLargeWin


# input='/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/03_Nhaplotypes/haplotypeGroups_min10SS_chr06.txt'
# winSize=1000000
# step_size=50000
# output=test
# chr=chr06

input=$1
winSize=$2
step_size=$3
output=$4
chr=$5

# calculate pairwise difference between large windows:
python3 difference_perLargeWin_withGaps.py $input $winSize $step_size $output $chr

```



job:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/06_mergedHaplotypes_PanGenome

mamba activate env_others

# window size 1 Mb: step 50 kb
for i in $(seq -w 1 12)
do echo $i
sbatch difference_perLargeWin_withGaps.sh haplotypeGroups_min50SS_10KbWin_chr$i.txt 1000000 50000 NDiff_1Mb_step50kb_haplotypeGroups_min50SS_chr$i chr$i
done
# 

# window size 500 kb: step 50 kb
for i in $(seq -w 1 12)
do echo $i
sbatch difference_perLargeWin_withGaps.sh haplotypeGroups_min50SS_10KbWin_chr$i.txt 500000 50000 NDiff_500kb_step50kb_haplotypeGroups_min50SS_chr$i chr$i
done
# 

# window size 1 Mb: step 250 kb
for i in $(seq -w 1 12)
do echo $i
sbatch difference_perLargeWin_withGaps.sh haplotypeGroups_min50SS_10KbWin_chr$i.txt 1000000 250000 NDiff_1Mb_step250kb_haplotypeGroups_min50SS_chr$i chr$i
done
# 

# window size 500 kb: step 250 kb
for i in $(seq -w 1 12)
do echo $i
sbatch difference_perLargeWin_withGaps.sh haplotypeGroups_min50SS_10KbWin_chr$i.txt 500000 250000 NDiff_500kb_step250kb_haplotypeGroups_min50SS_chr$i chr$i
done
# 

# For all above: 06.04.2024: 2456484..2456531



# window size 1 Mb: step 50 kb
grep SS NDiff_1Mb_step50kb_haplotypeGroups_min50SS_chr01.txt > NDiff_1Mb_step50kb_haplotypeGroups_min50SS_allChr.txt
grep -hv SS NDiff_1Mb_step50kb_haplotypeGroups_min50SS_chr*.txt >> NDiff_1Mb_step50kb_haplotypeGroups_min50SS_allChr.txt

# window size 500 kb: step 50 kb
grep SS NDiff_500kb_step50kb_haplotypeGroups_min50SS_chr01.txt > NDiff_500kb_step50kb_haplotypeGroups_min50SS_allChr.txt
grep -hv SS NDiff_500kb_step50kb_haplotypeGroups_min50SS_chr*.txt >> NDiff_500kb_step50kb_haplotypeGroups_min50SS_allChr.txt

# window size 1 Mb: step 250 kb
grep SS NDiff_1Mb_step250kb_haplotypeGroups_min50SS_chr01.txt > NDiff_1Mb_step250kb_haplotypeGroups_min50SS_allChr.txt
grep -hv SS NDiff_1Mb_step250kb_haplotypeGroups_min50SS_chr*.txt >> NDiff_1Mb_step250kb_haplotypeGroups_min50SS_allChr.txt

# window size 500 kb: step 250 kb
grep SS NDiff_500kb_step250kb_haplotypeGroups_min50SS_chr01.txt > NDiff_500kb_step250kb_haplotypeGroups_min50SS_allChr.txt
grep -hv SS NDiff_500kb_step250kb_haplotypeGroups_min50SS_chr*.txt >> NDiff_500kb_step250kb_haplotypeGroups_min50SS_allChr.txt











# window size 1 Mb: step 50 kb
for i in $(seq -w 1 12)
do echo $i
sbatch difference_perLargeWin_withGaps.sh haplotypeGroups_min20SS_10KbWin_chr$i.txt 1000000 50000 NDiff_1Mb_step50kb_haplotypeGroups_min20SS_chr$i chr$i
done
# 

# window size 500 kb: step 50 kb
for i in $(seq -w 1 12)
do echo $i
sbatch difference_perLargeWin_withGaps.sh haplotypeGroups_min20SS_10KbWin_chr$i.txt 500000 50000 NDiff_500kb_step50kb_haplotypeGroups_min20SS_chr$i chr$i
done
# 

# window size 1 Mb: step 250 kb
for i in $(seq -w 1 12)
do echo $i
sbatch difference_perLargeWin_withGaps.sh haplotypeGroups_min20SS_10KbWin_chr$i.txt 1000000 250000 NDiff_1Mb_step250kb_haplotypeGroups_min20SS_chr$i chr$i
done
# 

# window size 500 kb: step 250 kb
for i in $(seq -w 1 12)
do echo $i
sbatch difference_perLargeWin_withGaps.sh haplotypeGroups_min20SS_10KbWin_chr$i.txt 500000 250000 NDiff_500kb_step250kb_haplotypeGroups_min20SS_chr$i chr$i
done
# 06.04.2024: with clustering change 2456537..2456584



# window size 1 Mb: step 50 kb
grep SS NDiff_1Mb_step50kb_haplotypeGroups_min20SS_chr01.txt > NDiff_1Mb_step50kb_haplotypeGroups_min20SS_allChr.txt
grep -hv SS NDiff_1Mb_step50kb_haplotypeGroups_min20SS_chr*.txt >> NDiff_1Mb_step50kb_haplotypeGroups_min20SS_allChr.txt

# window size 500 kb: step 50 kb
grep SS NDiff_500kb_step50kb_haplotypeGroups_min20SS_chr01.txt > NDiff_500kb_step50kb_haplotypeGroups_min20SS_allChr.txt
grep -hv SS NDiff_500kb_step50kb_haplotypeGroups_min20SS_chr*.txt >> NDiff_500kb_step50kb_haplotypeGroups_min20SS_allChr.txt

# window size 1 Mb: step 250 kb
grep SS NDiff_1Mb_step250kb_haplotypeGroups_min20SS_chr01.txt > NDiff_1Mb_step250kb_haplotypeGroups_min20SS_allChr.txt
grep -hv SS NDiff_1Mb_step250kb_haplotypeGroups_min20SS_chr*.txt >> NDiff_1Mb_step250kb_haplotypeGroups_min20SS_allChr.txt

# window size 500 kb: step 250 kb
grep SS NDiff_500kb_step250kb_haplotypeGroups_min20SS_chr01.txt > NDiff_500kb_step250kb_haplotypeGroups_min20SS_allChr.txt
grep -hv SS NDiff_500kb_step250kb_haplotypeGroups_min20SS_chr*.txt >> NDiff_500kb_step250kb_haplotypeGroups_min20SS_allChr.txt

# done

```





### Plot Distribution of Difference between Large Windows:



```R
# .libPaths("/dss/dsshome1/lxc03/di36guz2/R")
# mamba activate env_others
.libPaths()

setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/06_mergedHaplotypes_PanGenome")

library(ggplot2)
library(tidyverse)


minSS <- "50SS"
minSS <- "20SS"

# winSize<-"1Mb"
# step<-"50kb"
# # 
# winSize<-"1Mb"
# step<-"250kb"
# #
# winSize<-"500kb"
# step<-"50kb"
# #
# winSize<-"500kb"
# step<-"250kb"
# 


winSize_list <- c("1Mb", "1Mb", "500kb", "500kb")
step_list <- c("50kb", "250kb", "50kb", "250kb")

for (index in seq(1,4)){
  print(index)
  winSize<-winSize_list[index]
  step<-step_list[index]
  
  table_hapdiff <- read_table(paste0("NDiff_", winSize, "_step", step, "_haplotypeGroups_min", minSS, "_allChr.txt"), T)
  
  table_hapdiff %>% 
    ggplot(aes(SS))+
    geom_histogram(fill="black") +
    labs(x = "SS",
         y = "Count") +
    facet_grid(. ~ chr)+
    theme_classic()+
    theme(panel.grid.major = element_line(color = "grey80"),
          panel.grid.minor = element_line(color = "grey80"))
  
  ggsave(paste0("01_NumberDiff_", winSize, "_step", step, "_", minSS, ".png"), height=7, width=14)
  ggsave(paste0("01_NumberDiff_", winSize, "_step", step, "_", minSS, ".pdf"), height=7, width=14)
  
  
  table_hapdiff %>% 
    mutate(frac=SS/N10kbWin) %>%
    ggplot(aes(frac))+
    geom_histogram(fill="black") +
    labs(x = "frac SS",
         y = "Count") +
    facet_grid(. ~ chr)+
    theme_classic()+
    theme(panel.grid.major = element_line(color = "grey80"),
          panel.grid.minor = element_line(color = "grey80"))
  
  ggsave(paste0("02_FracDiff_", winSize, "_step", step, "_", minSS, ".png"), height=7, width=14)
  ggsave(paste0("02_FracDiff_", winSize, "_step", step, "_", minSS, ".pdf"), height=7, width=14)
}




library(ggplot2)
library(dplyr)
library(hrbrthemes)

# winSize_list <- c("1Mb", "1Mb", "500kb", "500kb")
# step_list <- c("50kb", "250kb", "50kb", "250kb")


winSize<-"1Mb"
step<-"50kb"

winSize<-"1Mb"
step<-"250kb"

winSize<-"500kb"
step<-"50kb"

winSize<-"500kb"
step<-"250kb"


minSS <- "50SS"
table_hapdiff <- read_table(paste0("NDiff_", winSize, "_step", step, "_haplotypeGroups_min", minSS, "_allChr.txt"), T) %>%
  mutate(minSS="50SS")

minSS <- "20SS"
table_hapdiff2 <- read_table(paste0("NDiff_", winSize, "_step", step, "_haplotypeGroups_min", minSS, "_allChr.txt"), T) %>%
  mutate(minSS="20SS")


table_hapdiff %>% 
  rbind(table_hapdiff2) %>%
  # filter(SS<0)
  ggplot(aes(SS, fill=minSS))+
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  # scale_fill_manual(values=c("#69b3a2", "#404080")) +
  scale_fill_manual(values=c("blue", "red")) +
  scale_x_continuous(expand = c(0,0))+
  # geom_histogram() +
  # geom_histogram(position="dodge") +
  labs(x = "SS",
       y = "Count") +
  # facet_grid(. ~ chr)+
  facet_grid(chr ~ .)+
  theme_classic()+
  # theme_ipsum() +
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"))

ggsave(paste0("03_NumberDiff_", winSize, "_step", step, ".png"), height=7, width=14)
ggsave(paste0("03_NumberDiff_", winSize, "_step", step, ".pdf"), height=7, width=14)


```





## Clustering of large haplotypes



script: SecondClusteringLargeHaplotypes_overlappingWin.py

```python
import os
import sys
import numpy as np
import pandas as pd

import scipy.cluster.hierarchy as shc
from scipy.spatial.distance import pdist, squareform
import math

# input = 'haplotypeGroups_min50SS_10KbWin_chr06.txt'
# winSize=int(1000000)
# step_size = int(50000)
# output="test"
# chr="chr06"
# pro_max_d = float(0.25)

input = str(sys.argv[1])
winSize = int(sys.argv[2])
step_size = int(sys.argv[3])
output = str(sys.argv[4])
chr = str(sys.argv[5])
pro_max_d = float(sys.argv[6])

table = pd.read_table(input, sep="\t", names=['pos', 'sample', 'haplotype', 'missingP'])  
tableSpread = table[['pos', 'sample', 'haplotype']].pivot(index = 'pos', columns = "sample", values = "haplotype").reset_index().rename_axis("", axis = 1)
tableSpread_missingP = table[['pos', 'sample', 'missingP']].pivot(index = 'pos', columns = "sample", values = "missingP").reset_index().rename_axis("", axis = 1)

output_file = open(f'{output}.txt', "w")
output_file.write(f'chr\tpos\tsample\thaplotype\tMean_NAN_haplo\n')
output_file.close()

end = max(tableSpread['pos'])
list_haplotypes = list(tableSpread.columns)[1:]
NHap = len(list_haplotypes)

max_d = int((winSize/10000)*pro_max_d)
table_hap_groups_pre = pd.DataFrame()


for start in range(0,end,step_size):
  tableSpreadWin = tableSpread[(tableSpread['pos'] >= start) & (tableSpread['pos'] < start+winSize)]
  tableSpreadWin_missingP = tableSpread_missingP[(tableSpread_missingP['pos'] >= start) & (tableSpread_missingP['pos'] < start+winSize)]
  total_missingP_win = [x.mean() for x in tableSpreadWin_missingP[list(tableSpreadWin_missingP.columns)[1:]].to_numpy().T]
  SS_array = np.zeros((NHap, NHap))
  pairDistance = [row[:,None] != row for row in tableSpreadWin[list(tableSpreadWin.columns)[1:]].to_numpy()]
  for i in pairDistance:
    SS_array += i
  condensed_dist = squareform(SS_array)
  linkresult = shc.linkage(condensed_dist, method='complete')
  clusters = shc.fcluster(linkresult, max_d, criterion='distance')
  if not table_hap_groups_pre.empty:
    max_cluster_pre = table_hap_groups_pre['haplotype_group_pre'].max()
    table_hap_groups = pd.DataFrame(data={'labs': list(list_haplotypes), 'haplotype_group': list(clusters+max_cluster_pre), 'Mean_NAN_haplo': list(total_missingP_win)})
    comparisonDF = table_hap_groups.merge(table_hap_groups_pre, on='labs', how='left').dropna().groupby(['haplotype_group', 'haplotype_group_pre']).agg({'count'})['labs'].reset_index()
    idx = comparisonDF.groupby('haplotype_group')['count'].transform(max) == comparisonDF['count']
    idx2 = comparisonDF[idx].groupby('haplotype_group_pre')['count'].transform(max) == comparisonDF[idx]['count']
    previous_hap_groups = comparisonDF[idx][idx2].groupby('haplotype_group_pre').first().reset_index().groupby('haplotype_group').first().reset_index()
    table_hap_groups_merge = table_hap_groups.merge( previous_hap_groups[['haplotype_group', 'haplotype_group_pre']], on='haplotype_group', how='left')
    table_hap_groups_merge['haplotype_group_ed'] = [int(table_hap_groups_merge['haplotype_group'][x]) if math.isnan(table_hap_groups_merge['haplotype_group_pre'][x]) else int(table_hap_groups_merge['haplotype_group_pre'][x]) for x in range(0, table_hap_groups_merge.shape[0])]
    table_hap_groups_merge['chr'] = chr
    table_hap_groups_merge['pos'] = start
    table_hap_groups_merge[['chr', 'pos', 'labs', 'haplotype_group_ed', 'Mean_NAN_haplo']].to_csv(f"{output}.txt", sep="\t", index=False, header=False, mode='a')
    table_hap_groups_pre = table_hap_groups_merge[['labs', 'haplotype_group_ed']].rename(columns={"labs": "labs", "haplotype_group_ed": "haplotype_group_pre"})
  else:
    table_hap_groups = pd.DataFrame(data={'labs': list(list_haplotypes), 'haplotype_group': list(clusters), 'Mean_NAN_haplo': list(total_missingP_win)})
    table_hap_groups['chr'] = chr
    table_hap_groups['pos'] = start
    table_hap_groups[['chr', 'pos', 'labs', 'haplotype_group', 'Mean_NAN_haplo']].to_csv(f"{output}.txt", sep="\t", index=False, header=False, mode='a')
    table_hap_groups_pre = table_hap_groups[['labs', 'haplotype_group']].rename(columns={"labs": "labs", "haplotype_group": "haplotype_group_pre"})




```



script removing gaps: SecondClusteringLargeHaplotypes_overlappingWin.py

```python
import os
import sys
import numpy as np
import pandas as pd

import scipy.cluster.hierarchy as shc
from scipy.spatial.distance import pdist, squareform
import math

# input = 'haplotypeGroups_min50SS_10KbWin_chr07.txt'
# winSize=int(1000000)
# step_size = int(50000)
# output="test"
# chr="chr07"
# pro_max_d = float(0.25)

input = str(sys.argv[1])
winSize = int(sys.argv[2])
step_size = int(sys.argv[3])
output = str(sys.argv[4])
chr = str(sys.argv[5])
pro_max_d = float(sys.argv[6])

table = pd.read_table(input, sep="\t", names=['pos', 'sample', 'haplotype', 'missingP'])  
tableSpread = table[['pos', 'sample', 'haplotype']].pivot(index = 'pos', columns = "sample", values = "haplotype").reset_index().rename_axis("", axis = 1)
tableSpread_missingP = table[['pos', 'sample', 'missingP']].pivot(index = 'pos', columns = "sample", values = "missingP").reset_index().rename_axis("", axis = 1)

output_file = open(f'{output}.txt', "w")
output_file.write(f'chr\tpos\tsample\thaplotype\tMean_NAN_haplo\n')
output_file.close()

end = max(tableSpread['pos'])
list_haplotypes = list(tableSpread.columns)[1:]
NHap = len(list_haplotypes)

max_d = int((winSize/10000)*pro_max_d)
table_hap_groups_pre = pd.DataFrame()


for start in range(0,end,step_size):
  tableSpreadWin = tableSpread[(tableSpread['pos'] >= start) & (tableSpread['pos'] < start+winSize)]
  tableSpreadWin_missingP = tableSpread_missingP[(tableSpread_missingP['pos'] >= start) & (tableSpread_missingP['pos'] < start+winSize)]
  if tableSpreadWin.shape[0] > 0:
    total_missingP_win = [x.mean() for x in tableSpreadWin_missingP[list(tableSpreadWin_missingP.columns)[1:]].to_numpy().T]
    SS_array = np.zeros((NHap, NHap))
    tableSpreadWin_missingP_np = tableSpreadWin_missingP[list(tableSpreadWin.columns)[1:]].to_numpy()>0.9
    pairDistance = [row[:,None] != row for row in tableSpreadWin[list(tableSpreadWin.columns)[1:]].to_numpy()]
    for i in range(len(pairDistance)):
      new_dis = pairDistance[i]
      new_dis[tableSpreadWin_missingP_np[i], :] = False
      new_dis[:, tableSpreadWin_missingP_np[i]] = False
      SS_array += new_dis
    condensed_dist = squareform(SS_array)
    linkresult = shc.linkage(condensed_dist, method='complete')
    clusters = shc.fcluster(linkresult, max_d, criterion='distance')
    if not table_hap_groups_pre.empty:
      max_cluster_pre = table_hap_groups_pre['haplotype_group_pre'].max()
      table_hap_groups = pd.DataFrame(data={'labs': list(list_haplotypes), 'haplotype_group': list(clusters+max_cluster_pre), 'Mean_NAN_haplo': list(total_missingP_win)})
      comparisonDF = table_hap_groups.merge(table_hap_groups_pre, on='labs', how='left').dropna().groupby(['haplotype_group', 'haplotype_group_pre']).agg({'count'})['labs'].reset_index()
      idx = comparisonDF.groupby('haplotype_group')['count'].transform(max) == comparisonDF['count']
      idx2 = comparisonDF[idx].groupby('haplotype_group_pre')['count'].transform(max) == comparisonDF[idx]['count']
      previous_hap_groups = comparisonDF[idx][idx2].groupby('haplotype_group_pre').first().reset_index().groupby('haplotype_group').first().reset_index()
      table_hap_groups_merge = table_hap_groups.merge( previous_hap_groups[['haplotype_group', 'haplotype_group_pre']], on='haplotype_group', how='left')
      table_hap_groups_merge['haplotype_group_ed'] = [int(table_hap_groups_merge['haplotype_group'][x]) if math.isnan(table_hap_groups_merge['haplotype_group_pre'][x]) else int(table_hap_groups_merge['haplotype_group_pre'][x]) for x in range(0, table_hap_groups_merge.shape[0])]
      table_hap_groups_merge['chr'] = chr
      table_hap_groups_merge['pos'] = start
      table_hap_groups_merge[['chr', 'pos', 'labs', 'haplotype_group_ed', 'Mean_NAN_haplo']].to_csv(f"{output}.txt", sep="\t", index=False, header=False, mode='a')
      table_hap_groups_pre = table_hap_groups_merge[['labs', 'haplotype_group_ed']].rename(columns={"labs": "labs", "haplotype_group_ed": "haplotype_group_pre"})
    else:
      table_hap_groups = pd.DataFrame(data={'labs': list(list_haplotypes), 'haplotype_group': list(clusters), 'Mean_NAN_haplo': list(total_missingP_win)})
      table_hap_groups['chr'] = chr
      table_hap_groups['pos'] = start
      table_hap_groups[['chr', 'pos', 'labs', 'haplotype_group', 'Mean_NAN_haplo']].to_csv(f"{output}.txt", sep="\t", index=False, header=False, mode='a')
      table_hap_groups_pre = table_hap_groups[['labs', 'haplotype_group']].rename(columns={"labs": "labs", "haplotype_group": "haplotype_group_pre"})
  else:
    table_hap_groups['chr'] = chr
    table_hap_groups['pos'] = start
    table_hap_groups[['chr', 'pos', 'labs', 'haplotype_group', 'Mean_NAN_haplo']].to_csv(f"{output}.txt", sep="\t", index=False, header=False, mode='a')



```



Bash script: SecondClusteringLargeHaplotypes_overlappingWin.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=12:00:00
#SBATCH -J Second_Clustering

# input = 'haplotypeGroups_min50SS_10KbWin_chr06.txt'
# winSize=int(1000000)
# step_size = int(50000)
# output="test"
# chr="chr06"
# pro_max_d = float(0.25)

input=$1
winSize=$2
step_size=$3
output=$4
chr=$5
pro_max_d=$6

# calculate pairwise difference between large windows:
python3 SecondClusteringLargeHaplotypes_overlappingWin.py $input $winSize $step_size $output $chr $pro_max_d


```



job:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/06_mergedHaplotypes_PanGenome

mamba activate env_others

# window size 1 Mb: step 50 kb
for i in $(seq -w 1 12)
do echo $i
sbatch SecondClusteringLargeHaplotypes_overlappingWin.sh haplotypeGroups_min50SS_10KbWin_chr$i.txt 1000000 50000 clusteringLargeHaplotypes_1Mb_step50kb_haplotypeGroups_min50SS_chr$i chr$i 0.25
done
# 2223791..2223802

# window size 500kb: step 50 kb
for i in $(seq -w 1 12)
do echo $i
sbatch SecondClusteringLargeHaplotypes_overlappingWin.sh haplotypeGroups_min50SS_10KbWin_chr$i.txt 500000 50000 clusteringLargeHaplotypes_500kb_step50kb_haplotypeGroups_min50SS_chr$i chr$i 0.25
done
# 2223969..2223980

# window size 250kb: step 50 kb
for i in $(seq -w 1 12)
do echo $i
sbatch SecondClusteringLargeHaplotypes_overlappingWin.sh haplotypeGroups_min50SS_10KbWin_chr$i.txt 250000 50000 clusteringLargeHaplotypes_250kb_step50kb_haplotypeGroups_min50SS_chr$i chr$i 0.25
done
# 2223983..2223994

# window size 200kb: step 50 kb
for i in $(seq -w 1 12)
do echo $i
sbatch SecondClusteringLargeHaplotypes_overlappingWin.sh haplotypeGroups_min50SS_10KbWin_chr$i.txt 200000 50000 clusteringLargeHaplotypes_200kb_step50kb_haplotypeGroups_min50SS_chr$i chr$i 0.25
done
# old: 23.11:2244696..2244707
# 31.01.2024: 2313149..2313160



# window size 200kb: step 100kb
for i in $(seq -w 1 12)
do echo $i
sbatch SecondClusteringLargeHaplotypes_overlappingWin.sh haplotypeGroups_min50SS_10KbWin_chr$i.txt 200000 100000 clusteringLargeHaplotypes_200kb_step100kb_haplotypeGroups_min50SS_chr$i chr$i 0.25
done
# 2225965..2225976

# 06.04.2024: 2456586..2456586
# for all the above


chr=chr06 ; awk '{print $2"\t"$4}' clusteringLargeHaplotypes_500kb_step50kb_haplotypeGroups_min50SS_$chr.txt | grep -v pos | sort | uniq | awk '{print $1}' | sort | uniq -c | sort -n -k 2 > tem_500_$chr.txt ; awk '{print $2"\t"$4}' clusteringLargeHaplotypes_1Mb_step50kb_haplotypeGroups_min50SS_$chr.txt | grep -v pos | sort | uniq | awk '{print $1}' | sort | uniq -c | sort -n -k 2 > tem_1000_$chr.txt ; awk '{print $2"\t"$4}' clusteringLargeHaplotypes_250kb_step50kb_haplotypeGroups_min50SS_$chr.txt | grep -v pos | sort | uniq | awk '{print $1}' | sort | uniq -c | sort -n -k 2 > tem_250_$chr.txt ; awk '{print $2"\t"$4}' clusteringLargeHaplotypes_200kb_step50kb_haplotypeGroups_min50SS_$chr.txt | grep -v pos | sort | uniq | awk '{print $1}' | sort | uniq -c | sort -n -k 2 > tem_200_$chr.txt ; awk '{print $2"\t"$4}' clusteringLargeHaplotypes_200kb_step100kb_haplotypeGroups_min50SS_$chr.txt | grep -v pos | sort | uniq | awk '{print $1}' | sort | uniq -c | sort -n -k 2 > tem_200_100_$chr.txt ; paste tem_1000_$chr.txt tem_500_$chr.txt tem_250_$chr.txt tem_200_$chr.txt tem_200_100_$chr.txt | less -S 

```



script: NormaliseHapNums_to40_withOffset.py

```python
#Script to take one of Sergio's haplotypeGroup files and
import sys

inFile = sys.argv[1] #/Users/craig/Documents/Sergio_HapVis/haplotypeGroups_min10SS_chr06
outFile= inFile+"_normed40_withOffset.tsv"
out=open(outFile,"w+")
#windowSize=10000 #update this
current_start="0"
HapList=[] #List of haplotypes seen in a given window
Lines=[] #Store lines of the file for one window.
hapOffset=0.02 #An offset added to each subsequent member of the same haplotype cluster, so that the lines in R don't all draw on top of eachother
maxGap = float(sys.argv[2])

#Genomes= ["A_hap21","A_hap22","A_hap23","A_hap24","B_hap21","B_hap22","B_hap23","B_hap24","C_hap21","C_hap22","C_hap23","C_hap24","D_hap21","D_hap22","D_hap23","D_hap24","E_hap21","E_hap22","E_hap23","E_hap24","F_hap21","F_hap22","F_hap23","F_hap24","G_hap21","G_hap22","G_hap23","G_hap24","H_hap21","H_hap22","H_hap23","H_hap24","I_hap21","I_hap22","I_hap23","I_hap24","J_hap21","J_hap22","J_hap23","J_hap24"]

#Go through lines of the haplotypeGroup file
for line in open(inFile+".txt","r"):

        chr, window_start, genome, hap, missingFraction = line.split("\t")

        #If we have some new window, output the old stuff
        if window_start != current_start:
                Sorted_haps = sorted(set(HapList)) #sort the list of haplotype groups observed for this window
                seenHaps={} #Dictionary which will record the offset needed for drawing lines later, given how often we've seen a particular haplotype group in a window.

                for l in Lines: #for each line seen in this window
                        tmp_chr, tmp_window_start, tmp_genome, tmp_hap, tmp_missingFraction = l.rstrip().split("\t")

                        if float(tmp_missingFraction)>=0.9: # assign a -1 hap to missing values
                                tmp_hap="-1"
                        else: # otherwise replace the hap number with the ordinal position of the haplotype among observed haplotypes
                                tmp_hap = str(Sorted_haps.index(tmp_hap)+1) #eg if observed haplotypes are [3,6,24], then 3->1, 6->2, 24->3 respectively

                        #Update the list of haplotypes we have seen.
                        if tmp_hap not in seenHaps:
                                seenHaps[tmp_hap]=hapOffset # If this is the first time we've seen this haplotype, add the offset ready for the next one.
                        else: #if we have seen this before, increase the offset
                                new_tmp_hap = str(float(tmp_hap) + seenHaps[tmp_hap])
                                seenHaps[tmp_hap]=seenHaps[tmp_hap]+hapOffset
                                tmp_hap=new_tmp_hap

                        #output
                        out.write("\t".join([tmp_chr, tmp_window_start, tmp_genome, tmp_hap, tmp_missingFraction])+"\n")

                        #Reset some vars
                        current_start=window_start
                        HapList=[]
                        Lines=[]
        # else, if we are still getting info for this same window
        if window_start == current_start:

                Lines.append(line) #add this line to the list of lines seen.
                if float(missingFraction)<0.9: #only add haps that aren't missing
                        HapList.append(hap) #build the list of haplotypes



#catch the last window here
Sorted_haps = sorted(HapList)
for l in Lines:
        tmp_chr, tmp_window_start, tmp_genome, tmp_hap, tmp_missingFraction = l.rstrip().split("\t")
        if float(tmp_missingFraction)>=0.9: # assign a -1 hap to missing values
                tmp_hap="-1"
        else: # replace the hap nurmber with the ordinal position of the haplotype among observed haplotypes
                tmp_hap = str(Sorted_haps.index(tmp_hap)+1)
        out.write("\t".join([tmp_chr, tmp_window_start, tmp_genome, tmp_hap, tmp_missingFraction])+"\n")


```



produce normalized tables:

```bash
mamba activate env_others

# maxGap = 0.9
for i in $(seq -w 1 12)
do echo $i
grep -v pos clusteringLargeHaplotypes_200kb_step50kb_haplotypeGroups_min50SS_chr$i.txt > tem_clusteringLargeHaplotypes_200kb_step50kb_haplotypeGroups_min50SS_chr$i.txt
python3 NormaliseHapNums_to40_withOffset.py tem_clusteringLargeHaplotypes_200kb_step50kb_haplotypeGroups_min50SS_chr$i 0.9
mv tem_clusteringLargeHaplotypes_200kb_step50kb_haplotypeGroups_min50SS_chr$i"_"normed40_withOffset.tsv clusteringLargeHaplotypes_200kb_step50kb_haplotypeGroups_min50SS_chr$i"_"normed40_withOffset.tsv 
rm tem_clusteringLargeHaplotypes_200kb_step50kb_haplotypeGroups_min50SS_chr$i.txt
done
# 


# maxGap = 0.6
for i in $(seq -w 1 12)
do echo $i
grep -v pos clusteringLargeHaplotypes_200kb_step50kb_haplotypeGroups_min50SS_chr$i.txt > tem_clusteringLargeHaplotypes_200kb_step50kb_haplotypeGroups_min50SS_chr$i.txt
python3 NormaliseHapNums_to40_withOffset.py tem_clusteringLargeHaplotypes_200kb_step50kb_haplotypeGroups_min50SS_chr$i 0.6
mv tem_clusteringLargeHaplotypes_200kb_step50kb_haplotypeGroups_min50SS_chr$i"_"normed40_withOffset.tsv clusteringLargeHaplotypes_200kb_step50kb_haplotypeGroups_min50SS_chr$i"_MaxGap06_"normed40_withOffset.tsv 
rm tem_clusteringLargeHaplotypes_200kb_step50kb_haplotypeGroups_min50SS_chr$i.txt
done
# 

```



### Plots:



```R
# .libPaths("/dss/dsshome1/lxc03/di36guz2/R")
# mamba activate env_others
.libPaths()

setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/06_mergedHaplotypes_PanGenome")

library(ggplot2)
library(tidyverse)

library(dplyr)
library(stringi)
options(scipen = 999) 

winSize="200kb"
stepSize="50kb"
minSS="50SS"
chrID="01"
MaxGap="06"


chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")


for (chrID in chr_list) {
  #Read in haplotypeNumberFile, with haplotype numbers reduced by the python script. 
  
  hapNums.in <- read.table(paste0("clusteringLargeHaplotypes_", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_chr", chrID, 
                                  "_MaxGap", MaxGap, "_normed40_withOffset.tsv"), 
                           header=F, sep="\t")
  
  names(hapNums.in) = c("chr", "pos","Genome", "HapPos", "MissingFraction")
  
  hapNums.in %>%
    mutate(pos=pos/1000000, 
           posEnd=pos+(50000/1000000)) %>%
    ggplot()+
    geom_segment(aes(x=pos,xend=posEnd, y=HapPos,yend=HapPos, colour=Genome)) +
    scale_x_continuous(breaks=seq(0,100,2), expand = c(0,0))+
    scale_y_continuous(breaks=seq(0,40,1))+
    labs(title=paste0("Chr", chrID),
         x ="Pos (Mb)", 
         y = "Hap Sample")+
    theme_classic()+
    theme(legend.position = "none",
          # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          panel.grid.minor = element_blank())
  
  
  ggsave(filename=paste0("04_clusteringLargeHaplotypes_chr", chrID, "_Win", winSize, "_step", stepSize, 
                         "_haplotypeGroups_min", minSS, "_MaxGap", MaxGap, ".png"), 
         device = png, units="cm", 
         width=25, height=14)
  ggsave(filename=paste0("04_clusteringLargeHaplotypes_chr", chrID, "_Win", winSize, "_step", stepSize, 
                         "_haplotypeGroups_min", minSS, "_MaxGap", MaxGap, ".pdf"), 
         device = pdf, units="cm", 
         width=25, height=14)
  
}



# Particular examples:

winSize="200kb"
stepSize="50kb"
minSS="50SS"
chrID="06"
MaxGap="06"
StartSet=10000000
EndSet=30000000

hapNums.in <- read.table(paste0("clusteringLargeHaplotypes_", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_chr", chrID, 
                                "_MaxGap", MaxGap, "_normed40_withOffset.tsv"), 
                         header=F, sep="\t")

names(hapNums.in) = c("chr", "pos","Genome", "HapPos", "MissingFraction")

hapNums.in %>%
  filter(pos>StartSet) %>%
  filter(pos<EndSet) %>%
  mutate(pos=pos/1000000, 
         posEnd=pos+(50000/1000000)) %>%
  ggplot()+
  geom_segment(aes(x=pos,xend=posEnd, y=HapPos,yend=HapPos, colour=Genome)) +
  scale_x_continuous(breaks=seq(0,100,2), expand = c(0,0))+
  scale_y_continuous(breaks=seq(0,40,1))+
  labs(title=paste0("Chr", chrID),
       x ="Pos (Mb)", 
       y = "Hap Sample")+
  theme_classic()+
  theme(legend.position = "none",
        # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.grid.minor = element_blank())


ggsave(filename=paste0("05_example_clusteringLargeHaplotypes_chr", chrID, "_Win", winSize, "_step", stepSize, 
                       "_haplotypeGroups_min", minSS, "_MaxGap", MaxGap, ".png"), 
       device = png, units="cm", 
       width=25, height=14)
ggsave(filename=paste0("05_example_clusteringLargeHaplotypes_chr", chrID, "_Win", winSize, "_step", stepSize, 
                       "_haplotypeGroups_min", minSS, "_MaxGap", MaxGap, ".pdf"), 
       device = pdf, units="cm", 
       width=25, height=14)










winSize="200kb"
stepSize="50kb"
minSS="50SS"
chrID="01"
MaxGap="06"
StartSet=50000000
EndSet=65000000

hapNums.in <- read.table(paste0("clusteringLargeHaplotypes_", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_chr", chrID, 
                                "_MaxGap", MaxGap, "_normed40_withOffset.tsv"), 
                         header=F, sep="\t")

names(hapNums.in) = c("chr", "pos","Genome", "HapPos", "MissingFraction")

hapNums.in %>%
  filter(pos>StartSet) %>%
  filter(pos<EndSet) %>%
  mutate(pos=pos/1000000, 
         posEnd=pos+(50000/1000000)) %>%
  ggplot()+
  geom_segment(aes(x=pos,xend=posEnd, y=HapPos,yend=HapPos, colour=Genome)) +
  scale_x_continuous(breaks=seq(0,100,2), expand = c(0,0))+
  scale_y_continuous(breaks=seq(0,40,1))+
  labs(title=paste0("Chr", chrID),
       x ="Pos (Mb)", 
       y = "Hap Sample")+
  theme_classic()+
  theme(legend.position = "none",
        # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.grid.minor = element_blank())


ggsave(filename=paste0("05_example_clusteringLargeHaplotypes_chr", chrID, "_Win", winSize, "_step", stepSize, 
                       "_haplotypeGroups_min", minSS, "_MaxGap", MaxGap, ".png"), 
       device = png, units="cm", 
       width=25, height=14)
ggsave(filename=paste0("05_example_clusteringLargeHaplotypes_chr", chrID, "_Win", winSize, "_step", stepSize, 
                       "_haplotypeGroups_min", minSS, "_MaxGap", MaxGap, ".pdf"), 
       device = pdf, units="cm", 
       width=25, height=14)














winSize="200kb"
stepSize="50kb"
minSS="50SS"
chrID="10"
MaxGap="06"
StartSet=32000000
EndSet=45000000

hapNums.in <- read.table(paste0("clusteringLargeHaplotypes_", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_chr", chrID, 
                                "_MaxGap", MaxGap, "_normed40_withOffset.tsv"), 
                         header=F, sep="\t")

names(hapNums.in) = c("chr", "pos","Genome", "HapPos", "MissingFraction")

hapNums.in %>%
  filter(pos>StartSet) %>%
  filter(pos<EndSet) %>%
  mutate(pos=pos/1000000, 
         posEnd=pos+(50000/1000000)) %>%
  ggplot()+
  geom_segment(aes(x=pos,xend=posEnd, y=HapPos,yend=HapPos, colour=Genome)) +
  scale_x_continuous(breaks=seq(0,100,2), expand = c(0,0))+
  scale_y_continuous(breaks=seq(0,40,1))+
  labs(title=paste0("Chr", chrID),
       x ="Pos (Mb)", 
       y = "Hap Sample")+
  theme_classic()+
  theme(legend.position = "none",
        # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.grid.minor = element_blank())


ggsave(filename=paste0("05_example_clusteringLargeHaplotypes_chr", chrID, "_Win", winSize, "_step", stepSize, 
                       "_haplotypeGroups_min", minSS, "_MaxGap", MaxGap, ".png"), 
       device = png, units="cm", 
       width=25, height=14)
ggsave(filename=paste0("05_example_clusteringLargeHaplotypes_chr", chrID, "_Win", winSize, "_step", stepSize, 
                       "_haplotypeGroups_min", minSS, "_MaxGap", MaxGap, ".pdf"), 
       device = pdf, units="cm", 
       width=25, height=14)









# haplotype IDs


winSize="200kb"
stepSize="50kb"
minSS="50SS"
chrID="10"
MaxGap="06"

chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")


for (chrID in chr_list) {
  #Read in haplotypeNumberFile, with haplotype numbers reduced by the python script. 
  hapNums.in <- read.table(paste0("clusteringLargeHaplotypes_", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_chr", chrID, 
                                  ".txt"), 
                           header=T, sep="\t")
  n=10
  hapNums.in %>%
    filter(Mean_NAN_haplo<0.6) %>%
    ggplot(aes(pos/1000000, sample, fill=factor(haplotype) , colour=factor(haplotype))) +
    geom_tile()+
    scale_color_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
                                      hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
    scale_fill_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
                                     hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
    scale_x_continuous(breaks=seq(0,100,5), name = "Position (Mb)", expand = c(0,0))+
    ggtitle(paste0("Chr:", chrID))+
    theme_classic()+
    theme_classic()+
    theme(panel.grid.major = element_line(color = "grey80"),
          panel.grid.minor = element_line(color = "grey80"),
          axis.text = element_text(size=18),
          axis.title = element_text(size=18),
          legend.position = "none")
  
  ggsave(filename=paste0("06_hapDis_clusteringLargeHaplotypes_chr", chrID, "_Win", winSize, "_step", stepSize, 
                         "_haplotypeGroups_min", minSS, ".png"), 
         device = png, units="cm", 
         width=25, height=18)
  ggsave(filename=paste0("06_hapDis_clusteringLargeHaplotypes_chr", chrID, "_Win", winSize, "_step", stepSize, 
                         "_haplotypeGroups_min", minSS, ".pdf"), 
         device = pdf, units="cm", 
         width=25, height=18)
}










# examples with haplotype colours:
winSize="200kb"
stepSize="50kb"
minSS="50SS"
chrID="01"
MaxGap="06"
StartSet=50000000
EndSet=65000000


#Read in haplotypeNumberFile, with haplotype numbers reduced by the python script. 
hapNums.in <- read.table(paste0("clusteringLargeHaplotypes_", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_chr", chrID, 
                                ".txt"), 
                         header=T, sep="\t")
n=10
hapNums.in %>%
  filter(pos>StartSet) %>%
  filter(pos<EndSet) %>%
  filter(Mean_NAN_haplo<MaxGap) %>%
  ggplot(aes(pos/1000000, sample, fill=factor(haplotype) , colour=factor(haplotype))) +
  geom_tile()+
  scale_color_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
                                    hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
  scale_fill_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
                                   hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
  scale_x_continuous(breaks=seq(0,100,2), name = "Position (Mb)", expand = c(0,0))+
  ggtitle(paste0("Chr:", chrID))+
  theme_classic()+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"),
        axis.text = element_text(size=18),
        axis.title = element_text(size=18),
        legend.position = "none")

ggsave(filename=paste0("07_example_hapDis_clusteringLargeHaplotypes_chr", chrID, "_Win", winSize, "_step", stepSize, 
                       "_haplotypeGroups_min", minSS, ".png"), 
       device = png, units="cm", 
       width=25, height=18)
ggsave(filename=paste0("07_example_hapDis_clusteringLargeHaplotypes_chr", chrID, "_Win", winSize, "_step", stepSize, 
                       "_haplotypeGroups_min", minSS, ".pdf"), 
       device = pdf, units="cm", 
       width=25, height=18)

# 







# examples with haplotype colours:
winSize="200kb"
stepSize="50kb"
minSS="50SS"
chrID="06"
MaxGap="06"
StartSet=10000000
EndSet=30000000


#Read in haplotypeNumberFile, with haplotype numbers reduced by the python script. 
hapNums.in <- read.table(paste0("clusteringLargeHaplotypes_", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_chr", chrID, 
                                ".txt"), 
                         header=T, sep="\t")
n=10
hapNums.in %>%
  filter(pos>StartSet) %>%
  filter(pos<EndSet) %>%
  filter(Mean_NAN_haplo<MaxGap) %>%
  ggplot(aes(pos/1000000, sample, fill=factor(haplotype) , colour=factor(haplotype))) +
  geom_tile()+
  scale_color_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
                                    hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
  scale_fill_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
                                   hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
  scale_x_continuous(breaks=seq(0,100,2), name = "Position (Mb)", expand = c(0,0))+
  ggtitle(paste0("Chr:", chrID))+
  theme_classic()+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"),
        axis.text = element_text(size=18),
        axis.title = element_text(size=18),
        legend.position = "none")

ggsave(filename=paste0("07_example_hapDis_clusteringLargeHaplotypes_chr", chrID, "_Win", winSize, "_step", stepSize, 
                       "_haplotypeGroups_min", minSS, ".png"), 
       device = png, units="cm", 
       width=25, height=18)
ggsave(filename=paste0("07_example_hapDis_clusteringLargeHaplotypes_chr", chrID, "_Win", winSize, "_step", stepSize, 
                       "_haplotypeGroups_min", minSS, ".pdf"), 
       device = pdf, units="cm", 
       width=25, height=18)

# 












# examples with haplotype colours:
winSize="200kb"
stepSize="50kb"
minSS="50SS"
chrID="10"
MaxGap="06"
StartSet=32000000
EndSet=48000000


#Read in haplotypeNumberFile, with haplotype numbers reduced by the python script. 
hapNums.in <- read.table(paste0("clusteringLargeHaplotypes_", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_chr", chrID, 
                                ".txt"), 
                         header=T, sep="\t")
n=10
hapNums.in %>%
  filter(pos>StartSet) %>%
  filter(pos<EndSet) %>%
  filter(Mean_NAN_haplo<MaxGap) %>%
  ggplot(aes(pos/1000000, sample, fill=factor(haplotype) , colour=factor(haplotype))) +
  geom_tile()+
  scale_color_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
                                    hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
  scale_fill_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
                                   hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
  scale_x_continuous(breaks=seq(0,100,2), name = "Position (Mb)", expand = c(0,0))+
  ggtitle(paste0("Chr:", chrID))+
  theme_classic()+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"),
        axis.text = element_text(size=18),
        axis.title = element_text(size=18),
        legend.position = "none")

ggsave(filename=paste0("07_example_hapDis_clusteringLargeHaplotypes_chr", chrID, "_Win", winSize, "_step", stepSize, 
                       "_haplotypeGroups_min", minSS, ".png"), 
       device = png, units="cm", 
       width=25, height=18)
ggsave(filename=paste0("07_example_hapDis_clusteringLargeHaplotypes_chr", chrID, "_Win", winSize, "_step", stepSize, 
                       "_haplotypeGroups_min", minSS, ".pdf"), 
       device = pdf, units="cm", 
       width=25, height=18)

# 


















# Removing short un-shared blocks:

# examples with haplotype colours:
winSize="200kb"
stepSize="50kb"
minSS="50SS"
chrID="01"
MaxGap="06"
n=10

#Read in haplotypeNumberFile, with haplotype numbers reduced by the python script. 
hapNums.in <- read.table(paste0("clusteringLargeHaplotypes_", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_chr", chrID, 
                                ".txt"), 
                         header=T, sep="\t")


StartSet=42000000
EndSet=62000000
Min_NWin=10
Min_SharedHap=3

#hapNums.in %>%
#  filter(Mean_NAN_haplo<0.8) %>%
#  group_by(sample, haplotype) %>%
#  summarise(NWin=n()) %>%
#  ggplot(aes(NWin))+
#  geom_histogram() 

MinLenHaplotypeID <- hapNums.in %>%
  filter(Mean_NAN_haplo<0.8) %>%
  group_by(sample, haplotype) %>%
  summarise(NWin=n()) %>%
  filter(NWin>Min_NWin) %>%
  ungroup() %>%
  pull(haplotype) %>%
  unique()

#hapNums.in %>%
#  filter(Mean_NAN_haplo<0.8) %>%
#  filter(haplotype %in% MinLenHaplotypeID) %>%
#  group_by(haplotype) %>%
#  summarise(NSamples=length(unique(sample))) %>%
#  ggplot(aes(NSamples))+
#  geom_histogram() 

MinSharedLenHaplotypeID <- hapNums.in %>%
  filter(Mean_NAN_haplo<0.8) %>%
  filter(haplotype %in% MinLenHaplotypeID) %>%
  group_by(haplotype) %>%
  summarise(NSamples=length(unique(sample))) %>%
  filter(NSamples>Min_SharedHap) %>%
  ungroup() %>%
  pull(haplotype) %>%
  unique()


hapNums.in %>%
  filter(pos>StartSet) %>%
  filter(pos<EndSet) %>%
  filter(Mean_NAN_haplo<MaxGap) %>%
  filter(haplotype %in% MinSharedLenHaplotypeID) %>%
  ggplot(aes(pos/1000000, sample, fill=factor(haplotype) , colour=factor(haplotype))) +
  geom_tile()+
  scale_color_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
                                    hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
  scale_fill_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
                                   hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
  scale_x_continuous(breaks=seq(0,100,2), name = "Position (Mb)", expand = c(0,0))+
  ggtitle(paste0("Chr:", chrID))+
  theme_classic()+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"),
        axis.text = element_text(size=12),
        axis.title = element_text(size=18),
        legend.position = "none")


ggsave(filename=paste0("08_example_Simplified_hapDis_clusteringLargeHaplotypes_chr", chrID,"_", StartSet/1000000, "_", EndSet/1000000, "_WinSize", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_MinNWin_", Min_NWin, "_Min_SharedHap_", Min_SharedHap, ".png"), 
       device = png, units="cm", 
       width=25, height=18)
ggsave(filename=paste0("08_example_Simplified_hapDis_clusteringLargeHaplotypes_chr", chrID,"_", StartSet/1000000, "_", EndSet/1000000, "_WinSize", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_MinNWin_", Min_NWin, "_Min_SharedHap_", Min_SharedHap, ".pdf"), 
       device = pdf, units="cm", 
       width=25, height=18)








# examples with haplotype colours:
winSize="200kb"
stepSize="50kb"
minSS="50SS"
chrID="01"
MaxGap="06"
n=10

#Read in haplotypeNumberFile, with haplotype numbers reduced by the python script. 
hapNums.in <- read.table(paste0("clusteringLargeHaplotypes_", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_chr", chrID, 
                                ".txt"), 
                         header=T, sep="\t")


StartSet=15000000
EndSet=30000000
Min_NWin=10
Min_SharedHap=3

#hapNums.in %>%
#  filter(Mean_NAN_haplo<0.8) %>%
#  group_by(sample, haplotype) %>%
#  summarise(NWin=n()) %>%
#  ggplot(aes(NWin))+
#  geom_histogram() 

MinLenHaplotypeID <- hapNums.in %>%
  filter(Mean_NAN_haplo<0.8) %>%
  group_by(sample, haplotype) %>%
  summarise(NWin=n()) %>%
  filter(NWin>Min_NWin) %>%
  ungroup() %>%
  pull(haplotype) %>%
  unique()

#hapNums.in %>%
#  filter(Mean_NAN_haplo<0.8) %>%
#  filter(haplotype %in% MinLenHaplotypeID) %>%
#  group_by(haplotype) %>%
#  summarise(NSamples=length(unique(sample))) %>%
#  ggplot(aes(NSamples))+
#  geom_histogram() 

MinSharedLenHaplotypeID <- hapNums.in %>%
  filter(Mean_NAN_haplo<0.8) %>%
  filter(haplotype %in% MinLenHaplotypeID) %>%
  group_by(haplotype) %>%
  summarise(NSamples=length(unique(sample))) %>%
  filter(NSamples>Min_SharedHap) %>%
  ungroup() %>%
  pull(haplotype) %>%
  unique()


hapNums.in %>%
  filter(pos>StartSet) %>%
  filter(pos<EndSet) %>%
  filter(Mean_NAN_haplo<MaxGap) %>%
  filter(haplotype %in% MinSharedLenHaplotypeID) %>%
  ggplot(aes(pos/1000000, sample, fill=factor(haplotype) , colour=factor(haplotype))) +
  geom_tile()+
  scale_color_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
                                    hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
  scale_fill_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
                                   hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
  scale_x_continuous(breaks=seq(0,100,2), name = "Position (Mb)", expand = c(0,0))+
  ggtitle(paste0("Chr:", chrID))+
  theme_classic()+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"),
        axis.text = element_text(size=12),
        axis.title = element_text(size=18),
        legend.position = "none")


ggsave(filename=paste0("08_example_Simplified_hapDis_clusteringLargeHaplotypes_chr", chrID,"_", StartSet/1000000, "_", EndSet/1000000, "_WinSize", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_MinNWin_", Min_NWin, "_Min_SharedHap_", Min_SharedHap, ".png"), 
       device = png, units="cm", 
       width=25, height=18)
ggsave(filename=paste0("08_example_Simplified_hapDis_clusteringLargeHaplotypes_chr", chrID,"_", StartSet/1000000, "_", EndSet/1000000, "_WinSize", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_MinNWin_", Min_NWin, "_Min_SharedHap_", Min_SharedHap, ".pdf"), 
       device = pdf, units="cm", 
       width=25, height=18)










# examples with haplotype colours:
winSize="200kb"
stepSize="50kb"
minSS="50SS"
chrID="03"
MaxGap="06"
n=10

#Read in haplotypeNumberFile, with haplotype numbers reduced by the python script. 
hapNums.in <- read.table(paste0("clusteringLargeHaplotypes_", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_chr", chrID, 
                                ".txt"), 
                         header=T, sep="\t")


StartSet=7000000
EndSet=30000000
Min_NWin=10
Min_SharedHap=3

#hapNums.in %>%
#  filter(Mean_NAN_haplo<0.8) %>%
#  group_by(sample, haplotype) %>%
#  summarise(NWin=n()) %>%
#  ggplot(aes(NWin))+
#  geom_histogram() 

MinLenHaplotypeID <- hapNums.in %>%
  filter(Mean_NAN_haplo<0.8) %>%
  group_by(sample, haplotype) %>%
  summarise(NWin=n()) %>%
  filter(NWin>Min_NWin) %>%
  ungroup() %>%
  pull(haplotype) %>%
  unique()

#hapNums.in %>%
#  filter(Mean_NAN_haplo<0.8) %>%
#  filter(haplotype %in% MinLenHaplotypeID) %>%
#  group_by(haplotype) %>%
#  summarise(NSamples=length(unique(sample))) %>%
#  ggplot(aes(NSamples))+
#  geom_histogram() 

MinSharedLenHaplotypeID <- hapNums.in %>%
  filter(Mean_NAN_haplo<0.8) %>%
  filter(haplotype %in% MinLenHaplotypeID) %>%
  group_by(haplotype) %>%
  summarise(NSamples=length(unique(sample))) %>%
  filter(NSamples>Min_SharedHap) %>%
  ungroup() %>%
  pull(haplotype) %>%
  unique()


hapNums.in %>%
  filter(pos>StartSet) %>%
  filter(pos<EndSet) %>%
  filter(Mean_NAN_haplo<MaxGap) %>%
  filter(haplotype %in% MinSharedLenHaplotypeID) %>%
  ggplot(aes(pos/1000000, sample, fill=factor(haplotype) , colour=factor(haplotype))) +
  geom_tile()+
  scale_color_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
                                    hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
  scale_fill_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
                                   hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
  scale_x_continuous(breaks=seq(0,100,2), name = "Position (Mb)", expand = c(0,0))+
  ggtitle(paste0("Chr:", chrID))+
  theme_classic()+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"),
        axis.text = element_text(size=12),
        axis.title = element_text(size=18),
        legend.position = "none")


ggsave(filename=paste0("08_example_Simplified_hapDis_clusteringLargeHaplotypes_chr", chrID,"_", StartSet/1000000, "_", EndSet/1000000, "_WinSize", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_MinNWin_", Min_NWin, "_Min_SharedHap_", Min_SharedHap, ".png"), 
       device = png, units="cm", 
       width=25, height=18)
ggsave(filename=paste0("08_example_Simplified_hapDis_clusteringLargeHaplotypes_chr", chrID,"_", StartSet/1000000, "_", EndSet/1000000, "_WinSize", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_MinNWin_", Min_NWin, "_Min_SharedHap_", Min_SharedHap, ".pdf"), 
       device = pdf, units="cm", 
       width=25, height=18)









# examples with haplotype colours:
winSize="200kb"
stepSize="50kb"
minSS="50SS"
chrID="05"
MaxGap="06"
n=10

#Read in haplotypeNumberFile, with haplotype numbers reduced by the python script. 
hapNums.in <- read.table(paste0("clusteringLargeHaplotypes_", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_chr", chrID, 
                                ".txt"), 
                         header=T, sep="\t")


StartSet=25000000
EndSet=42000000
Min_NWin=10
Min_SharedHap=3

#hapNums.in %>%
#  filter(Mean_NAN_haplo<0.8) %>%
#  group_by(sample, haplotype) %>%
#  summarise(NWin=n()) %>%
#  ggplot(aes(NWin))+
#  geom_histogram() 

MinLenHaplotypeID <- hapNums.in %>%
  filter(Mean_NAN_haplo<0.8) %>%
  group_by(sample, haplotype) %>%
  summarise(NWin=n()) %>%
  filter(NWin>Min_NWin) %>%
  ungroup() %>%
  pull(haplotype) %>%
  unique()

#hapNums.in %>%
#  filter(Mean_NAN_haplo<0.8) %>%
#  filter(haplotype %in% MinLenHaplotypeID) %>%
#  group_by(haplotype) %>%
#  summarise(NSamples=length(unique(sample))) %>%
#  ggplot(aes(NSamples))+
#  geom_histogram() 

MinSharedLenHaplotypeID <- hapNums.in %>%
  filter(Mean_NAN_haplo<0.8) %>%
  filter(haplotype %in% MinLenHaplotypeID) %>%
  group_by(haplotype) %>%
  summarise(NSamples=length(unique(sample))) %>%
  filter(NSamples>Min_SharedHap) %>%
  ungroup() %>%
  pull(haplotype) %>%
  unique()


hapNums.in %>%
  filter(pos>StartSet) %>%
  filter(pos<EndSet) %>%
  filter(Mean_NAN_haplo<MaxGap) %>%
  filter(haplotype %in% MinSharedLenHaplotypeID) %>%
  ggplot(aes(pos/1000000, sample, fill=factor(haplotype) , colour=factor(haplotype))) +
  geom_tile()+
  scale_color_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
                                    hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
  scale_fill_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
                                   hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
  scale_x_continuous(breaks=seq(0,100,2), name = "Position (Mb)", expand = c(0,0))+
  ggtitle(paste0("Chr:", chrID))+
  theme_classic()+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"),
        axis.text = element_text(size=12),
        axis.title = element_text(size=18),
        legend.position = "none")


ggsave(filename=paste0("08_example_Simplified_hapDis_clusteringLargeHaplotypes_chr", chrID,"_", StartSet/1000000, "_", EndSet/1000000, "_WinSize", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_MinNWin_", Min_NWin, "_Min_SharedHap_", Min_SharedHap, ".png"), 
       device = png, units="cm", 
       width=25, height=18)
ggsave(filename=paste0("08_example_Simplified_hapDis_clusteringLargeHaplotypes_chr", chrID,"_", StartSet/1000000, "_", EndSet/1000000, "_WinSize", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_MinNWin_", Min_NWin, "_Min_SharedHap_", Min_SharedHap, ".pdf"), 
       device = pdf, units="cm", 
       width=25, height=18)















# examples with haplotype colours:
winSize="200kb"
stepSize="50kb"
minSS="50SS"
chrID="06"
MaxGap="06"
n=10

#Read in haplotypeNumberFile, with haplotype numbers reduced by the python script. 
hapNums.in <- read.table(paste0("clusteringLargeHaplotypes_", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_chr", chrID, 
                                ".txt"), 
                         header=T, sep="\t")


StartSet=10000000
EndSet=36000000
Min_NWin=5
Min_SharedHap=3

#hapNums.in %>%
#  filter(Mean_NAN_haplo<0.8) %>%
#  group_by(sample, haplotype) %>%
#  summarise(NWin=n()) %>%
#  ggplot(aes(NWin))+
#  geom_histogram() 

MinLenHaplotypeID <- hapNums.in %>%
  filter(Mean_NAN_haplo<0.8) %>%
  group_by(sample, haplotype) %>%
  summarise(NWin=n()) %>%
  filter(NWin>Min_NWin) %>%
  ungroup() %>%
  pull(haplotype) %>%
  unique()

#hapNums.in %>%
#  filter(Mean_NAN_haplo<0.8) %>%
#  filter(haplotype %in% MinLenHaplotypeID) %>%
#  group_by(haplotype) %>%
#  summarise(NSamples=length(unique(sample))) %>%
#  ggplot(aes(NSamples))+
#  geom_histogram() 

MinSharedLenHaplotypeID <- hapNums.in %>%
  filter(Mean_NAN_haplo<0.8) %>%
  filter(haplotype %in% MinLenHaplotypeID) %>%
  group_by(haplotype) %>%
  summarise(NSamples=length(unique(sample))) %>%
  filter(NSamples>Min_SharedHap) %>%
  ungroup() %>%
  pull(haplotype) %>%
  unique()


hapNums.in %>%
  filter(pos>StartSet) %>%
  filter(pos<EndSet) %>%
  filter(Mean_NAN_haplo<MaxGap) %>%
  filter(haplotype %in% MinSharedLenHaplotypeID) %>%
  ggplot(aes(pos/1000000, sample, fill=factor(haplotype) , colour=factor(haplotype))) +
  geom_tile()+
  scale_color_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
                                    hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
  scale_fill_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
                                   hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
  scale_x_continuous(breaks=seq(0,100,2), name = "Position (Mb)", expand = c(0,0))+
  ggtitle(paste0("Chr:", chrID))+
  theme_classic()+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"),
        axis.text = element_text(size=12),
        axis.title = element_text(size=18),
        legend.position = "none")


ggsave(filename=paste0("08_example_Simplified_hapDis_clusteringLargeHaplotypes_chr", chrID,"_", StartSet/1000000, "_", EndSet/1000000, "_WinSize", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_MinNWin_", Min_NWin, "_Min_SharedHap_", Min_SharedHap, ".png"), 
       device = png, units="cm", 
       width=25, height=18)
ggsave(filename=paste0("08_example_Simplified_hapDis_clusteringLargeHaplotypes_chr", chrID,"_", StartSet/1000000, "_", EndSet/1000000, "_WinSize", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_MinNWin_", Min_NWin, "_Min_SharedHap_", Min_SharedHap, ".pdf"), 
       device = pdf, units="cm", 
       width=25, height=18)












# examples with haplotype colours:
winSize="200kb"
stepSize="50kb"
minSS="50SS"
chrID="07"
MaxGap="06"
n=10

#Read in haplotypeNumberFile, with haplotype numbers reduced by the python script. 
hapNums.in <- read.table(paste0("clusteringLargeHaplotypes_", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_chr", chrID, 
                                ".txt"), 
                         header=T, sep="\t")


StartSet=10000000
EndSet=32000000
Min_NWin=5
Min_SharedHap=3

#hapNums.in %>%
#  filter(Mean_NAN_haplo<0.8) %>%
#  group_by(sample, haplotype) %>%
#  summarise(NWin=n()) %>%
#  ggplot(aes(NWin))+
#  geom_histogram() 

MinLenHaplotypeID <- hapNums.in %>%
  filter(Mean_NAN_haplo<0.8) %>%
  group_by(sample, haplotype) %>%
  summarise(NWin=n()) %>%
  filter(NWin>Min_NWin) %>%
  ungroup() %>%
  pull(haplotype) %>%
  unique()

#hapNums.in %>%
#  filter(Mean_NAN_haplo<0.8) %>%
#  filter(haplotype %in% MinLenHaplotypeID) %>%
#  group_by(haplotype) %>%
#  summarise(NSamples=length(unique(sample))) %>%
#  ggplot(aes(NSamples))+
#  geom_histogram() 

MinSharedLenHaplotypeID <- hapNums.in %>%
  filter(Mean_NAN_haplo<0.8) %>%
  filter(haplotype %in% MinLenHaplotypeID) %>%
  group_by(haplotype) %>%
  summarise(NSamples=length(unique(sample))) %>%
  filter(NSamples>Min_SharedHap) %>%
  ungroup() %>%
  pull(haplotype) %>%
  unique()


hapNums.in %>%
  filter(pos>StartSet) %>%
  filter(pos<EndSet) %>%
  filter(Mean_NAN_haplo<MaxGap) %>%
  filter(haplotype %in% MinSharedLenHaplotypeID) %>%
  ggplot(aes(pos/1000000, sample, fill=factor(haplotype) , colour=factor(haplotype))) +
  geom_tile()+
  scale_color_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
                                    hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
  scale_fill_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
                                   hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
  scale_x_continuous(breaks=seq(0,100,2), name = "Position (Mb)", expand = c(0,0))+
  ggtitle(paste0("Chr:", chrID))+
  theme_classic()+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"),
        axis.text = element_text(size=12),
        axis.title = element_text(size=18),
        legend.position = "none")


ggsave(filename=paste0("08_example_Simplified_hapDis_clusteringLargeHaplotypes_chr", chrID,"_", StartSet/1000000, "_", EndSet/1000000, "_WinSize", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_MinNWin_", Min_NWin, "_Min_SharedHap_", Min_SharedHap, ".png"), 
       device = png, units="cm", 
       width=25, height=18)
ggsave(filename=paste0("08_example_Simplified_hapDis_clusteringLargeHaplotypes_chr", chrID,"_", StartSet/1000000, "_", EndSet/1000000, "_WinSize", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_MinNWin_", Min_NWin, "_Min_SharedHap_", Min_SharedHap, ".pdf"), 
       device = pdf, units="cm", 
       width=25, height=18)










# examples with haplotype colours:
winSize="200kb"
stepSize="50kb"
minSS="50SS"
chrID="08"
MaxGap="06"
n=10

#Read in haplotypeNumberFile, with haplotype numbers reduced by the python script. 
hapNums.in <- read.table(paste0("clusteringLargeHaplotypes_", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_chr", chrID, 
                                ".txt"), 
                         header=T, sep="\t")


StartSet=10000000
EndSet=40000000
Min_NWin=5
Min_SharedHap=3

#hapNums.in %>%
#  filter(Mean_NAN_haplo<0.8) %>%
#  group_by(sample, haplotype) %>%
#  summarise(NWin=n()) %>%
#  ggplot(aes(NWin))+
#  geom_histogram() 

MinLenHaplotypeID <- hapNums.in %>%
  filter(Mean_NAN_haplo<0.8) %>%
  group_by(sample, haplotype) %>%
  summarise(NWin=n()) %>%
  filter(NWin>Min_NWin) %>%
  ungroup() %>%
  pull(haplotype) %>%
  unique()

#hapNums.in %>%
#  filter(Mean_NAN_haplo<0.8) %>%
#  filter(haplotype %in% MinLenHaplotypeID) %>%
#  group_by(haplotype) %>%
#  summarise(NSamples=length(unique(sample))) %>%
#  ggplot(aes(NSamples))+
#  geom_histogram() 

MinSharedLenHaplotypeID <- hapNums.in %>%
  filter(Mean_NAN_haplo<0.8) %>%
  filter(haplotype %in% MinLenHaplotypeID) %>%
  group_by(haplotype) %>%
  summarise(NSamples=length(unique(sample))) %>%
  filter(NSamples>Min_SharedHap) %>%
  ungroup() %>%
  pull(haplotype) %>%
  unique()


hapNums.in %>%
  filter(pos>StartSet) %>%
  filter(pos<EndSet) %>%
  filter(Mean_NAN_haplo<MaxGap) %>%
  filter(haplotype %in% MinSharedLenHaplotypeID) %>%
  ggplot(aes(pos/1000000, sample, fill=factor(haplotype) , colour=factor(haplotype))) +
  geom_tile()+
  scale_color_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
                                    hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
  scale_fill_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
                                   hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
  scale_x_continuous(breaks=seq(0,100,2), name = "Position (Mb)", expand = c(0,0))+
  ggtitle(paste0("Chr:", chrID))+
  theme_classic()+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"),
        axis.text = element_text(size=12),
        axis.title = element_text(size=18),
        legend.position = "none")


ggsave(filename=paste0("08_example_Simplified_hapDis_clusteringLargeHaplotypes_chr", chrID,"_", StartSet/1000000, "_", EndSet/1000000, "_WinSize", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_MinNWin_", Min_NWin, "_Min_SharedHap_", Min_SharedHap, ".png"), 
       device = png, units="cm", 
       width=25, height=18)
ggsave(filename=paste0("08_example_Simplified_hapDis_clusteringLargeHaplotypes_chr", chrID,"_", StartSet/1000000, "_", EndSet/1000000, "_WinSize", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_MinNWin_", Min_NWin, "_Min_SharedHap_", Min_SharedHap, ".pdf"), 
       device = pdf, units="cm", 
       width=25, height=18)














# examples with haplotype colours:
winSize="200kb"
stepSize="50kb"
minSS="50SS"
chrID="09"
MaxGap="06"
n=8

#Read in haplotypeNumberFile, with haplotype numbers reduced by the python script. 
hapNums.in <- read.table(paste0("clusteringLargeHaplotypes_", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_chr", chrID, 
                                ".txt"), 
                         header=T, sep="\t")


StartSet=14000000
EndSet=50000000
Min_NWin=5
Min_SharedHap=3

#hapNums.in %>%
#  filter(Mean_NAN_haplo<0.8) %>%
#  group_by(sample, haplotype) %>%
#  summarise(NWin=n()) %>%
#  ggplot(aes(NWin))+
#  geom_histogram() 

MinLenHaplotypeID <- hapNums.in %>%
  filter(Mean_NAN_haplo<0.8) %>%
  group_by(sample, haplotype) %>%
  summarise(NWin=n()) %>%
  filter(NWin>Min_NWin) %>%
  ungroup() %>%
  pull(haplotype) %>%
  unique()

#hapNums.in %>%
#  filter(Mean_NAN_haplo<0.8) %>%
#  filter(haplotype %in% MinLenHaplotypeID) %>%
#  group_by(haplotype) %>%
#  summarise(NSamples=length(unique(sample))) %>%
#  ggplot(aes(NSamples))+
#  geom_histogram() 

MinSharedLenHaplotypeID <- hapNums.in %>%
  filter(Mean_NAN_haplo<0.8) %>%
  filter(haplotype %in% MinLenHaplotypeID) %>%
  group_by(haplotype) %>%
  summarise(NSamples=length(unique(sample))) %>%
  filter(NSamples>Min_SharedHap) %>%
  ungroup() %>%
  pull(haplotype) %>%
  unique()


hapNums.in %>%
  filter(pos>StartSet) %>%
  filter(pos<EndSet) %>%
  filter(Mean_NAN_haplo<MaxGap) %>%
  filter(haplotype %in% MinSharedLenHaplotypeID) %>%
  ggplot(aes(pos/1000000, sample, fill=factor(haplotype) , colour=factor(haplotype))) +
  geom_tile()+
  scale_color_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
                                    hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
  scale_fill_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
                                   hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
  scale_x_continuous(breaks=seq(0,100,2), name = "Position (Mb)", expand = c(0,0))+
  ggtitle(paste0("Chr:", chrID))+
  theme_classic()+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"),
        axis.text = element_text(size=12),
        axis.title = element_text(size=18),
        legend.position = "none")


ggsave(filename=paste0("08_example_Simplified_hapDis_clusteringLargeHaplotypes_chr", chrID,"_", StartSet/1000000, "_", EndSet/1000000, "_WinSize", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_MinNWin_", Min_NWin, "_Min_SharedHap_", Min_SharedHap, ".png"), 
       device = png, units="cm", 
       width=25, height=18)
ggsave(filename=paste0("08_example_Simplified_hapDis_clusteringLargeHaplotypes_chr", chrID,"_", StartSet/1000000, "_", EndSet/1000000, "_WinSize", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_MinNWin_", Min_NWin, "_Min_SharedHap_", Min_SharedHap, ".pdf"), 
       device = pdf, units="cm", 
       width=25, height=18)
















# examples with haplotype colours:
winSize="200kb"
stepSize="50kb"
minSS="50SS"
chrID="10"
MaxGap="06"
n=8

#Read in haplotypeNumberFile, with haplotype numbers reduced by the python script. 
hapNums.in <- read.table(paste0("clusteringLargeHaplotypes_", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_chr", chrID, 
                                ".txt"), 
                         header=T, sep="\t")


StartSet=00000000
EndSet=70000000
Min_NWin=5
Min_SharedHap=3

#hapNums.in %>%
#  filter(Mean_NAN_haplo<0.8) %>%
#  group_by(sample, haplotype) %>%
#  summarise(NWin=n()) %>%
#  ggplot(aes(NWin))+
#  geom_histogram() 

MinLenHaplotypeID <- hapNums.in %>%
  filter(Mean_NAN_haplo<0.8) %>%
  group_by(sample, haplotype) %>%
  summarise(NWin=n()) %>%
  filter(NWin>Min_NWin) %>%
  ungroup() %>%
  pull(haplotype) %>%
  unique()

#hapNums.in %>%
#  filter(Mean_NAN_haplo<0.8) %>%
#  filter(haplotype %in% MinLenHaplotypeID) %>%
#  group_by(haplotype) %>%
#  summarise(NSamples=length(unique(sample))) %>%
#  ggplot(aes(NSamples))+
#  geom_histogram() 

MinSharedLenHaplotypeID <- hapNums.in %>%
  filter(Mean_NAN_haplo<0.8) %>%
  filter(haplotype %in% MinLenHaplotypeID) %>%
  group_by(haplotype) %>%
  summarise(NSamples=length(unique(sample))) %>%
  filter(NSamples>Min_SharedHap) %>%
  ungroup() %>%
  pull(haplotype) %>%
  unique()


hapNums.in %>%
  filter(pos>StartSet) %>%
  filter(pos<EndSet) %>%
  filter(Mean_NAN_haplo<MaxGap) %>%
  filter(haplotype %in% MinSharedLenHaplotypeID) %>%
  ggplot(aes(pos/1000000, sample, fill=factor(haplotype) , colour=factor(haplotype))) +
  geom_tile()+
  scale_color_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
                                    hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
  scale_fill_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
                                   hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
  scale_x_continuous(breaks=seq(0,100,2), name = "Position (Mb)", expand = c(0,0))+
  ggtitle(paste0("Chr:", chrID))+
  theme_classic()+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"),
        axis.text = element_text(size=12),
        axis.title = element_text(size=18),
        legend.position = "none")


ggsave(filename=paste0("08_example_Simplified_hapDis_clusteringLargeHaplotypes_chr", chrID,"_", StartSet/1000000, "_", EndSet/1000000, "_WinSize", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_MinNWin_", Min_NWin, "_Min_SharedHap_", Min_SharedHap, ".png"), 
       device = png, units="cm", 
       width=25, height=18)
ggsave(filename=paste0("08_example_Simplified_hapDis_clusteringLargeHaplotypes_chr", chrID,"_", StartSet/1000000, "_", EndSet/1000000, "_WinSize", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_MinNWin_", Min_NWin, "_Min_SharedHap_", Min_SharedHap, ".pdf"), 
       device = pdf, units="cm", 
       width=25, height=18)












# examples with haplotype colours:
winSize="200kb"
stepSize="50kb"
minSS="50SS"
chrID="11"
MaxGap="06"
n=8

#Read in haplotypeNumberFile, with haplotype numbers reduced by the python script. 
hapNums.in <- read.table(paste0("clusteringLargeHaplotypes_", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_chr", chrID, 
                                ".txt"), 
                         header=T, sep="\t")


StartSet=10000000
EndSet=47000000
Min_NWin=5
Min_SharedHap=3

#hapNums.in %>%
#  filter(Mean_NAN_haplo<0.8) %>%
#  group_by(sample, haplotype) %>%
#  summarise(NWin=n()) %>%
#  ggplot(aes(NWin))+
#  geom_histogram() 

MinLenHaplotypeID <- hapNums.in %>%
  filter(Mean_NAN_haplo<0.8) %>%
  group_by(sample, haplotype) %>%
  summarise(NWin=n()) %>%
  filter(NWin>Min_NWin) %>%
  ungroup() %>%
  pull(haplotype) %>%
  unique()

#hapNums.in %>%
#  filter(Mean_NAN_haplo<0.8) %>%
#  filter(haplotype %in% MinLenHaplotypeID) %>%
#  group_by(haplotype) %>%
#  summarise(NSamples=length(unique(sample))) %>%
#  ggplot(aes(NSamples))+
#  geom_histogram() 

MinSharedLenHaplotypeID <- hapNums.in %>%
  filter(Mean_NAN_haplo<0.8) %>%
  filter(haplotype %in% MinLenHaplotypeID) %>%
  group_by(haplotype) %>%
  summarise(NSamples=length(unique(sample))) %>%
  filter(NSamples>Min_SharedHap) %>%
  ungroup() %>%
  pull(haplotype) %>%
  unique()


hapNums.in %>%
  filter(pos>StartSet) %>%
  filter(pos<EndSet) %>%
  filter(Mean_NAN_haplo<MaxGap) %>%
  filter(haplotype %in% MinSharedLenHaplotypeID) %>%
  ggplot(aes(pos/1000000, sample, fill=factor(haplotype) , colour=factor(haplotype))) +
  geom_tile()+
  scale_color_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
                                    hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
  scale_fill_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
                                   hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
  scale_x_continuous(breaks=seq(0,100,2), name = "Position (Mb)", expand = c(0,0))+
  ggtitle(paste0("Chr:", chrID))+
  theme_classic()+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"),
        axis.text = element_text(size=12),
        axis.title = element_text(size=18),
        legend.position = "none")


ggsave(filename=paste0("08_example_Simplified_hapDis_clusteringLargeHaplotypes_chr", chrID,"_", StartSet/1000000, "_", EndSet/1000000, "_WinSize", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_MinNWin_", Min_NWin, "_Min_SharedHap_", Min_SharedHap, ".png"), 
       device = png, units="cm", 
       width=25, height=18)
ggsave(filename=paste0("08_example_Simplified_hapDis_clusteringLargeHaplotypes_chr", chrID,"_", StartSet/1000000, "_", EndSet/1000000, "_WinSize", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_MinNWin_", Min_NWin, "_Min_SharedHap_", Min_SharedHap, ".pdf"), 
       device = pdf, units="cm", 
       width=25, height=18)













# examples with haplotype colours:
winSize="200kb"
stepSize="50kb"
minSS="50SS"
chrID="12"
MaxGap="06"
n=8

#Read in haplotypeNumberFile, with haplotype numbers reduced by the python script. 
hapNums.in <- read.table(paste0("clusteringLargeHaplotypes_", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_chr", chrID, 
                                ".txt"), 
                         header=T, sep="\t")


StartSet=10000000
EndSet=50000000
Min_NWin=5
Min_SharedHap=3

#hapNums.in %>%
#  filter(Mean_NAN_haplo<0.8) %>%
#  group_by(sample, haplotype) %>%
#  summarise(NWin=n()) %>%
#  ggplot(aes(NWin))+
#  geom_histogram() 

MinLenHaplotypeID <- hapNums.in %>%
  filter(Mean_NAN_haplo<0.8) %>%
  group_by(sample, haplotype) %>%
  summarise(NWin=n()) %>%
  filter(NWin>Min_NWin) %>%
  ungroup() %>%
  pull(haplotype) %>%
  unique()

#hapNums.in %>%
#  filter(Mean_NAN_haplo<0.8) %>%
#  filter(haplotype %in% MinLenHaplotypeID) %>%
#  group_by(haplotype) %>%
#  summarise(NSamples=length(unique(sample))) %>%
#  ggplot(aes(NSamples))+
#  geom_histogram() 

MinSharedLenHaplotypeID <- hapNums.in %>%
  filter(Mean_NAN_haplo<0.8) %>%
  filter(haplotype %in% MinLenHaplotypeID) %>%
  group_by(haplotype) %>%
  summarise(NSamples=length(unique(sample))) %>%
  filter(NSamples>Min_SharedHap) %>%
  ungroup() %>%
  pull(haplotype) %>%
  unique()


hapNums.in %>%
  filter(pos>StartSet) %>%
  filter(pos<EndSet) %>%
  filter(Mean_NAN_haplo<MaxGap) %>%
  filter(haplotype %in% MinSharedLenHaplotypeID) %>%
  ggplot(aes(pos/1000000, sample, fill=factor(haplotype) , colour=factor(haplotype))) +
  geom_tile()+
  scale_color_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
                                    hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
  scale_fill_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
                                   hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
  scale_x_continuous(breaks=seq(0,100,2), name = "Position (Mb)", expand = c(0,0))+
  ggtitle(paste0("Chr:", chrID))+
  theme_classic()+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey80"),
        axis.text = element_text(size=12),
        axis.title = element_text(size=18),
        legend.position = "none")


ggsave(filename=paste0("08_example_Simplified_hapDis_clusteringLargeHaplotypes_chr", chrID,"_", StartSet/1000000, "_", EndSet/1000000, "_WinSize", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_MinNWin_", Min_NWin, "_Min_SharedHap_", Min_SharedHap, ".png"), 
       device = png, units="cm", 
       width=25, height=18)
ggsave(filename=paste0("08_example_Simplified_hapDis_clusteringLargeHaplotypes_chr", chrID,"_", StartSet/1000000, "_", EndSet/1000000, "_WinSize", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_MinNWin_", Min_NWin, "_Min_SharedHap_", Min_SharedHap, ".pdf"), 
       device = pdf, units="cm", 
       width=25, height=18)














# WHole chromosomes with haplotype colours:
winSize="200kb"
stepSize="50kb"
minSS="50SS"
# chrID="12"
MaxGap="06"
n=10

chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")


for (chrID in chr_list){
  
  
  #Read in haplotypeNumberFile, with haplotype numbers reduced by the python script. 
  hapNums.in <- read.table(paste0("clusteringLargeHaplotypes_", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_chr", chrID, 
                                  ".txt"), 
                           header=T, sep="\t")
  
  Min_NWin=5
  Min_SharedHap=3
  
  #hapNums.in %>%
  #  filter(Mean_NAN_haplo<0.8) %>%
  #  group_by(sample, haplotype) %>%
  #  summarise(NWin=n()) %>%
  #  ggplot(aes(NWin))+
  #  geom_histogram() 
  
  MinLenHaplotypeID <- hapNums.in %>%
    filter(Mean_NAN_haplo<0.8) %>%
    group_by(sample, haplotype) %>%
    summarise(NWin=n()) %>%
    filter(NWin>Min_NWin) %>%
    ungroup() %>%
    pull(haplotype) %>%
    unique()
  
  #hapNums.in %>%
  #  filter(Mean_NAN_haplo<0.8) %>%
  #  filter(haplotype %in% MinLenHaplotypeID) %>%
  #  group_by(haplotype) %>%
  #  summarise(NSamples=length(unique(sample))) %>%
  #  ggplot(aes(NSamples))+
  #  geom_histogram() 
  
  MinSharedLenHaplotypeID <- hapNums.in %>%
    filter(Mean_NAN_haplo<0.8) %>%
    filter(haplotype %in% MinLenHaplotypeID) %>%
    group_by(haplotype) %>%
    summarise(NSamples=length(unique(sample))) %>%
    filter(NSamples>Min_SharedHap) %>%
    ungroup() %>%
    pull(haplotype) %>%
    unique()
  
  
  hapNums.in %>%
    filter(Mean_NAN_haplo<MaxGap) %>%
    filter(haplotype %in% MinSharedLenHaplotypeID) %>%
    ggplot(aes(pos/1000000, sample, fill=factor(haplotype) , colour=factor(haplotype))) +
    geom_tile()+
    scale_color_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
                                      hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
    scale_fill_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
                                     hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
    scale_x_continuous(breaks=seq(0,100,2), name = "Position (Mb)", expand = c(0,0))+
    ggtitle(paste0("Chr:", chrID))+
    theme_classic()+
    theme_classic()+
    theme(panel.grid.major = element_line(color = "grey80"),
          panel.grid.minor = element_line(color = "grey80"),
          axis.text = element_text(size=12),
          axis.title = element_text(size=18),
          legend.position = "none")
  
  if (chrID=="01"){
    ggsave(filename=paste0("09_wholeChr_Simplified_hapDis_clusteringLargeHaplotypes_chr", chrID, "_WinSize", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_MinNWin_", Min_NWin, "_Min_SharedHap_", Min_SharedHap, ".png"), 
           device = png, units="cm", 
           width=30, height=18)
    ggsave(filename=paste0("09_wholeChr_Simplified_hapDis_clusteringLargeHaplotypes_chr", chrID, "_WinSize", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_MinNWin_", Min_NWin, "_Min_SharedHap_", Min_SharedHap, ".pdf"), 
           device = pdf, units="cm", 
           width=30, height=18)
  } else {
    ggsave(filename=paste0("09_wholeChr_Simplified_hapDis_clusteringLargeHaplotypes_chr", chrID, "_WinSize", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_MinNWin_", Min_NWin, "_Min_SharedHap_", Min_SharedHap, ".png"), 
           device = png, units="cm", 
           width=25, height=18)
    ggsave(filename=paste0("09_wholeChr_Simplified_hapDis_clusteringLargeHaplotypes_chr", chrID, "_WinSize", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_MinNWin_", Min_NWin, "_Min_SharedHap_", Min_SharedHap, ".pdf"), 
           device = pdf, units="cm", 
           width=25, height=18)
  }
  
}















# Whole chromosomes with haplotype colours, but clustering haplotypes:
winSize="200kb"
stepSize="50kb"
minSS="50SS"
# chrID="12"
MaxGap="06"
n=10

singleL <- c("A","B","C","D","E","F","G","H","I","J")
AcName <- c("WhR", "BdF", "EgH", "PrW", "Flo", "Ack", "Fla", "Kat", "Lum", "EdB")

4*12
singleL_list <- c()
for (i in singleL){
  singleL_list <- c(singleL_list, rep(i, 48))
}

AcName_list <- c()
for (i in AcName){
  AcName_list <- c(AcName_list, rep(i, 48))
}

label_table <- data.frame(hapN=seq(1,48), singleL_list, AcName_list, simpleHapNum=rep(seq(1,4), 120)) %>%
  mutate(sample=paste0(singleL_list, "_hap", hapN)) %>%
  mutate(label=paste0(AcName_list, "_", simpleHapNum)) %>%
  select(sample, label)


chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

chrID="01"

for (chrID in chr_list){
  
  #Read in haplotypeNumberFile, with haplotype numbers reduced by the python script. 
  hapNums.in <- read.table(paste0("clusteringLargeHaplotypes_", winSize, "_step", stepSize,
                                  "_haplotypeGroups_min", minSS, "_chr", chrID, ".txt"), 
                           header=T, sep="\t")
  
  hapNums.in2 <- hapNums.in %>%
    left_join(label_table, by="sample") 
  
  matrix_hapNums.in <- hapNums.in2 %>%
    filter(Mean_NAN_haplo<0.8) %>%
    select(pos, label, haplotype) %>%
    spread(key=pos, value=haplotype, fill = NA) %>%
    select(-label) %>%
    as.matrix()
  
  samples_ID <- hapNums.in2 %>%
    filter(Mean_NAN_haplo<0.8) %>%
    select(pos, label, haplotype) %>%
    spread(key=pos, value=haplotype, fill = NA) %>%
    select(label) %>% unlist() %>% as.vector()
  
  d <- dist(matrix_hapNums.in, method = "euclidean")
  H.fit <- hclust(d, method="ward.D")
  H.fit$labels<-samples_ID
  
  hapNums.in2$label <- factor(hapNums.in2$label, levels = H.fit$labels[H.fit$order])
  
  Min_NWin=5
  Min_SharedHap=3
  
  MinLenHaplotypeID <- hapNums.in2 %>%
    filter(Mean_NAN_haplo<0.8) %>%
    group_by(sample, haplotype) %>%
    summarise(NWin=n()) %>%
    filter(NWin>Min_NWin) %>%
    ungroup() %>%
    pull(haplotype) %>%
    unique()
  
  MinSharedLenHaplotypeID <- hapNums.in2 %>%
    filter(Mean_NAN_haplo<0.8) %>%
    filter(haplotype %in% MinLenHaplotypeID) %>%
    group_by(haplotype) %>%
    summarise(NSamples=length(unique(sample))) %>%
    filter(NSamples>Min_SharedHap) %>%
    ungroup() %>%
    pull(haplotype) %>%
    unique()
  
  
  hapNums.in2 %>%
    filter(Mean_NAN_haplo<MaxGap) %>%
    filter(haplotype %in% MinSharedLenHaplotypeID) %>%
    ggplot(aes(pos/1000000, label, fill=factor(haplotype) , colour=factor(haplotype))) +
    geom_tile()+
    scale_color_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
                                      hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
    scale_fill_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
                                     hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
    scale_x_continuous(breaks=seq(0,100,5), name = "Position (Mb)", expand = c(0,0), limits = c(0,90))+
    scale_y_discrete(name="Haplotype")+
    ggtitle(paste0("Chr:", chrID))+
    theme_classic()+
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=18),
          axis.line.x = element_blank(),
          legend.position = "none")
  
  ggsave(filename=paste0("10_wholeChr_Simplified_hapDis_clusteringLargeHaplotypes_chr", chrID, "_WinSize", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_MinNWin_", Min_NWin, "_Min_SharedHap_", Min_SharedHap, ".png"), 
         device = png, units="cm", 
         width=30, height=18)
  ggsave(filename=paste0("10_wholeChr_Simplified_hapDis_clusteringLargeHaplotypes_chr", chrID, "_WinSize", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_MinNWin_", Min_NWin, "_Min_SharedHap_", Min_SharedHap, ".pdf"), 
         device = pdf, units="cm", 
         width=30, height=18)
  
  # if (chrID=="01"){
  #   ggsave(filename=paste0("10_wholeChr_Simplified_hapDis_clusteringLargeHaplotypes_chr", chrID, "_WinSize", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_MinNWin_", Min_NWin, "_Min_SharedHap_", Min_SharedHap, ".png"), 
  #          device = png, units="cm", 
  #          width=30, height=18)
  #   ggsave(filename=paste0("10_wholeChr_Simplified_hapDis_clusteringLargeHaplotypes_chr", chrID, "_WinSize", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_MinNWin_", Min_NWin, "_Min_SharedHap_", Min_SharedHap, ".pdf"), 
  #          device = pdf, units="cm", 
  #          width=30, height=18)
  # } else {
  #   ggsave(filename=paste0("10_wholeChr_Simplified_hapDis_clusteringLargeHaplotypes_chr", chrID, "_WinSize", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_MinNWin_", Min_NWin, "_Min_SharedHap_", Min_SharedHap, ".png"), 
  #          device = png, units="cm", 
  #          width=30, height=18)
  #   ggsave(filename=paste0("10_wholeChr_Simplified_hapDis_clusteringLargeHaplotypes_chr", chrID, "_WinSize", winSize, "_step", stepSize, "_haplotypeGroups_min", minSS, "_MinNWin_", Min_NWin, "_Min_SharedHap_", Min_SharedHap, ".pdf"), 
  #          device = pdf, units="cm", 
  #          width=30, height=18)
  # }
}


```





# Genotype calling 40 haplotypes + 20 wild samples:



## Reference genomes per haplotype:



script: reference_bd.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00
#SBATCH -J reference_bd

sample=$1

samtools faidx whole_$sample.fa

# Index file - bwa 
bwa index whole_$sample.fa

```



Job:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/02_reference/Haplotype_resolved_potato_assemblies_v2_20240126/reference_perSample

ls -1 ../*fasta | awk -F_ '{print $1}' | sort | uniq | sed 's@../@@g' | grep -v "O" > list_samples.txt
sample=A

for sample in $(cat list_samples.txt); do echo $sample
cat ../$sample""*fasta | sed 's/\([ACTG]\)>/\1\n>/g' > whole_$sample.fa
done

# Creating the FASTA sequence dictionary file:

mamba deactivate
for sample in $(cat list_samples.txt); do echo $sample
/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/gatk-4.2.0.0/gatk CreateSequenceDictionary -R whole_$sample.fa
done


mamba activate env_others
for sample in $(cat list_samples.txt); do echo $sample
sbatch reference_bd.sh $sample
done
# 29.04.2024: 2509702..2509711

```





## Mapping reads to each assembled genomes to split reads per haplotype:



Script: bwa.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=4-00:00:00
#SBATCH -J BWA_Mapping

Threads=20

Reference=$1
sample=$2
Read1=$3

ReadGroup="@RG\tID:wHAIPI016410-33\tPL:PacBio\tPU:150211_I191_FCC6L7EANXX_L5_wHAIPI016410-33\tLB:wHAIPI016410-33\tSM:"$sample"\tCN:BGI"

echo "the files are:"
echo $Read1
echo $Reference

bwa mem -x pacbio -t $Threads -R $ReadGroup -M $Reference $Read1 | samtools view -@ 20 -bS - | samtools sort -@ 20 -o $sample.bam -

```

Job:

```bash

cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/01_mapping_toHap

mamba activate env_others

cp /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/02_reference/Haplotype_resolved_potato_assemblies_v2_20240126/reference_perSample/list_samples.txt .

for sample in $(cat list_samples.txt); do echo $sample
reference=/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/02_reference/Haplotype_resolved_potato_assemblies_v2_20240126/reference_perSample/whole_$sample.fa
sbatch ./bwa.sh $reference $sample /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/01_raw_data/Founders_10_samples_long_reads/sample_$sample"_"zcat_HiFi_trim.fastq.gz 
done
# 30.04.2024: 2510701..2510710

# for i in {2510701..2510710}; do stop $i ; done

```



### Mapping with minimap2

script: minimap2.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00
#SBATCH -J Minimap2

Reference=$1
sample=$2
Read1=$3

minimap2 -ax map-hifi -t 20 -N 1 --secondary=no $Reference $Read1 | samtools view -@ 20 -bS - | samtools sort -@ 20 -o $sample.bam -

samtools index -@ 20 $sample.bam

```

Job:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/01_mapping_toHap_minimap

mamba activate env_others

cp /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/02_reference/Haplotype_resolved_potato_assemblies_v2_20240126/reference_perSample/list_samples.txt .

for sample in $(cat list_samples.txt); do echo $sample
reference=/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/02_reference/Haplotype_resolved_potato_assemblies_v2_20240126/reference_perSample/whole_$sample.fa
sbatch ./minimap2.sh $reference $sample /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/01_raw_data/Founders_10_samples_long_reads/sample_$sample"_"zcat_HiFi_trim.fastq.gz 
done
# 30.04.2024: 2510689..2510698

# for i in {2509988..2509997}; do stop $i ; done
```







## Split reads by haplotype



Script: extractReads.sh 

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00
#SBATCH -J ExtractReads


bam_file=$1
list_chr_haplotype=$2
haplotype_name=$3

#samtools view -@ 4 -bL $list_chr_haplotype -o $haplotype_name.bam $bam_file
#samtools sort -@ 4 -n -o $haplotype_name.sort.bam $haplotype_name.bam
#samtools fastq -@ 4 $haplotype_name.sort.bam > $haplotype_name.fastq


samtools view -@ 5 -bL $list_chr_haplotype $bam_file | samtools sort -@ 5 - | samtools fastq -@ 5 - > $haplotype_name.fastq

gzip $haplotype_name.fastq


```

Job:

```bash
mamba activate env_others

cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/02_splitReadperHap

cp /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/02_reference/Haplotype_resolved_potato_assemblies_v2_20240126/reference_perSample/list_samples.txt .

for sample in $(cat list_samples.txt); do echo $sample
for hap in $(seq -w 1 4)
do echo $hap
grep ">" /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/02_reference/Haplotype_resolved_potato_assemblies_v2_20240126/$sample"_"hap$hap"_"genome.fasta | grep -v PGA | sed 's/>//g' > list_hap_$sample.$hap.txt
reference=/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/02_reference/Haplotype_resolved_potato_assemblies_v2_20240126/reference_perSample/whole_$sample.fa
grep -wf list_hap_$sample.$hap.txt $reference.fai | awk '{print $1"\t0\t"$2}' > list_$sample.hap$hap.bed
bam_file=/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/01_mapping_toHap_minimap/$sample.bam
sbatch extractReads.sh $bam_file list_$sample.hap$hap.bed readsSample$sample"_"Hap_$hap
done
done
# 01.05.2024: 2514208..2514247

for i in {2514167..2514206}; do stop $i ; done
```



## Mapping to reference DM



Script: bwa.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00
#SBATCH -J BWA_Mapping

Threads=20

Reference=$1
sample=$2
Read1=$3

ReadGroup="@RG\tID:wHAIPI016410-33\tPL:PacBio\tPU:150211_I191_FCC6L7EANXX_L5_wHAIPI016410-33\tLB:wHAIPI016410-33\tSM:"$sample"\tCN:BGI"

echo "the files are:"
echo $Read1
echo $Reference

bwa mem -x pacbio -t $Threads -R $ReadGroup -M $Reference $Read1 | samtools view -@ 20 -bS - | samtools sort -@ 20 -o $sample.bam -


```

Job:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/01_population_data/15_VariantsPerHaplotype/03_mapping_ReferenceDM

mamba activate env_others

cp /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/02_reference/Haplotype_resolved_potato_assemblies_v2_20240126/reference_perSample/list_samples.txt .

reference=/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/02_reference/DM_1-3_516_R44_potato_genome_assembly.v6.1/DM_1-3_516_R44_potato_genome_assembly.v6.1.sm.fa

# Cultivar haplotypes:
for sample in $(cat list_samples.txt); do echo $sample
for hap in $(seq -w 1 4)
do echo $hap
sbatch ./bwa.sh $reference Sample$sample"_"Hap_$hap /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/02_splitReadperHap/readsSample$sample"_"Hap_$hap.fastq.gz  
done
done
# 30.04.2024: 2516787..2516826



# Wild Samples:

ls -1 /dss/dsslegfs01/pn29fi/pn29fi-dss-0004/Potato_wildSamples_Hifi/ | sed 's/.fastq.gz//g' > list_wildSamples.txt

reference=/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/02_reference/DM_1-3_516_R44_potato_genome_assembly.v6.1/DM_1-3_516_R44_potato_genome_assembly.v6.1.sm.fa

for sample in $(cat list_wildSamples.txt); do echo $sample
sbatch ./bwa.sh $reference $sample /dss/dsslegfs01/pn29fi/pn29fi-dss-0004/Potato_wildSamples_Hifi/$sample.fastq.gz  
done
# 01.05.2024: 2515125..2515144

# for i in {2510333..2510342}; do stop $i ; done

```



### Mappint to DM using minimap2



Script: minimap2.sh

```bash

#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00
#SBATCH -J Minimap2

Reference=$1
sample=$2
Read1=$3

minimap2 -ax map-hifi -t 20 -R @RG\\tID:A00253_251_HTN2JDSXY.2\\tPL:PACBIO\tLB:LB1\\tSM:$sample -N 1 $Reference $Read1 | samtools view -@ 20 -bS - | samtools sort -@ 20 -o $sample.bam -

samtools index -@ 20 $sample.bam

```

Job: 

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/01_population_data/15_VariantsPerHaplotype/03_mapping_ReferenceDM

mamba activate env_others

cp /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/02_reference/Haplotype_resolved_potato_assemblies_v2_20240126/reference_perSample/list_samples.txt .

reference=/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/02_reference/DM_1-3_516_R44_potato_genome_assembly.v6.1/DM_1-3_516_R44_potato_genome_assembly.v6.1.sm.fa

# Cultivar haplotypes:
for sample in $(cat list_samples.txt); do echo $sample
for hap in $(seq -w 1 4)
do echo $hap
sbatch ./minimap2.sh $reference Sample$sample"_"Hap_$hap /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/02_splitReadperHap/readsSample$sample"_"Hap_$hap.fastq.gz  
done
done
# 01.05.2024: 2516916..2516955



# Wild Samples:

ls -1 /dss/dsslegfs01/pn29fi/pn29fi-dss-0004/Potato_wildSamples_Hifi/ | sed 's/.fastq.gz//g' > list_wildSamples.txt

reference=/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/02_reference/DM_1-3_516_R44_potato_genome_assembly.v6.1/DM_1-3_516_R44_potato_genome_assembly.v6.1.sm.fa

for sample in $(cat list_wildSamples.txt); do echo $sample
sbatch ./minimap2.sh $reference $sample /dss/dsslegfs01/pn29fi/pn29fi-dss-0004/Potato_wildSamples_Hifi/$sample.fastq.gz  
done
# 01.05.2024: 2516957..2516976


./minimap2.sh $reference test /dss/dsslegfs01/pn29fi/pn29fi-dss-0004/Potato_wildSamples_Hifi/$sample.fastq.gz  

# for i in {2515125..2515144}; do stop $i ; done

```



## HaplotypeCaller



### DeepVariant



script: deepvariant.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=15000mb
#SBATCH --time=20:00:00
#SBATCH -J deepvariant


module use /dss/dsslegfs02/pn73se/pn73se-dss-0000/spack/modules/$LRZ_INSTRSET/linux* ; module load user_spack
module load charliecloud

sample=$1
run=$2
region=$3
#regions.bed

mkdir -p ./vcf
mkdir -p ./vcf/$run


ch-run --bind $PWD:/bedfiles/ --bind $PWD/vcf/$run:/output/ --bind /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/02_reference/DM_1-3_516_R44_potato_genome_assembly.v6.1:/reference_file --bind /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/03_mapping_ReferenceDM:/bam_files -w /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/deepvariant/deepvariant_1_4_0_charliecloud_container -- /opt/deepvariant/bin/run_deepvariant \
    --model_type PACBIO \
    --ref /reference_file/DM_1-3_516_R44_potato_genome_assembly.v6.1.sm.fa \
    --reads /bam_files/$sample.bam \
    --output_vcf /output/$sample.$run.vcf.gz \
    --output_gvcf /output/$sample.$run.gvcf.gz \
    --regions /bedfiles/$region \
    --num_shards 8



```

job: 



```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling

module use /dss/dsslegfs02/pn73se/pn73se-dss-0000/spack/modules/$LRZ_INSTRSET/linux* ; module load user_spack
module load charliecloud

grep "^@SQ" /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/02_reference/DM_1-3_516_R44_potato_genome_assembly.v6.1/interval_files/0000-scattered.interval_list | grep chr | sed 's/SN://g' | sed 's/LN://g' | awk '{print $2"\t1\t"$3}' > chr_size.bed

# Cultivars

ls -1 /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/03_mapping_ReferenceDM/*bam | grep -v tmp | grep -v test | grep "/Sample" | sed 's@/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/03_mapping_ReferenceDM/@@g' | sed 's@.bam@@g' > list_samples_culivarHap.txt

for run in $(seq -w 3 3 12)
do echo $run
for sample in $( cat list_samples_culivarHap.txt ); do echo $sample
start=$(echo $run | sed 's/^0//g'); echo $start
head -n $start chr_size.bed | tail -n 3 > list_region.$run.bed
sbatch deepvariant.sh $sample $run list_region.$run.bed
done
done
# 06.05.2024: 2533487,2533488, 2534595..2534752


for run in $(seq -w 3 3 12)
do echo $run
for sample in $( cat list_samples_culivarHap.txt ); do echo $sample
start=$(echo $run | sed 's/^0//g'); echo $start
head -n $start chr_size.bed | tail -n 3 > list_region.$run.bed
if [ "$run" == "03" ]
then
if [ "$sample" == "SampleA_Hap_1" ]
then
echo "done"
else
if [ "$sample" == "SampleA_Hap_2" ]
then
echo "done"
else
sbatch deepvariant.sh $sample $run list_region.$run.bed
fi
fi
else
sbatch deepvariant.sh $sample $run list_region.$run.bed
fi
done
done


# wild samples:

ls -1 /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/03_mapping_ReferenceDM/*bam | grep -v tmp | grep -v test | grep "/SRR" | sed 's@/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/03_mapping_ReferenceDM/@@g' | sed 's@.bam@@g' > list_samples_wild.txt

for run in $(seq -w 3 3 12)
do echo $run
for sample in $( cat list_samples_wild.txt ); do echo $sample
start=$(echo $run | sed 's/^0//g'); echo $start
head -n $start chr_size.bed | tail -n 3 > list_region.$run.bed
sbatch deepvariant2.sh $sample $run list_region.$run.bed
done
done
# 06.05.2024: 4046644..4046723


#cm2
for sample in $( cat list_samples_wild.txt | tail -n 12); do echo $sample
for run in $(seq -w 3 3 12); do echo $run
start=$(echo $run | sed 's/^0//g'); echo $start
head -n $start chr_size.bed | tail -n 3 > list_region.$run.bed
sbatch deepvariant_cm2.sh $sample $run list_region.$run.bed
done
done
# 07.05.2024: 647473..647520


#cm2_std
for sample in $( cat list_samples_wild.txt | head -8 ); do echo $sample
for run in $(seq -w 3 3 12); do echo $run
start=$(echo $run | sed 's/^0//g'); echo $start
head -n $start chr_size.bed | tail -n 3 > list_region.$run.bed
sbatch deepvariant_cm2_st.sh $sample $run list_region.$run.bed
done
done
# 07.05.2024: 162519..162550


#mpp3 - didn't work
for sample in $( cat list_samples_wild.txt | head -8); do echo $sample
for run in $(seq -w 3 3 12); do echo $run
start=$(echo $run | sed 's/^0//g'); echo $start
head -n $start chr_size.bed | tail -n 3 > list_region.$run.bed
sbatch deepvariant_mpp3.sh $sample $run list_region.$run.bed
done
done
# 07.05.2024: 661212..661243

for i in {162407..162454}; do  rm slurm-$i.out ; done

#for i in {4045609..4045688}; do stopserial $i ; done
# for i in {162374..162406}; do  scancel --cluster=cm2 $i ; done


grep TIME slurm-* | awk -F: '{print $1}' > non_finishedJobs.txt

grep "" slurm-* | grep "out:time seq" | grep -f non_finishedJobs.txt | sed 's/^.*bam_files\///g' | sed 's/.bam.*bedfiles\/list_region./@/g'  | sed 's/.bed.*//g' > ToRun_non_finishedJobs.txt

for line in $( cat ToRun_non_finishedJobs.txt ); do echo $line
sample=$(echo $line | awk -F@ '{print $1}'); echo $sample
run=$(echo $line | awk -F@ '{print $2}'); echo $run
sbatch deepvariant2.sh $sample $run list_region.$run.bed
done
# 09.05.2024: 4047387..4047411


```

### Merge VCF

script: glnexus.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=12
#SBATCH --mem=40000mb
#SBATCH --time=20:00:00
#SBATCH -J glnexus


module use /dss/dsslegfs02/pn73se/pn73se-dss-0000/spack/modules/$LRZ_INSTRSET/linux* ; module load user_spack
module load charliecloud

export OMP_NUM_THREADS=12

run=$1

mkdir -p ./merged_vcf
mkdir -p ./merged_vcf/$run
mkdir -p ./temporal_glnexus
mkdir -p ./temporal_glnexus/tem_$run

export TMPDIR=$PWD"/"temporal_glnexus/tem_$run


ch-run --bind $PWD/temporal_glnexus/tem_$run:/temporal --bind $PWD/vcf/$run:/gvcf/ --bind $PWD/merged_vcf/$run:/output_mergedvcf/ -w /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_raul/glnexus/glnexus_image \
    -- bash -c "/usr/local/bin/glnexus_cli \
    --config DeepVariant \
    --mem-gbytes 40 \
    --threads 12 \
    --dir /temporal/$run \
    /gvcf/*.g.vcf.gz | bcftools view -Oz - > /output_mergedvcf/ChrGroup_$run.vcf.gz && \
    tabix -p vcf /output_mergedvcf/ChrGroup_$run.vcf.gz"



```

Job:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling

module use /dss/dsslegfs02/pn73se/pn73se-dss-0000/spack/modules/$LRZ_INSTRSET/linux* ; module load user_spack
module load charliecloud


ls -1 vcf/[0-1]*/*gvcf.g* > listFiles.txt

for file in $( cat listFiles.txt ); do echo $file
output=$(echo $file | sed 's/gvcf/g\.vcf/g')
echo $output
mv $file $output
done

for run in $(seq -w 3 3 12)
do echo $run
sbatch glnexus2.sh $run 
done
# 11.05.2024: 2538145..2538148

# 2538100 (run=03)
sbatch glnexus2.sh $run 

```





### Filter low quality variants

script: QUA_val.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --cpus-per-task=1
#SBATCH --mem=4700mb
#SBATCH --time=2:00:00
#SBATCH -J QUAL_Val

#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=4700mb
#SBATCH --time=2:00:00
#SBATCH -J QUAL_Val

#run=03
run=$1
output=$2

echo "chr,QUAL,NVar" > $output
zcat merged_vcf/$run"/"ChrGroup_$run.vcf.gz | grep -v "^#" | awk '{print $1"\t"$6}' | sort | uniq -c | awk '{print $2","$3","$1}' >> $output


```

Job: 

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling

for run in $(seq -w 3 3 12)
do echo $run
sbatch ./QUA_val.sh $run merged_vcf/QUAL_val.$run.csv
done
# 16.05.2024: 4054340..4054343


for run in $(seq -w 3 3 12)
do echo $run
sbatch ./QUA_val2.sh $run merged_vcf/QUALval.$run.csv
done
# 16.05.2024: 2543059..2543062

for i in {2543356..2543367}; do stop $i ; done
for i in {4054336..4054339}; do stopserial $i ; done

```

Plot:

```R

# .libPaths("/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/SITG_packages/R")
.libPaths()

setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling")

library(ggplot2)
library(tidyverse)

table_qual_03 <- read_csv("./merged_vcf/QUALval.03.csv", T)
table_qual_06 <- read_csv("./merged_vcf/QUALval.06.csv", T)
table_qual_09 <- read_csv("./merged_vcf/QUALval.09.csv", T)
table_qual_12 <- read_csv("./merged_vcf/QUALval.12.csv", T)

table_qual <- rbind(table_qual_03, table_qual_06) %>%
  rbind(table_qual_09) %>%
  rbind(table_qual_12) 



# Allele Count:

table_qual %>%
  # mutate(intervalQ=floor(QUAL/5)*5) %>%
  mutate(intervalQ=ifelse(QUAL>30, "Included", "Excluded")) %>%
  group_by(chr, intervalQ) %>%
  summarise(NVar_int=sum(NVar)) %>%
  ggplot(aes(chr, NVar_int))+
  geom_bar(stat="identity")+
  scale_y_continuous(breaks=seq(0,10000000,1000000)) +
  #labs(x = "Allele Fq.",
  #     y = "Num. Variant Sites (e6)") +
  facet_grid(intervalQ~., scale="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor.x = element_line(color = "grey80"))

# ggsave("plots/01_FilteringVar_QUAL_min30.png", height=5, width=8)

```



Script: filterVar_QUAL30.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=4700mb
#SBATCH --time=4:00:00
#SBATCH -J QUAL_Val

#run=03
run=$1
chr=$2

mkdir -p merged_vcf/chr

bcftools view -O z --targets $chr -o merged_vcf/chr"/"$chr"_"filterQUAL30.vcf.gz -e 'QUAL<=30' merged_vcf/$run"/"ChrGroup_$run.vcf.gz

tabix -p vcf merged_vcf/chr"/"$chr"_"filterQUAL30.vcf.gz


```

Job:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling

mamba activate env_others

for run in $(seq -w 3 3 12)
do echo $run
grep -v "QUAL" merged_vcf/QUALval.$run.csv | awk -F, '{print $1}' | sort | uniq > list_Chr.$run.txt
for chr in $( cat list_Chr.$run.txt ); do echo $chr
sbatch filterVar_QUAL30.sh $run $chr
done
done
# 16.05.2024: 2543157..2543168

```



## Split VCF per Windows



script: splitvcf.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=8:00:00
#SBATCH -J splitvcf

vcf=$1
output=$2
winSize=$3
nSnps=$4

mkdir -p splitVCF
mkdir -p splitVCF/$output

zcat $vcf | grep -v "##contig=<ID=scaffold" | awk -F"\t" 'BEGIN{ FS=OFS="\t" }{if($1 ~ /^chr/) $3="." ; print $0}' > splitVCF/tem_$output.vcf

java -jar /dss/dsshome1/lxc03/di36guz2/miniconda3/envs/myjvarkit/share/jvarkit-2024.04.20-0/jvarkit.jar biostar497922 -n $nSnps -o splitVCF/$output -D $winSize ./splitVCF/tem_$output.vcf

rm ./splitVCF/tem_$output.vcf

```



job: 



```bash
/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling

mamba activate myjvarkit

for chr in $(seq -w 1 12)
do echo $chr
vcf=/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/merged_vcf/chr/chr$chr"_"filterQUAL30.vcf.gz
sbatch splitvcf.sh $vcf $chr 100k 100000 
done
# 20.05.2024: 2544615..2544626



```



## Phylogeny



### **Phylogeny without cultivars**



Script: filterSamples.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00

chr=$1

mkdir -p splitVCF
mkdir -p splitVCF/$chr

ls -1 /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/splitVCF/$chr"/"*.vcf.gz > list_vcf_$chr.txt


for vcf in $(cat list_vcf_$chr.txt); do echo $win
vcf_file=${vcf##*/}
vcf_file_base=${vcf_file%.vcf.gz}
vcftools --gzvcf $vcf --recode --recode-INFO-all --keep list_WildSamples.txt  --out splitVCF/$chr"/"$vcf_file_base
gzip -c splitVCF/$chr"/"$vcf_file_base.recode.vcf > splitVCF/$chr"/"$vcf_file_base.vcf.gz
rm splitVCF/$chr"/"$vcf_file_base.recode.vcf
done



```



Job:

```bash

cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/06_phylogeny_WildSamples

mamba activate env_others

vcf=/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/splitVCF/01/split.000852.vcf.gz

bcftools query -l $vcf | grep -v Sample > list_WildSamples.txt


for chr in $(seq -w 1 12)
do echo $chr
sbatch filterSamples.sh $chr 
done
# 20.05.2024: 2544647..2544658

for i in {2544635..2544646}; do stop $i ; rm slurm-$i.out ; done

```



script: vcf2phylip.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00
#SBATCH -J vcf2phylip


chr=$1

mkdir -p splitVCF_phy
mkdir -p splitVCF_phy/$chr

ls -1 /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/splitVCF/$chr"/"*.vcf.gz > list_vcf_$chr.txt

for vcf in $(cat list_vcf_$chr.txt); do echo $win
vcf_file=${vcf##*/}
vcf_file_base=${vcf_file%.vcf.gz}
python /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/vcf2phylip/vcf2phylip.py --input ./splitVCF/$chr"/"$vcf_file --fasta --nexus --nexus-binary --min-samples-locus 10 --output-folder ./splitVCF_phy/$chr 
done


"done finished"

```



Job:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/06_phylogeny_WildSamples

mamba activate mypython3

for i in $(seq -w 1 12); do echo $i
sbatch ./vcf2phylip.sh $i
done
# 20.05.2024: 2544660..2544671


```



script: IQTreeperWin.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=10-00:00:00
#SBATCH -J IQTree


chr=$1

mkdir -p IQTree
mkdir -p IQTree/$chr

ls -1 /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/06_phylogeny_WildSamples/splitVCF_phy/$chr"/"*.phy > list_phy_$chr.txt

for phy in $(cat list_phy_$chr.txt); do echo $phy
phy_file=${phy##*/}
phy_file_base=${phy_file%.phy}
iqtree -s $phy -m GTR -B 1000 -T 4 --prefix ./IQTree/$chr"/"$phy_file_base
done

"done finished"


```



Job:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/06_phylogeny_WildSamples

mamba activate env_others

for i in $(seq -w 1 12); do echo $i
sbatch ./IQTreeperWin.sh $i
done
# 20.05.2024: 2544675..2544686


```



script: combine_trees_outgroup.py

```python
#!/usr/bin/env python

"""                                                                                 
python merge_seqs.py list_samples.txt 'alig_*'
"""                                                                                 

from ete3 import Tree
from Bio import SeqIO
from Bio.Seq import Seq

import glob
import sys

# list of trees
listTrees = [line.strip() for line in open(sys.argv[1])]

newtrees = []
for line in listTrees:
	# line_tree=open(line, 'r')
	tree_line = Tree(line)
	# set outgroup:
	tree_line.set_outgroup(str(sys.argv[3]))
	newtrees += [tree_line.write()]
	# print(newtrees)

# Print trees in a new file
rootedtrees = open(sys.argv[2], 'w')
for tree in newtrees:
	rootedtrees.write(tree + "\n")


```

Job:

```bash
# conda create -n ete3 -c bioconda
mamba activate ete3
#mamba install -c etetoolkit ete3 ete_toolchain
#mamba install -c bioconda biopython

cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/06_phylogeny_WildSamples

#chr=01
for i in $(seq -w 1 12); do echo $i
ls -1 /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/06_phylogeny_WildSamples/IQTree/$i"/"*treefile > list_trees.chr$i.txt
python combine_trees_outgroup.py list_trees.chr$i.txt merged_trees.chr$i.treefile SRR15458990
nw_ed merged_trees.chr$i.treefile 'i & b<=25' o > merged_trees.chr$i.min25.treefile
done

```





script: astral.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00
#SBATCH -J Astral


genetrees=$1.treefile
spptree=astral.$1.treefile
output1=astral.$1.scored.treefile
output2=astral.$1.poly.treefile
log1=astral.$1.log
log2=astral.$1.scored.log
log3=astral.$1.poly.log


# Astral
java -Xmx15G -D"java.library.path=/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/SITG_packages/miniconda3/envs/ete3/share/astral-tree-5.7.8-1/lib/" -jar /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/SITG_packages/miniconda3/envs/ete3/share/astral-tree-5.7.8-1/astral.5.7.8.jar -i $genetrees -o $spptree 2> $log1

# scoring_astral:
java -Xmx15G -D"java.library.path=/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/SITG_packages/miniconda3/envs/ete3/share/astral-tree-5.7.8-1/lib/" -jar /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/SITG_packages/miniconda3/envs/ete3/share/astral-tree-5.7.8-1/astral.5.7.8.jar -q $spptree -i $genetrees -o $output1 -t 8 2> $log2

# polytomy_test_astral:
java -Xmx15G -D"java.library.path=/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/SITG_packages/miniconda3/envs/ete3/share/astral-tree-5.7.8-1/lib/" -jar /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/SITG_packages/miniconda3/envs/ete3/share/astral-tree-5.7.8-1/astral.5.7.8.jar -q $spptree -i $genetrees -o $output2 -t 10 2> $log3


```



Job:



```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/06_phylogeny_WildSamples

mamba activate ete3

for i in $(seq -w 1 12); do echo $i
sbatch astral.sh merged_trees.chr$i.min25
done
# 20.05.2024: 2544732..2544743

# for i in {2545278..2545289}; do stop $i ; rm slurm-$i* ; done

```



### **Phylogeny with cultivars**



script: vcf2phylip.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00
#SBATCH -J vcf2phylip


chr=$1

mkdir -p splitVCF_phy
mkdir -p splitVCF_phy/$chr

ls -1 /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/splitVCF/$chr"/"*.vcf.gz > list_vcf_$chr.txt

for vcf in $(cat list_vcf_$chr.txt); do echo $win
python /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/vcf2phylip/vcf2phylip.py --input $vcf --fasta --nexus --nexus-binary --min-samples-locus 10 --output-folder ./splitVCF_phy/$chr 
done


"done finished"

```



per chr script: vcf2phylip_chr.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00
#SBATCH -J vcf2phylip

#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000mb
#SBATCH --time=12:00:00
#SBATCH -J Example_script

vcf=$1

mkdir -p chr_phy

python /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/vcf2phylip/vcf2phylip.py --input $vcf --fasta --nexus --nexus-binary --min-samples-locus 10 --output-folder ./chr_phy

"done finished"

```



Job:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/07_phylogeny_withCulWildSamples

mamba activate mypython3

for i in $(seq -w 1 12); do echo $i
sbatch ./vcf2phylip.sh $i
done
# 21.05.2024: 2545265..2545276


# Per Chromosome:
for i in $(seq -w 1 12); do echo $i
sbatch ./vcf2phylip_chr.sh /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/merged_vcf/chr/chr$i"_"filterQUAL30.vcf.gz
done
# 21.05.2024: 2545278..2545289

for i in $(seq -w 1 12); do echo $i
sbatch ./vcf2phylip_chr2.sh /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/merged_vcf/chr/chr$i"_"filterQUAL30.vcf.gz
done
# 21.05.2024: 4059739..4059750

```



script: IQTreeperWin.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=4-00:00:00
#SBATCH -J IQTree


chr=$1

mkdir -p IQTree
mkdir -p IQTree/$chr

ls -1 /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/07_phylogeny_withCulWildSamples/splitVCF_phy/$chr"/"*.phy > list_phy_$chr.txt

for phy in $(cat list_phy_$chr.txt); do echo $phy
phy_file=${phy##*/}
phy_file_base=${phy_file%.phy}
iqtree -s $phy -m GTR -B 1000 -T 4 --prefix ./IQTree/$chr"/"$phy_file_base
done

"done finished"


```

script: IQTreeperChr.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=10-00:00:00
#SBATCH -J IQTree_Chr

#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00
#SBATCH -J IQTree_Chr

phy=$1

mkdir -p IQTree_Chr

phy_file=${phy##*/}
phy_file_base=${phy_file%.phy}
iqtree -s $phy -m GTR -B 1000 -T 4 --prefix ./IQTree_Chr/$phy_file_base

echo "done finished"


```



Job:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/07_phylogeny_withCulWildSamples

mamba activate env_others

for i in $(seq -w 1 12); do echo $i
sbatch ./IQTreeperWin.sh $i
done
# 21.05.2024: 4060628..4060639

# for i in {2545389..2545400}; do stop $i ; rm slurm-$i* ; done

# per chromosome:
for i in $(seq -w 1 12); do echo $i
sbatch ./IQTreeperChr.sh /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/07_phylogeny_withCulWildSamples/chr_phy_serial/chr$i"_"filterQUAL30.min10.phy
done
# 21.05.2024: 4060450..4060461 (01 is the last one)


```



script: combine_trees_outgroup.py

```python
#!/usr/bin/env python

"""                                                                                 
python merge_seqs.py list_samples.txt 'alig_*'
"""                                                                                 

from ete3 import Tree
from Bio import SeqIO
from Bio.Seq import Seq

import glob
import sys

# list of trees
listTrees = [line.strip() for line in open(sys.argv[1])]

newtrees = []
for line in listTrees:
	# line_tree=open(line, 'r')
	tree_line = Tree(line)
	# set outgroup:
	tree_line.set_outgroup(str(sys.argv[3]))
	newtrees += [tree_line.write()]
	# print(newtrees)

# Print trees in a new file
rootedtrees = open(sys.argv[2], 'w')
for tree in newtrees:
	rootedtrees.write(tree + "\n")


```

Job:

```bash
# conda create -n ete3 -c bioconda
mamba activate ete3
#mamba install -c etetoolkit ete3 ete_toolchain
#mamba install -c bioconda biopython

cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/07_phylogeny_withCulWildSamples

# Using as outgroup: SRR15458990_S_morelliforme_C1_2
for i in $(seq -w 1 12); do echo $i
ls -1 /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/07_phylogeny_withCulWildSamples/IQTree/$i"/"*treefile > list_trees.chr$i.txt
python combine_trees_outgroup.py list_trees.chr$i.txt merged_trees.chr$i.treefile SRR15458990
nw_ed merged_trees.chr$i.treefile 'i & b<=25' o > merged_trees.chr$i.min25.treefile
done

#whole Genome:
# Using as outgroup: SRR15458990_S_morelliforme_C1_2
ls -1 /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/07_phylogeny_withCulWildSamples/IQTree/*"/"*treefile > list_trees.WGenome.txt
python combine_trees_outgroup.py list_trees.WGenome.txt merged_trees.WGenome.treefile SRR15458990
nw_ed merged_trees.WGenome.treefile 'i & b<=25' o > merged_trees.WGenome.min25.treefile


```





script: astral.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00
#SBATCH -J Astral


genetrees=$1.treefile
spptree=astral.$1.treefile
output1=astral.$1.scored.treefile
output2=astral.$1.poly.treefile
log1=astral.$1.log
log2=astral.$1.scored.log
log3=astral.$1.poly.log


# Astral
java -Xmx15G -D"java.library.path=/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/SITG_packages/miniconda3/envs/ete3/share/astral-tree-5.7.8-1/lib/" -jar /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/SITG_packages/miniconda3/envs/ete3/share/astral-tree-5.7.8-1/astral.5.7.8.jar -i $genetrees -o $spptree 2> $log1

# scoring_astral:
java -Xmx15G -D"java.library.path=/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/SITG_packages/miniconda3/envs/ete3/share/astral-tree-5.7.8-1/lib/" -jar /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/SITG_packages/miniconda3/envs/ete3/share/astral-tree-5.7.8-1/astral.5.7.8.jar -q $spptree -i $genetrees -o $output1 -t 8 2> $log2

# polytomy_test_astral:
java -Xmx15G -D"java.library.path=/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/SITG_packages/miniconda3/envs/ete3/share/astral-tree-5.7.8-1/lib/" -jar /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/SITG_packages/miniconda3/envs/ete3/share/astral-tree-5.7.8-1/astral.5.7.8.jar -q $spptree -i $genetrees -o $output2 -t 10 2> $log3


```



Job:



```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/07_phylogeny_withCulWildSamples

mamba activate ete3

sbatch astral.sh merged_trees.WGenome.min25
# 2546172

for i in $(seq -w 1 12); do echo $i
sbatch astral.sh merged_trees.chr$i.min25
done
# 23.05.2024: 2546160..2546171

# for i in {2545278..2545289}; do stop $i ; rm slurm-$i* ; done

```



### GATK - Not Used



*Script:* haplotypecaller.sh

```bash
#!/bin/bash -l
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=48:00:00
#SBATCH -J HaplotypeCaller


#!/bin/bash -l
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=48:00:00
#SBATCH -J HaplotypeCaller

# in this script you need the location of the bam file, the interval file produced in the previous step, and the vcf file folder for the output

bamdir=/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/03_mapping_ReferenceDM
intdir=/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/02_reference/DM_1-3_516_R44_potato_genome_assembly.v6.1/interval_files

genomedir=/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/02_reference/DM_1-3_516_R44_potato_genome_assembly.v6.1
genome=${genomedir}/DM_1-3_516_R44_potato_genome_assembly.v6.1.sm.fa

#Sample and interval(0000-0029) as positional argument
SM=$1
RUN=$2
PLOIDY=$3

echo $SM
echo $RUN
echo $PLOIDY
echo $genome

# Used Only native pairhmm
#/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/gatk-4.2.0.0/gatk --java-options "-Xmx8g" HaplotypeCaller \
java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -Xmx9g -jar /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/gatk-4.2.0.0/gatk-package-4.2.0.0-local.jar HaplotypeCaller \
-R ${genome} -I ${bamdir}/${SM}.bam -pairHMM AVX_LOGLESS_CACHING \
-ploidy ${PLOIDY} \
-L ${intdir}/${RUN}-scattered.interval_list \
-O ${PWD}/vcf/scatter/${RUN}/${SM}.${RUN}.gvcf.gz -ERC GVCF  \
--output-mode EMIT_ALL_CONFIDENT_SITES --min-base-quality-score 0


```



**Job**

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/01_population_data/06_genotypeCalling

mkdir vcf
mkdir ./vcf/scatter/

ls -1 /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/02_reference/DM_1-3_516_R44_potato_genome_assembly.v6.1/interval_files/* | sed 's@/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/02_reference/DM_1-3_516_R44_potato_genome_assembly.v6.1/interval_files/@@g' | sed 's@-scattered.interval_list@@g' > list_intervals.txt

ls -1 /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/03_mapping_ReferenceDM/*bam | grep -v tmp | grep -v test | grep "/Sample" | sed 's@/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/03_mapping_ReferenceDM/@@g' | sed 's@.bam@@g' > list_samples_culivarHap.txt

moduleL
module load openjdk

# serial
# Cultivar haplotypes:
for sample in $( cat list_samples_culivarHap.txt | grep B); do echo $sample
for run in $( cat list_intervals.txt ); do echo $run
sbatch haplotypecaller.sh $sample $run 1
done
done
# 02.05.2024: 4039434..4039553 (A) 4039554..4039673 (B)
# 4040406, 4040407 for sample=SampleA_Hap_1, run=0024, run=0025





#biohpc
# Cultivar haplotypes:
for sample in $( cat list_samples_culivarHap.txt | grep D); do echo $sample
for run in $( cat list_intervals.txt ); do echo $run
sbatch haplotypecaller2.sh $sample $run 1
done
done
# 02.05.2024: 2518737..2518856 (C) 2518857..2518976 (D)

for sample in $( cat list_samples_culivarHap2.txt ); do echo $sample
for run in $( cat list_intervals.txt ); do echo $run
sbatch haplotypecaller2.sh $sample $run 1
done
done
# 02.05.2024: 2519856..2520545



#cm2
# Cultivar haplotypes:
for sample in $( cat list_samples_culivarHap.txt | grep E); do echo $sample
for run in $( cat list_intervals.txt ); do echo $run
sbatch haplotypecaller3.sh $sample $run 1
done
done
# 02.05.2024: 643780..643829 (E partly)


for i in {2517783..2518529}; do stop $i ; done





# re running jobs: 

ls -l vcf/scatter/*gz | awk '{if ($5<1000000) print $9}' | sed 's/.gvcf.gz//' | sed 's@vcf/scatter/@@' | sed 's/\./@/g' > unfinished_jobs.txt


ls -1 vcf/scatter/*tbi | sed 's/.gvcf.gz.tbi//' | sed 's@vcf/scatter/@@' | sed 's/\./@/g' > finishedJobs.txt
ls -1 vcf/scatter/*gz | sed 's/.gvcf.gz//' | sed 's@vcf/scatter/@@' | sed 's/\./@/g' | grep -vf finishedJobs.txt > unfinished_jobs.txt


jobinfo | grep Haplotyp | awk '{print $1}' > runningJobs.txt
jobinfoserial | grep Haplotyp | awk '{print $1}' >> runningJobs.txt
jobinfocm2 | grep Haplotyp | awk '{print $1}' >> runningJobs.txt

grep -w -B 1 "^00[0-9][0-9]" slurm-* | grep -v "\-\-" | sed 's@slurm-.*.out:@@g' | sed ':a;N;/\nslurm/!s/\n/ /;ta;P;D' | sed 's/slurm-//g' | sed 's/.out-/ /g' | sed 's/ /@/g' > list_runJobs.txt

grep -f unfinished_jobs.txt list_runJobs.txt | grep -vf runningJobs.txt > jobs_toRepeat.txt

for line in $( cat jobs_toRepeat.txt | head -250); do echo $line
jobid=$(echo $line | awk -F@ '{print $1}') ; echo $jobid
sample=$(echo $line | awk -F@ '{print $2}') ; echo $sample
run=$(echo $line | awk -F@ '{print $3}') ; echo $run
mv slurm-$jobid.out file_slurm/
mv vcf/scatter/$sample.$run.gvcf.gz temporal/
sbatch haplotypecaller2.sh $sample $run 1
done

# 2521388..2521637

for line in $( cat jobs_toRepeat.txt | head -450 | tail -n 200); do echo $line
jobid=$(echo $line | awk -F@ '{print $1}') ; echo $jobid
sample=$(echo $line | awk -F@ '{print $2}') ; echo $sample
run=$(echo $line | awk -F@ '{print $3}') ; echo $run
mv slurm-$jobid.out file_slurm/
mv vcf/scatter/$sample.$run.gvcf.gz temporal/
sbatch haplotypecaller.sh $sample $run 1
done
# 4041834..4042033

for line in $( cat jobs_toRepeat.txt | tail -560 ); do echo $line
jobid=$(echo $line | awk -F@ '{print $1}') ; echo $jobid
sample=$(echo $line | awk -F@ '{print $2}') ; echo $sample
run=$(echo $line | awk -F@ '{print $3}') ; echo $run
mv slurm-$jobid.out file_slurm/
mv vcf/scatter/$sample.$run.gvcf.gz temporal/
sbatch haplotypecaller2.sh $sample $run 1
done
# 2523002..2523561



for i in {2530503..2531702}; do stop $i ; done

for run in $( cat list_intervals.txt ); do echo $run
mkdir vcf/scatter/$run
done

mamba activate gatk4

#biohpc
# Cultivar haplotypes:
for sample in $( cat list_samples_culivarHap.txt); do echo $sample
for run in $( cat list_intervals.txt ); do echo $run
sbatch haplotypecaller2.sh $sample $run 1
done
done
# 04.05.2024: 


sample=SampleA_Hap_2
run=0001
sbatch haplotypecaller2.sh $sample $run 1

https://gatk.broadinstitute.org/hc/en-us/community/posts/360077582952-Haplotypecaller-produce-unfinished-vcf-but-no-idx-no-error

```



## Admixture



### Produce input - Plink

script: plink.sh 

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=3:00:00
#SBATCH -J Plink


chr=$1

VCF=/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/merged_vcf/chr/chr$chr"_"filterQUAL30.vcf.gz

FILE=chr$chr"_"filterQUAL30.LD0.3

#zcat $VCF | grep -v "##contig=<ID=scaffold" | awk -F"\t" 'BEGIN{ FS=OFS="\t" }{if($1 ~ /^chr/) $3="." ; print $0}' > tem_$FILE.vcf

vcftools --gzvcf $VCF --remove-indels --recode --recode-INFO-all --stdout | gzip -c > tem_$FILE.vcf.gz

plink --vcf tem_$FILE.vcf.gz --double-id --allow-extra-chr \
        --snps-only \
        --set-missing-var-ids @:# \
        --indep-pairwise 50 25 0.3 \
        --out $FILE --vcf-half-call 'missing'

plink --vcf tem_$FILE.vcf.gz --double-id --allow-extra-chr \
        --set-missing-var-ids @:# \
        --snps-only \
        --extract $FILE.prune.in \
        --vcf-half-call 'missing' \
        --make-bed --out $FILE

# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
cp $FILE.bim $FILE.bim_original
awk '{$1=0;print $0}' $FILE.bim > $FILE.bim.tmp
mv $FILE.bim.tmp $FILE.bim


```

Jobs:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/10_admixture

mamba activate env_others

for i in $(seq -w 1 12); do echo $i
sbatch plink.sh $i
done
# 31.05.2024: 2563099..2563110

#for i in {2563083..2563094}; do stop $i ; rm slurm-$i* ; done 
```



### Admixture runs

script: admixture_cv.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=1-00:00:00
#SBATCH -J admixture

# this environment contains plink, admixture
# conda activate env_others

chr=$1
k_value=$2
FILE=chr$chr"_"filterQUAL30.LD0.3

admixture --cv $FILE.bed $k_value > log${k_value}.$FILE.out


```

job:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/10_admixture

for rep in {1..5}; do echo $rep
mkdir -p replicate_$rep
cd replicate_$rep
cp ../chr* .
for k in {2..10}; do echo $k
for chr in $(seq -w 1 12); do echo $chr
sbatch ../admixture_cv.sh $chr $k 
done
done
cd ..
done
# 31.05.2024: 2564658..2565202

# to extract cross-validation error:
for chr in $(seq -w 1 12); do echo $chr
grep "CV" replicate_*/log*.chr$chr"_"*.out | awk -v Chr=$chr '{$5="chr"Chr ; print $5,$3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//' >> values.cv.error.txt
done


awk '{print $2"\t"$1}'  wildSamples.txt > SampleInformation.list
awk '{print "Cultivar\t"$1}' chr01_filterQUAL30.LD0.3.fam | grep Sample >> SampleInformation.list

for rep in {1..5}; do echo $rep
cp SampleInformation.list replicate_$rep"/"
done


#for i in {2564096..2564637}; do stop $i ; rm slurm-$i* ; done 
```

Plots:

Covariance error plot:

```R
# .libPaths("/dss/dsshome1/lxc03/di36guz2/R")
# mamba activate env_others
.libPaths()

setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/10_admixture")

library(ggplot2)
library(tidyverse)


error_values <- read_table("values.cv.error.txt",F)
names(error_values) <- c("chr", "k", "CVerror")
head(error_values)



error_values %>%
  ggplot(aes(factor(k), CVerror, group=chr, colour=chr)) +
  geom_point()+
  geom_line()+
  facet_grid(.~chr)+
  theme_classic()+
  theme(legend.position = "none")

ggsave(filename=paste0("plots/01_CVerror.png"), 
       device = png, units="cm", 
       width=35, height=4)
ggsave(filename=paste0("plots/01_CVerror.pdf"), 
       device = png, units="cm", 
       width=35, height=4)


```



Admixture plots:

Script: plotADMIXTURE.r

```R
#!/usr/bin/Rscript

# Usage: plotADMIXTURE.r -p <prefix> -i <info file, 2-column file with ind name and population/species name> 
#                        -k <max K value> -pop <comma-separated list of populations/species in the order to be plotted>
# This R script makes barplots for K=2 and all other K values until max K (specified with -k). It labels the individuals 
# and splits them into populations or species according to the individual and population/species names in the 2-column file specified with -i.
# The order of populations/species follows the list of populations/species given with -pop.
# Usage example: plotADMIXTURE.r -p fileXY -i file.ind.pop.txt -k 4 -pop pop1,pop2,pop3
# In this example, the script would use the files fileXY.2.Q, fileXY.3.Q, fileXY.4.Q to make barplots for the three populations.
# file.ind.pop.txt should contain one line for each individual in the same order as in the admixture files e.g.
# ind1 pop1
# ind2 pop1
# ind3 pop2
# ind4 pop3

# Author: Joana Meier, September 2019

# install.packages('optparse', lib="/home/sergio/R/libraries_Rackham/4.0")

# Read in the arguments
library("optparse")
library("dplyr")
library("wesanderson")

option_list = list(
  make_option(c("-p", "--prefix"), type="character", default=NULL, 
              help="prefix name (with path if not in the current directory)", metavar="character"),
  make_option(c("-i", "--infofile"), type="character", default=NULL, 
              help="info text file containing for each individual and extra", metavar="character"),
  make_option(c("-a", "--annotation"), type="character", default=NULL, 
              help="info text file containing for each individual the population/species information", metavar="character"),
  make_option(c("-e", "--excludeSamples"), type="character", default=NULL, 
              help="info text file containing for each individual the population/species information", metavar="character"),
  make_option(c("-k", "--maxK"), type="integer", default=NULL, 
              help="maximum K value", metavar="integer"),
  make_option(c("-m", "--minK"), type="integer", default=2, 
              help="minimum K value", metavar="integer"),
  make_option(c("-l", "--populations"), type="character", default=NULL, 
              help="comma-separated list of populations/species in the order to be plotted", metavar="character"),
  make_option(c("-o", "--outPrefix"), type="character", default="default", 
              help="output prefix (default: name provided with prefix)", metavar="character")
) 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Check that all required arguments are provided
if (is.null(opt$prefix)){
  print_help(opt_parser)
  stop("Please provide the prefix", call.=FALSE)
}else if (is.null(opt$infofile)){
  print_help(opt_parser)
  stop("Please provide the info file", call.=FALSE)
}else if (is.null(opt$maxK)){
  print_help(opt_parser)
  stop("Please provide the maximum K value to plot", call.=FALSE)
}else if (is.null(opt$populations)){
  print_help(opt_parser)
  stop("Please provide a comma-separated list of populations/species", call.=FALSE)
}

# If no output prefix is given, use the input prefix
if(opt$outPrefix=="default") opt$outPrefix=opt$prefix

# Assign the first argument to prefix
prefix=opt$prefix

# Get individual names in the correct order
labels_1<-read.table(opt$infofile)

# Name the columns
#names(labels)<-c("ind","pop")
names(labels_1)<-c("ind","extra")

annotation<-read.table(opt$annotation)
names(annotation)<-c("pop","ind")


#excludeSamples<-unlist(strsplit(opt$excludeSamples,","))

labels<-labels_1 %>%
  left_join(annotation, by="ind") %>%
  select(ind, pop) %>%
  data.frame()


# Add a column with population indices to order the barplots
# Use the order of populations provided as the fourth argument (list separated by commas)
# C_amarus,C_mucosospermus,C_lanatusE,C_cordophanus,C_lanatus_cordophanus,C_lanatus,C_lanatus_landrace,C_lanatus_cultivar,Other
labels$n<-factor(labels$pop,levels=unlist(strsplit(opt$populations,",")))
levels(labels$n)<-c(1:length(levels(labels$n)))
labels$n<-as.integer(as.character(labels$n))

labels2 <- labels %>%
  filter(!(ind %in% excludeSamples)) %>%
  data.frame()

# read in the different admixture output files
minK=opt$minK
maxK=opt$maxK
tbl<-lapply(minK:maxK, function(x) read.table(paste0(prefix,".",x,".Q")))

# # Prepare spaces to separate the populations/species
# rep<-as.vector(table(labels$n))
# spaces<-0
# for(i in 1:length(rep)){spaces=c(spaces,rep(0,rep[i]-1),0.5)}
# spaces<-spaces[-length(spaces)]

rep<-as.vector(table(labels2$n))
spaces<-0
for(i in 1:length(rep)){spaces=c(spaces,rep(0,rep[i]-1),0.5)}
spaces<-spaces[-length(spaces)]

# Plot the cluster assignments as a single bar for each individual for each K as a separate row
tiff(file=paste0(opt$outPrefix,".tiff"),width = 2000, height = 1200,res=200)
par(mfrow=c(maxK-1,1),mar=c(0,1,0,0),oma=c(2,1,9,1),mgp=c(0,0.2,0),xaxs="i",cex.lab=1.2,cex.axis=0.8)
# Plot minK
#bp<-barplot(t(as.matrix(tbl[[1]][order(labels$n),][!(labels$ind %in% excludeSamples),])), col=rainbow(n=minK),xaxt="n", border=NA,ylab=paste0("K=",minK),yaxt="n",space=spaces)
bp<-barplot(tbl[[1]] %>%
              cbind(labels) %>%
              #filter(!(labels$ind %in% excludeSamples)) %>%
              arrange(n) %>%
              select(-ind, -pop, -n) %>%
              as.matrix() %>%
              t(), 
            col=rainbow(n=minK), 
            #col=wes_palette("Darjeeling1", minK, type = "discrete"),
            xaxt="n", border=NA,ylab=paste0("K=",minK),yaxt="n",space=spaces)
axis(3,at=bp,labels= tbl[[1]] %>%
       cbind(labels) %>%
       filter(!(labels$ind %in% excludeSamples)) %>%
       arrange(n) %>%
       pull(ind) ,las=2,tick=F,cex=0.6)
# Plot higher K values
#if(maxK>minK)lapply(2:(maxK-1), function(x) barplot(t(as.matrix(tbl[[x]][order(labels$n),][!(labels$ind %in% excludeSamples),])), col=rainbow(n=x+1),xaxt="n", border=NA,ylab=paste0("K=",x+1),yaxt="n",space=spaces))
if(maxK>minK)lapply(2:(maxK-1), function(x) barplot(tbl[[x]] %>%
                                                      cbind(labels) %>%
                                                      filter(!(labels$ind %in% excludeSamples)) %>%
                                                      arrange(n) %>%
                                                      select(-ind, -pop, -n) %>%
                                                      as.matrix() %>%
                                                      t(), 
                                                    col=rainbow(n=x+1), 
                                                    #col=wes_palette("Darjeeling1", x+1, type = "discrete"),
                                                    xaxt="n", border=NA,ylab=paste0("K=",x+1),yaxt="n",space=spaces))
axis(1,at=c(which(spaces==0.5),bp[length(bp)])-diff(c(1,which(spaces==0.5),bp[length(bp)]))/2,
     labels=unlist(strsplit(opt$populations,",")))
dev.off()


```

Bash:

```bash

cd 

mamba activate env_others

chr=01
awk '{print $1"\t"$2}' chr$chr"_"filterQUAL30.LD0.3.fam > chr$chr"_"filterQUAL30.LD0.3.listSamples
Rscript ../plotADMIXTURE.r -p chr$chr"_"filterQUAL30.LD0.3 -i chr$chr"_"filterQUAL30.LD0.3.listSamples -a SampleInformation.list -k 2 -l Cultivar,C4_N,C4_S,C3,C1_2 -e "_"

for chr in $(seq -w 1 12); do echo $chr
awk '{print $1"\t"$2}' chr$chr"_"filterQUAL30.LD0.3.fam > chr$chr"_"filterQUAL30.LD0.3.listSamples
Rscript ../plotADMIXTURE.r -p chr$chr"_"filterQUAL30.LD0.3 -i chr$chr"_"filterQUAL30.LD0.3.listSamples -a SampleInformation.list -k 10 -l Cultivar,C4_N,C4_S,C3,C1_2 -e "_"
done



```

Plot V2:

```R

setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/06_mergedHaplotypes_PanGenome")

# Whole chromosomes with haplotype colours, but clustering haplotypes:
winSize="200kb"
stepSize="50kb"
minSS="50SS"
# chrID="12"
MaxGap="06"
n=10

singleL <- c("A","B","C","D","E","F","G","H","I","J")
AcName <- c("WhR", "BdF", "EgH", "PrW", "Flo", "Ack", "Fla", "Kat", "Lum", "EdB")

singleL_list <- c()
for (i in singleL){
  singleL_list <- c(singleL_list, rep(i, 48))
}

AcName_list <- c()
for (i in AcName){
  AcName_list <- c(AcName_list, rep(i, 48))
}

label_table <- data.frame(hapN=seq(1,48), singleL_list, AcName_list, simpleHapNum=rep(seq(1,4), 120)) %>%
  mutate(sample=paste0(singleL_list, "_hap", hapN)) %>%
  mutate(label=paste0(AcName_list, "_", simpleHapNum)) %>%
  select(sample, label)


singleL <- c("A","B","C","D","E","F","G","H","I","J")
AcName <- c("WhR", "BdF", "EgH", "PrW", "Flo", "Ack", "Fla", "Kat", "Lum", "EdB")

sample <- c()
for (i in singleL){
  for (hapN in seq(1,4)){
    sample <- c(sample, paste0(i, "_hap", hapN))
  }
}

labelS <- c()
for (i in AcName){
  for (hapN in seq(1,4)){
    labelS <- c(labelS, paste0(i, "_", hapN))
  }
}

label_tableD <- data.frame(sample, labelS)


singleL <- c("A","B","C","D","E","F","G","H","I","J")
Cultivar_samples <- c()
for (i in singleL){
  for (hapN in seq(1,4)){
    Cultivar_samples <- c(Cultivar_samples, paste0("Sample", i, "_Hap_", hapN))
  }
}


wildSamplesInfo <- read_tsv("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/08_distanceMatrix_perWin/WildSample_info_sim.txt", T)
wildSamplesInfo$Clade <- factor(wildSamplesInfo$Clade, levels=c("C4_N", "C4_S", "C4_S_C3", "C3_C4N", "C3", "C1_2"))
wildSamplesInfo

Cultivar_samples

head(label_table)

chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

chrID="12"

for (chrID in chr_list){
  
  # first we load the information of clustering on haplotypes to have the order of samples:
  
  setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/06_mergedHaplotypes_PanGenome")
  
  #Read in haplotypeNumberFile, with haplotype numbers reduced by the python script. 
  hapNums.in <- read.table(paste0("clusteringLargeHaplotypes_", winSize, "_step", stepSize,
                                  "_haplotypeGroups_min", minSS, "_chr", chrID, ".txt"), 
                           header=T, sep="\t")
  
  hapNums.in2 <- hapNums.in %>%
    left_join(label_table, by="sample") 
  
  matrix_hapNums.in <- hapNums.in2 %>%
    filter(Mean_NAN_haplo<0.8) %>%
    select(pos, label, haplotype) %>%
    spread(key=pos, value=haplotype, fill = NA) %>%
    select(-label) %>%
    as.matrix()
  
  samples_ID <- hapNums.in2 %>%
    filter(Mean_NAN_haplo<0.8) %>%
    select(pos, label, haplotype) %>%
    spread(key=pos, value=haplotype, fill = NA) %>%
    select(label) %>% unlist() %>% as.vector()
  
  d <- dist(matrix_hapNums.in, method = "euclidean")
  H.fit <- hclust(d, method="ward.D")
  H.fit$labels<-samples_ID
  
  # this line sets the order of samples:
  hapNums.in2$label <- factor(hapNums.in2$label, levels = H.fit$labels[H.fit$order])
  
  # Sequence similarity blocks:
  
  setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/08_distanceMatrix_perWin")
  
  values_WSp <- read_tsv(paste0("./distance_values/top90quantileSamples_dis_chr", chrID, ".txt"), F)
  names(values_WSp) <- c("chr", "win", "CulSample", "WildSample", "NTotalVar", "NVar", "support", "supportfq")
  # names(values_WSp) <- c("windowID","chr","win","end","mid","sites", "CulSample", "WildSample", "dis")
  
  # change Sample labels:
  
  # values_WSp %>%
  #   mutate(sample = str_replace_all(CulSample, c("^Sample" = "", "_Hap_" = "_hap"))) %>%
  #   left_join(label_tableD, by="sample") %>%
  #   mutate(label_or=factor(labelS, levels = H.fit$labels[H.fit$order])) %>%
  #   select(-labelS) %>%
  #   left_join(wildSamplesInfo, by="WildSample") %>%
  #   filter(supportfq>0.80) %>%
  #   group_by(chr, win, CulSample, label_or, Clade) %>%
  #   summarise(NSamples=n()) %>%
  #   group_by(chr, win, CulSample, label_or) %>%
  #   mutate(MaxNSamples=max(NSamples)) %>%
  #   filter(NSamples==MaxNSamples) %>%
  #   group_by(chr, win, CulSample, label_or) %>%
  #   summarise(NClades=length(unique(Clade)), 
  #             CladeG=Clade[1]) %>%
  #   filter(NClades==1) %>%
  #   # filter(CladeG!="C4_N") %>%
  #   ggplot(aes(x = win/1000000, y = label_or, fill = factor(CladeG))) +
  #   geom_tile() +
  #   # scale_fill_manual(values = c("#72b680","#7e73b6","#ae4f3b"))+
  #   scale_x_continuous(breaks=seq(0,100,5), expand = c(0,0), limits = c(0,90))+
  #   labs(title = paste0("Sequence similarity: Chr",chrID),
  #        x = "Position (Mb)",
  #        y = "Haplotype",
  #        fill = "Clade")+
  #   theme_classic()+
  #   theme(axis.text = element_text(size=12),
  #         axis.title = element_text(size=18),
  #         axis.line.x = element_blank(), 
  #         axis.title.x = element_text(hjust = 0.3))
  
  # ggsave(filename=paste0("plots/02_fromVCF_SeqSmilarity_perClade_chr", chrID, "_WinSize100kb_maxDis085.png"), 
  #        device = png, units="cm", 
  #        width=30, height=18)
  # ggsave(filename=paste0("plots/02_fromVCF_SeqSmilarity_perClade_chr", chrID, "_WinSize100kb_maxDis085.pdf"), 
  #        device = pdf, units="cm", 
  #        width=30, height=18)
  
  NSp_Cldes_table <- values_WSp %>%
    mutate(sample = str_replace_all(CulSample, c("^Sample" = "", "_Hap_" = "_hap"))) %>%
    left_join(label_tableD, by="sample") %>%
    mutate(label_or=factor(labelS, levels = H.fit$labels[H.fit$order])) %>%
    select(-labelS) %>%
    left_join(wildSamplesInfo, by="WildSample") %>%
    filter(supportfq>0.80) %>%
    group_by(chr, win, CulSample, label_or) %>%
    summarise(NSp=length(unique(Species)), NClades=length(unique(Clade))) 
  
  values_WSp %>%
    mutate(sample = str_replace_all(CulSample, c("^Sample" = "", "_Hap_" = "_hap"))) %>%
    left_join(label_tableD, by="sample") %>%
    mutate(label_or=factor(labelS, levels = H.fit$labels[H.fit$order])) %>%
    select(-labelS) %>%
    left_join(wildSamplesInfo, by="WildSample") %>%
    filter(supportfq>0.80) %>%
    left_join(NSp_Cldes_table, by=c("chr", "win", "CulSample", "label_or")) %>%
    filter(NSp<5) %>%
    filter(NClades<2) %>%
    group_by(chr, win, CulSample, label_or, Clade) %>%
    summarise(NSamples=n()) %>%
    group_by(chr, win, CulSample, label_or) %>%
    mutate(MaxNSamples=max(NSamples)) %>%
    filter(NSamples==MaxNSamples) %>%
    group_by(chr, win, CulSample, label_or) %>%
    summarise(NClades=length(unique(Clade)), 
              CladeG=Clade[1]) %>%
    filter(NClades==1) %>%
    # filter(CladeG!="C4_N") %>%
    ggplot(aes(x = win/1000000, y = label_or, fill = factor(CladeG), colour = factor(CladeG))) +
    geom_tile() +
    scale_fill_manual(values = c("#92e7a4", "#72b680","#7e73b6","#ae4f3b"))+
    scale_colour_manual(values = c("#92e7a4", "#72b680","#7e73b6","#ae4f3b"))+
    scale_x_continuous(breaks=seq(0,100,10), expand = c(0,0), limits = c(0,90))+
    labs(title = paste0("Sequence similarity: Chr",chrID),
         x = "Position (Mb)",
         y = "Haplotype",
         fill = "Clade")+
    theme_classic()+
    theme(axis.text = element_text(size=18, colour="black"),
          axis.title = element_text(size=20),
          title = element_text(size=20),
          legend.text = element_text(size=18),
          axis.line.x = element_blank(), 
          axis.title.x = element_text(hjust = 0.3))
  
  ggsave(filename=paste0("plots/02b_fromVCF_SeqSmilarity_perClade_chr", chrID, "_WinSize100kb_maxDis085.png"), 
         device = png, units="cm", 
         width=25, height=23)
  ggsave(filename=paste0("plots/02b_fromVCF_SeqSmilarity_perClade_chr", chrID, "_WinSize100kb_maxDis085.pdf"), 
         device = pdf, units="cm", 
         width=25, height=23)
  
}

```



### Whole genome run

script: plink_WG.sh 

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=3:00:00
#SBATCH -J Plink


VCF=$1

FILE=WG_filterQUAL30.LD0.3

#vcftools --gzvcf $VCF --remove-indels --recode --recode-INFO-all --stdout | gzip -c > tem_$FILE.vcf.gz

plink --vcf tem_$FILE.vcf.gz --double-id --allow-extra-chr \
        --snps-only \
        --set-missing-var-ids @:#:'$1':'$2' \
        --indep-pairwise 50 25 0.3 \
        --out $FILE --vcf-half-call 'missing'

plink --vcf tem_$FILE.vcf.gz --double-id --allow-extra-chr \
        --set-missing-var-ids @:#:'$1':'$2' \
        --snps-only \
        --extract $FILE.prune.in \
        --vcf-half-call 'missing' \
        --make-bed --out $FILE

# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
cp $FILE.bim $FILE.bim_original
awk '{$1=0;print $0}' $FILE.bim > $FILE.bim.tmp
mv $FILE.bim.tmp $FILE.bim

```

script: plink_WG_0.5.sh 

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=3:00:00
#SBATCH -J Plink


VCF=$1

FILE=WG_filterQUAL30.LD0.3
outFile=$2

#vcftools --gzvcf $VCF --remove-indels --recode --recode-INFO-all --stdout | gzip -c > tem_$FILE.vcf.gz

plink --vcf tem_$FILE.vcf.gz --double-id --allow-extra-chr \
        --snps-only \
        --geno 0.1 \
        --mac 2 \
        --set-missing-var-ids @:#:'$1':'$2' \
        --indep-pairwise 10kb 3000 0.5 \
        --out $outFile --vcf-half-call 'missing'

plink --vcf tem_$FILE.vcf.gz --double-id --allow-extra-chr \
        --set-missing-var-ids @:#:'$1':'$2' \
        --snps-only \
        --geno 0.1 \
        --mac 2 \
        --extract $outFile.prune.in \
        --vcf-half-call 'missing' \
        --make-bed --out $outFile

# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
cp $outFile.bim $outFile.bim_original
awk '{$1=0;print $0}' $outFile.bim > $outFile.bim.tmp
mv $outFile.bim.tmp $outFile.bim

```





Jobs:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/10_admixture

mamba activate env_others

sbatch plink_WG.sh /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/merged_vcf/wholegenome/WholeGenome_filterQUAL30.vcf.gz
# 07.06.2024: 2576055, 2576124


sbatch plink_WG_0.5.sh /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/merged_vcf/wholegenome/WholeGenome_filterQUAL30.vcf.gz WG_filterQUAL30.LD0.5
# 09.06.2024: 2576488

#for i in {2563083..2563094}; do stop $i ; rm slurm-$i* ; done 
```



script: admixture_cv_WG.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=2-00:00:00
#SBATCH -J admixture

# this environment contains plink, admixture
# conda activate env_others

k_value=$1
FILE=$2

admixture --cv $FILE.bed $k_value > log${k_value}.$FILE.out


```

Job:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/10_admixture

mamba activate env_others

for rep in {1..3}; do echo $rep
mkdir -p WG_replicate_$rep
cd WG_replicate_$rep
cp ../WG_* .
for k in {2..10}; do echo $k
sbatch ../admixture_cv_WG.sh $k WG_filterQUAL30.LD0.3
done
cd ..
done
# 08.07.2024:2576096..2576122


for rep in {1..3}; do echo $reps
mkdir -p LD0.5_WG_replicate_$rep
cd LD0.5_WG_replicate_$rep
cp ../WG_*.LD0.5* .
for k in {2..10}; do echo $k
sbatch ../admixture_cv_WG.sh $k WG_filterQUAL30.LD0.5
done
cd ..
done
# 09.07.2024:2576511..2576537



# to extract cross-validation error:
grep "CV" WG_replicate_*/log*.out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//' > WG_values.0.3.cv.error.txt

grep "CV" LD0.5_WG_replicate*/log*.out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//' > WG_values.0.5.cv.error.txt


awk '{print $2"\t"$1}'  wildSamples.txt > SampleInformation.list
awk '{print "Cultivar\t"$1}' WG_filterQUAL30.LD0.3.fam | grep Sample >> SampleInformation.list
cp SampleInformation.list WG_replicate_1/

awk '{print $2"\t"$1}'  ../wildSamples.txt > SampleInformation.list
awk '{print "Cultivar\t"$1}' WG_filterQUAL30.LD0.5.fam | grep Sample >> SampleInformation.list

#for i in {2564096..2564637}; do stop $i ; rm slurm-$i* ; done 
```

Plots:

```bash

cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/10_admixture/WG_replicate_1

awk '{print $1"\t"$2}' WG_filterQUAL30.LD0.3.fam > WG_filterQUAL30.LD0.3.listSamples
Rscript ../plotADMIXTURE.r -p WG_filterQUAL30.LD0.3 -i WG_filterQUAL30.LD0.3.listSamples -a SampleInformation.list -k 10 -l Cultivar,C4_N,C4_S,C3,C1_2 -e "_"


cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/10_admixture/LD0.5_WG_replicate_1

awk '{print $1"\t"$2}' WG_filterQUAL30.LD0.5.fam >  WG_filterQUAL30.LD0.5.listSamples
Rscript ../plotADMIXTURE.r -p  WG_filterQUAL30.LD0.5 -i WG_filterQUAL30.LD0.5.listSamples -a SampleInformation.list -k 3 -l Cultivar,C4_N,C4_S,C3,C1_2 -e "_"

```



## ABBA BABA



### DSuite



[Sources]: https://github.com/millanek/tutorials/tree/master/analysis_of_introgression_with_snp_data



#### **Test:**



```bash

cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/05_ABBA_BABA

mamba activate env_others


vcf=/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/splitVCF/01/split.008561.vcf.gz


Z_group="SampleA_Hap_1"
Y_group="SRR15458991"
W_group="SRR15458984"
X_group="SRR15458987"


bcftools query -l $vcf | awk '{print $1"\txxx"}' | sed 's/SampleE_Hap_4.*xx/SampleE_Hap_4\tSpecies1/g' | sed 's/SRR15458991.*xx/SRR15458991\tSpecies2/g' | sed 's/SRR15458984.*xx/SRR15458984\tSpecies3/g' | sed 's/SRR15458987.*xx/SRR15458987\tOutgroup/g' > species_sets.txt

#awk '{ if (NR <= 40) {sp=substr($1,1,3); print $1"\t"sp} else {print $1"\tOutgroup";}}' > species_sets.txt

/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/Dsuite/Build/Dsuite Dtrios -c -n testDSuite $vcf species_sets.txt 




vcf=/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/merged_vcf/chr/chr12_filterQUAL30.vcf.gz

#SampleE_Hap_4\tSpecies1
#SRR15458991.*xx/SRR15458991\tSpecies2
#SRR15458984.*xx/SRR15458984\tSpecies3
#SRR15458989.*xx/SRR15458989\tOutgroup
bcftools query -l $vcf | awk '{print $1"\txxx"}' | sed 's/SampleE_Hap_4.*xx/SampleE_Hap_4\tSampleE_Hap_4/g' | sed 's/SRR15458991.*xx/SRR15458991\tSRR15458991/g' | sed 's/SRR15458984.*xx/SRR15458984\tSRR15458984/g' | sed 's/SRR15458989.*xx/SRR15458989\tOutgroup/g' > species_sets.txt

/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/Dsuite/Build/Dsuite Dtrios -c -n chr12_E4_pir_Com_cho -t tree_Chr01_edCul.nwk $vcf species_sets.txt 

#SampleE_Hap_3\tSpecies1
#SRR15458991.*xx/SRR15458991\tSpecies2
#SRR15458984.*xx/SRR15458984\tSpecies3
#SRR15458989.*xx/SRR15458989\tOutgroup
bcftools query -l $vcf | awk '{print $1"\txxx"}' | sed 's/SampleE_Hap_3.*xx/SampleE_Hap_3\tSampleE_Hap_3/g' | sed 's/SRR15458991.*xx/SRR15458991\tSRR15458991/g' | sed 's/SRR15458984.*xx/SRR15458984\tSRR15458984/g' | sed 's/SRR15458989.*xx/SRR15458989\tOutgroup/g' > species_sets.txt

/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/Dsuite/Build/Dsuite Dtrios -c -n chr12_E3_pir_Com_cho -t tree_Chr01_edCul.nwk $vcf species_sets.txt 

#SampleE_Hap_2\tSpecies1
#SRR15458991.*xx/SRR15458991\tSpecies2
#SRR15458984.*xx/SRR15458984\tSpecies3
#SRR15458989.*xx/SRR15458989\tOutgroup
bcftools query -l $vcf | awk '{print $1"\txxx"}' | sed 's/SampleE_Hap_2.*xx/SampleE_Hap_2\tSampleE_Hap_2/g' | sed 's/SRR15458991.*xx/SRR15458991\tSRR15458991/g' | sed 's/SRR15458984.*xx/SRR15458984\tSRR15458984/g' | sed 's/SRR15458989.*xx/SRR15458989\tOutgroup/g' > species_sets.txt

/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/Dsuite/Build/Dsuite Dtrios -c -n chr12_E2_pir_Com_cho -t tree_Chr01_edCul.nwk $vcf species_sets.txt 

#SampleE_Hap_1\tSpecies1
#SRR15458991.*xx/SRR15458991\tSpecies2
#SRR15458984.*xx/SRR15458984\tSpecies3
#SRR15458989.*xx/SRR15458989\tOutgroup
bcftools query -l $vcf | awk '{print $1"\txxx"}' | sed 's/SampleE_Hap_1.*xx/SampleE_Hap_1\tSampleE_Hap_1/g' | sed 's/SRR15458991.*xx/SRR15458991\tSRR15458991/g' | sed 's/SRR15458984.*xx/SRR15458984\tSRR15458984/g' | sed 's/SRR15458989.*xx/SRR15458989\tOutgroup/g' > species_sets.txt

/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/Dsuite/Build/Dsuite Dtrios -c -n chr12_E1_pir_Com_cho -t tree_Chr01_edCul.nwk $vcf species_sets.txt 





Script=abba_baba.sh
chr=$1

vcf=/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/merged_vcf/chr/chr$chr"_"filterQUAL30.vcf.gz

#SampleE_Hap_4\tSpecies1
#SRR15458991.*xx/SRR15458991\tSpecies2
#SRR15458984.*xx/SRR15458984\tSpecies3
#SRR15458989.*xx/SRR15458989\tOutgroup
bcftools query -l $vcf | awk '{print $1"\txxx"}' | sed 's/SampleE_Hap_4.*xx/SampleE_Hap_4\tSampleE_Hap_4/g' | sed 's/SRR15458991.*xx/SRR15458991\tSRR15458991/g' | sed 's/SRR15458984.*xx/SRR15458984\tSRR15458984/g' | sed 's/SRR15458989.*xx/SRR15458989\tOutgroup/g' > species_sets.txt

/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/Dsuite/Build/Dsuite Dtrios -c -n chr$chr""_E4_pir_Com_cho -t tree_Chr01_edCul.nwk $vcf species_sets.txt 


bcftools query -l $vcf | awk '{print $1"\txxx"}' | sed 's/SampleE_Hap_4.*xx/SampleE_Hap_4\tSampleE_Hap_4/g' | sed 's/SRR15458991.*xx/SRR15458991\tSRR15458991/g' | sed 's/SRR15458984.*xx/SRR15458984\tSRR15458984/g' | sed 's/SRR15458989.*xx/SRR15458989\tOutgroup/g' > species_sets.txt
for chr in $(seq -w 1 11)
do echo $chr
sbatch abba_baba.sh $chr
done
# 2544602..2544612



SRR18082527
SRR5349609
SRR18082529









bcftools query -l $vcf | awk '{print $1"\txxx"}' | sed 's/SampleE_Hap_4.*xx/SampleE_Hap_4\tSpecies1/g' | sed 's/SRR15458991.*xx/SRR15458956\tSpecies2/g' | sed 's/SRR15458984.*xx/SRR15458984\tSpecies3/g' | sed 's/SRR15458989.*xx/SRR15458989\tOutgroup/g' > species_sets.txt

/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/Dsuite/Build/Dsuite Dtrios -c -n testDSuite_chr12_E4_buesii $vcf species_sets.txt 

/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/Dsuite/Build/Dsuite Dtrios -c -n chr12_E4_Com_lig_cho -t tree_Chr01_edCul.nwk $vcf species_sets.txt 


```

test with astral tree

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/05_ABBA_BABA/temporal3

mamba activate env_others

script:

vcf=/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/merged_vcf/chr/chr12_filterQUAL30.vcf.gz

bcftools query -l $vcf | awk '{print $1"\txxx"}' | sed 's/SampleE_Hap_4.*xx/SampleE_Hap_4\tSampleE_Hap_4/g' | sed 's/SRR15458957.*xx/SRR15458957\tSRR15458957/g' | sed 's/SRR15458984.*xx/SRR15458984\tSRR15458984/g' | sed 's/SRR15458989.*xx/SRR15458989\tOutgroup/g' > species_sets.txt

/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/Dsuite/Build/Dsuite Dtrios -c -n chr12_E4_lig_Com_cho -t tree_Chr01_edCul.nwk $vcf species_sets.txt 

#SampleE_Hap_3\tSpecies1
#SRR15458991.*xx/SRR15458991\tSpecies2
#SRR15458984.*xx/SRR15458984\tSpecies3
#SRR15458989.*xx/SRR15458989\tOutgroup
bcftools query -l $vcf | awk '{print $1"\txxx"}' | sed 's/SampleE_Hap_3.*xx/SampleE_Hap_3\tSampleE_Hap_3/g' | sed 's/SRR15458957.*xx/SRR15458957\tSRR15458957/g' | sed 's/SRR15458984.*xx/SRR15458984\tSRR15458984/g' | sed 's/SRR15458989.*xx/SRR15458989\tOutgroup/g' > species_sets.txt

/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/Dsuite/Build/Dsuite Dtrios -c -n chr12_E3_lig_Com_cho -t tree_Chr01_edCul.nwk $vcf species_sets.txt 

#SampleE_Hap_2\tSpecies1
#SRR15458991.*xx/SRR15458991\tSpecies2
#SRR15458984.*xx/SRR15458984\tSpecies3
#SRR15458989.*xx/SRR15458989\tOutgroup
bcftools query -l $vcf | awk '{print $1"\txxx"}' | sed 's/SampleE_Hap_2.*xx/SampleE_Hap_2\tSampleE_Hap_2/g' | sed 's/SRR15458957.*xx/SRR15458957\tSRR15458957/g' | sed 's/SRR15458984.*xx/SRR15458984\tSRR15458984/g' | sed 's/SRR15458989.*xx/SRR15458989\tOutgroup/g' > species_sets.txt

/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/Dsuite/Build/Dsuite Dtrios -c -n chr12_E2_lig_Com_cho -t tree_Chr01_edCul.nwk $vcf species_sets.txt 

#SampleE_Hap_1\tSpecies1
#SRR15458991.*xx/SRR15458991\tSpecies2
#SRR15458984.*xx/SRR15458984\tSpecies3
#SRR15458989.*xx/SRR15458989\tOutgroup
bcftools query -l $vcf | awk '{print $1"\txxx"}' | sed 's/SampleE_Hap_1.*xx/SampleE_Hap_1\tSampleE_Hap_1/g' | sed 's/SRR15458957.*xx/SRR15458957\tSRR15458957/g' | sed 's/SRR15458984.*xx/SRR15458984\tSRR15458984/g' | sed 's/SRR15458989.*xx/SRR15458989\tOutgroup/g' > species_sets.txt

/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/Dsuite/Build/Dsuite Dtrios -c -n chr12_E1_lig_Com_cho -t tree_Chr01_edCul.nwk $vcf species_sets.txt 



# new test: 24.05.2024:

cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/05_ABBA_BABA/temporal_withWGtree

mamba activate env_others

for chr in $(seq -w 1 12)
do echo $chr
sbatch script_Dsuite.sh $chr
done
#24.05.2024: 2547248..2547259

for job in $(seq 2 6); do echo $job
for chr in $(seq -w 1 12)
do echo $chr
sbatch script_Dsuite$job.sh $chr
done
done
#24.05.2024: 2547260..2547319

for chr in $(seq -w 1 12)
do echo $chr
sbatch script_Dsuite7.sh $chr
done
#24.05.2024: 2547358..2547369


for i in {2551486..2551965}; do stop $i ; rm slurm-$i.out ; done

2547171..2547244
```





#### **Per Chromosome:**

Script: ABBA_BABA_perChr.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=12:00:00


#chr=12
#sample1=SampleE_Hap_4
#sample2=SRR15458993
#sample3=SRR15458986
#sample4=SRR15458989
#output_prefix=bre_vem_Cul_cho_E4

chr=$1
sample1=$2
sample2=$3
sample3=$4
sample4=$5
output_prefix=$6
folder=$7

mkdir -p $folder

vcf=/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/merged_vcf/chr/chr$chr"_"filterQUAL30.vcf.gz

bcftools query -l $vcf | awk '{print $1"\txxx"}' | \
    sed "s/${sample1}.*xx/${sample1}\t${sample1}/g" | \
    sed "s/${sample2}.*xx/${sample2}\t${sample2}/g" | \
    sed "s/${sample3}.*xx/${sample3}\t${sample3}/g" | \
    sed "s/${sample4}.*xx/${sample4}\tOutgroup/g" > species_sets_${chr}_${output_prefix}.txt

echo -e "chr\tP1\tP2\tP3\tDstatistic\tZ-score\tp-value\tf4-ratio\tBBAA\tABBA\tBABA" > Values_Chr$chr"_"$output_prefix.txt

# delete files from previous window
if test -f ./species_sets_${chr}_${output_prefix}_chr$chr"_"$output_prefix"_tree.txt" ; then
  rm species_sets_${chr}_${output_prefix}_chr$chr"_"$output_prefix"_"*.txt
fi

# calculate D values
/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/Dsuite/Build/Dsuite Dtrios -c -n chr$chr"_"$output_prefix -t astral.merged_trees.WGenome.min25.nwk $vcf species_sets_${chr}_${output_prefix}.txt

# add values to output file:
grep "" species_sets_${chr}_${output_prefix}_chr$chr"_"$output_prefix"_tree.txt" | grep -v "ABBA" | awk -v chrN="$chr" 'BEGIN {OFS="\t"} {print chrN, $0}' >> Values_Chr$chr"_"$output_prefix.txt

rm species_sets_${chr}_${output_prefix}.txt

mv Values_Chr$chr"_"$output_prefix.txt $folder"/"

```

Job:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/05_ABBA_BABA/DSuite_perChr

mamba activate env_others


vcf=/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/merged_vcf/chr/chr12_filterQUAL30.vcf.gz

bcftools query -l $vcf | grep Sample > list_Cul_samples.txt

for chr in $(seq -w 1 12)
do echo $chr
for cul in $( cat list_Cul_samples.txt ); do echo $cul
shortCul=$(echo $cul | sed 's/Sample//g' | sed 's/_Hap_//g')
sample2=SRR15458993
sample3=SRR15458986
sample4=SRR15458989
folder=bre_ver_Cul_cho
output_prefix=$folder"_"$shortCul
sbatch ./ABBA_BABA_perChr.sh $chr $cul $sample2 $sample3 $sample4 $output_prefix $folder
done
done
# 29.05.2024: 2554861..2555340

for chr in $(seq -w 1 12)
do echo $chr
for cul in $( cat list_Cul_samples.txt ); do echo $cul
shortCul=$(echo $cul | sed 's/Sample//g' | sed 's/_Hap_//g')
sample2=SRR15458993
sample3=SRR15458982
sample4=SRR15458989
folder=bre_bol_Cul_cho
output_prefix=$folder"_"$shortCul
sbatch ./ABBA_BABA_perChr.sh $chr $cul $sample2 $sample3 $sample4 $output_prefix $folder
done
done
# 29.05.2024: 2555341..2555870


for chr in $(seq -w 1 12)
do echo $chr
for cul in $( cat list_Cul_samples.txt ); do echo $cul
shortCul=$(echo $cul | sed 's/Sample//g' | sed 's/_Hap_//g')
sample2=SRR15458993
sample3=SRR15458984
sample4=SRR15458989
folder=bre_com_Cul_cho
output_prefix=$folder"_"$shortCul
sbatch ./ABBA_BABA_perChr.sh $chr $cul $sample2 $sample3 $sample4 $output_prefix $folder
done
done
# 29.05.2024: 2556431..2556910

for chr in $(seq -w 1 12)
do echo $chr
for cul in $( cat list_Cul_samples.txt ); do echo $cul
shortCul=$(echo $cul | sed 's/Sample//g' | sed 's/_Hap_//g')
sample2=SRR15458984
sample3=SRR15458982
sample4=SRR15458989
folder=com_bol_Cul_cho
output_prefix=$folder"_"$shortCul
sbatch ./ABBA_BABA_perChr.sh $chr $cul $sample2 $sample3 $sample4 $output_prefix $folder
done
done
# 29.05.2024: 



for chr in $(seq -w 1 12)
do echo $chr
for cul in $( cat list_Cul_samples.txt ); do echo $cul
shortCul=$(echo $cul | sed 's/Sample//g' | sed 's/_Hap_//g')
sample2=SRR15458993
sample3=SRR15458989
sample4=SRR15458990
folder=Cul_bre_choma_more
output_prefix=$folder"_"$shortCul
sbatch ./ABBA_BABA_perChr.sh $chr $cul $sample2 $sample3 $sample4 $output_prefix $folder
done
done
# 29.05.2024: 



for chr in $(seq -w 1 12)
do echo $chr
for cul in $( cat list_Cul_samples.txt ); do echo $cul
shortCul=$(echo $cul | sed 's/Sample//g' | sed 's/_Hap_//g')
sample2=SRR15458993
sample3=SRR15458980
sample4=SRR15458990
folder=Cul_bre_caj_more
output_prefix=$folder"_"$shortCul
sbatch ./ABBA_BABA_perChr.sh $chr $cul $sample2 $sample3 $sample4 $output_prefix $folder
done
done
# 29.05.2024: 



for chr in $(seq -w 1 12)
do echo $chr
for cul in $( cat list_Cul_samples.txt ); do echo $cul
shortCul=$(echo $cul | sed 's/Sample//g' | sed 's/_Hap_//g')
sample2=SRR15458993
sample3=SRR15458959
sample4=SRR15458990
folder=Cul_bre_and_more
output_prefix=$folder"_"$shortCul
sbatch ./ABBA_BABA_perChr.sh $chr $cul $sample2 $sample3 $sample4 $output_prefix $folder
done
done
# 29.05.2024: 2555882..2556362

```



##### Plots:

```R

setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/05_ABBA_BABA/DSuite_perChr")

# # Whole chromosomes with haplotype colours, but clustering haplotypes:
# winSize="200kb"
# stepSize="50kb"
# minSS="50SS"
# # chrID="12"
# MaxGap="06"
# n=10

# singleL <- c("A","B","C","D","E","F","G","H","I","J")
# AcName <- c("WhR", "BdF", "EgH", "PrW", "Flo", "Ack", "Fla", "Kat", "Lum", "EdB")
# 
# singleL_list <- c()
# for (i in singleL){
#   singleL_list <- c(singleL_list, rep(i, 48))
# }
# 
# AcName_list <- c()
# for (i in AcName){
#   AcName_list <- c(AcName_list, rep(i, 48))
# }
# 
# label_table <- data.frame(hapN=seq(1,48), singleL_list, AcName_list, simpleHapNum=rep(seq(1,4), 120)) %>%
#   mutate(sample=paste0(singleL_list, "_hap", hapN)) %>%
#   mutate(label=paste0(AcName_list, "_", simpleHapNum)) %>%
#   select(sample, label)
# 

singleL <- c("A","B","C","D","E","F","G","H","I","J")
AcName <- c("WhR", "BdF", "EgH", "PrW", "Flo", "Ack", "Fla", "Kat", "Lum", "EdB")

sample <- c()
for (i in singleL){
  for (hapN in seq(1,4)){
    sample <- c(sample, paste0(i, "_hap", hapN))
  }
}

labelS <- c()
for (i in AcName){
  for (hapN in seq(1,4)){
    labelS <- c(labelS, paste0(i, "_", hapN))
  }
}

label_tableD <- data.frame(sample, labelS)


singleL <- c("A","B","C","D","E","F","G","H","I","J")
Cultivar_samples <- c()
for (i in singleL){
  for (hapN in seq(1,4)){
    Cultivar_samples <- c(Cultivar_samples, paste0("Sample", i, "_Hap_", hapN))
  }
}


Cultivar_samples

head(label_tableD)



### Introgression from C4_S

case="bre_ver_Cul_cho"

testDone <- c("bre_bol_Cul_cho", 
              "bre_com_Cul_cho", 
              "bre_ver_Cul_cho", 
              "com_bol_Cul_cho")

chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

chrID="12"

if (exists("merged_combined_data2")) {
  remove(merged_combined_data2)
}


for (case in testDone){
  
  for (chrID in chr_list){
    
    # setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/06_mergedHaplotypes_PanGenome")
    # 
    # #Read in haplotypeNumberFile, with haplotype numbers reduced by the python script. 
    # hapNums.in <- read.table(paste0("clusteringLargeHaplotypes_", winSize, "_step", stepSize,
    #                                 "_haplotypeGroups_min", minSS, "_chr", chrID, ".txt"), 
    #                          header=T, sep="\t")
    # 
    # hapNums.in2 <- hapNums.in %>%
    #   left_join(label_table, by="sample") 
    # 
    # matrix_hapNums.in <- hapNums.in2 %>%
    #   filter(Mean_NAN_haplo<0.8) %>%
    #   select(pos, label, haplotype) %>%
    #   spread(key=pos, value=haplotype, fill = NA) %>%
    #   select(-label) %>%
    #   as.matrix()
    # 
    # samples_ID <- hapNums.in2 %>%
    #   filter(Mean_NAN_haplo<0.8) %>%
    #   select(pos, label, haplotype) %>%
    #   spread(key=pos, value=haplotype, fill = NA) %>%
    #   select(label) %>% unlist() %>% as.vector()
    # 
    # d <- dist(matrix_hapNums.in, method = "euclidean")
    # H.fit <- hclust(d, method="ward.D")
    # H.fit$labels<-samples_ID
    # 
    # # this line sets the order of samples:
    # hapNums.in2$label <- factor(hapNums.in2$label, levels = H.fit$labels[H.fit$order])
    # 

    
    # Introgression blocks:
    
    setwd(paste0("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/05_ABBA_BABA/DSuite_perChr/", case, "/"))
    
    # Define the pattern to match the files
    file_pattern <- paste0("Values_Chr", chrID)
    
    # List all files matching the pattern
    file_list <- list.files(pattern = file_pattern)
    
    # Read and combine all files into a single dataframe
    data_list <- lapply(file_list, read_tsv)
    combined_data <- bind_rows(data_list)
    
    combined_data2 <- combined_data %>%
      mutate(case=case) %>%
      mutate(Clade="C4_S")
    
    if (exists("merged_combined_data2")) {
      merged_combined_data2 <- merged_combined_data2 %>%
        rbind(combined_data2)
    } else {
      merged_combined_data2 <- combined_data2
    }
  }
}
    

setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/05_ABBA_BABA/DSuite_perChr")

# write.table(merged_combined_data2, "Values_perChr_FStats_C4S.txt", sep = "\t", 
#             col.names = T, row.names = F, quote = F)

merged_combined_data2 <- read_table("Values_perChr_FStats_C4S.txt", T)

# change Sample labels:
    
merged_combined_data2 %>%
  mutate(sample = str_replace_all(P3, c("^Sample" = "", "_Hap_" = "_hap"))) %>%
  left_join(label_tableD, by="sample") %>%
  mutate(sig=ifelse(`p-value`< 0.001, "***", ifelse(`p-value`< 0.01,"**", ifelse(`p-value`< 0.05,"*", "")))) %>%
  # mutate(label_or=factor(labelS, levels = H.fit$labels[H.fit$order])) %>%
  # select(-labelS) %>%
  # filter((ABBA/BBAA)>0.2) %>%
  # filter((`f4-ratio`)>0.3) %>%
  # filter(`Z-score`>3) %>%
  ggplot(aes(x = labelS, y = `f4-ratio`, fill=case)) +
  geom_bar(stat="identity", position="dodge")+
  geom_text(aes(label = sig, y=0.4), 
            vjust = 0.9, 
            colour="gray50",
            position = position_dodge(.9), angle=90)+
  # scale_fill_viridis_c(direction= -1) +
  # scale_x_continuous(breaks=seq(0,100,5), expand = c(0,0), limits = c(0,90))+
  scale_y_continuous(breaks=seq(0,1,0.2), expand = c(0,0), limits = c(0,0.6))+
  labs(x = "Haplotype",
       y = "f4",
       fill = "f4")+
  facet_grid(chr~.)+
  theme_minimal()+
  # theme_classic()+
  theme(axis.text = element_text(size=12), 
        axis.text.x = element_text(angle = 90),
        axis.title = element_text(size=18),
        axis.line.x = element_blank())
# legend.position = "none")
    
ggsave(filename=paste0("./plots/01_F4Stat_perChr_Cases_C4S.png"),
           device = png, units="cm",
           width=30, height=30)
ggsave(filename=paste0("./plots/01_F4Stat_perChr_Cases_C4S.pdf"),
           device = pdf, units="cm",
           width=30, height=30)
    
    



### Introgression from C3

case="Cul_bre_and_more"

testDone <- c("Cul_bre_and_more", 
              "Cul_bre_caj_more", 
              "Cul_bre_choma_more")

chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

chrID="12"
# {

if (exists("merged_combined_data3")) {
  remove(merged_combined_data3)
}


for (case in testDone){
  
  for (chrID in chr_list){
    
    # setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/06_mergedHaplotypes_PanGenome")
    # 
    # #Read in haplotypeNumberFile, with haplotype numbers reduced by the python script. 
    # hapNums.in <- read.table(paste0("clusteringLargeHaplotypes_", winSize, "_step", stepSize,
    #                                 "_haplotypeGroups_min", minSS, "_chr", chrID, ".txt"), 
    #                          header=T, sep="\t")
    # 
    # hapNums.in2 <- hapNums.in %>%
    #   left_join(label_table, by="sample") 
    # 
    # matrix_hapNums.in <- hapNums.in2 %>%
    #   filter(Mean_NAN_haplo<0.8) %>%
    #   select(pos, label, haplotype) %>%
    #   spread(key=pos, value=haplotype, fill = NA) %>%
    #   select(-label) %>%
    #   as.matrix()
    # 
    # samples_ID <- hapNums.in2 %>%
    #   filter(Mean_NAN_haplo<0.8) %>%
    #   select(pos, label, haplotype) %>%
    #   spread(key=pos, value=haplotype, fill = NA) %>%
    #   select(label) %>% unlist() %>% as.vector()
    # 
    # d <- dist(matrix_hapNums.in, method = "euclidean")
    # H.fit <- hclust(d, method="ward.D")
    # H.fit$labels<-samples_ID
    # 
    # # this line sets the order of samples:
    # hapNums.in2$label <- factor(hapNums.in2$label, levels = H.fit$labels[H.fit$order])
    # 
    
    
    # Introgression blocks:
    
    setwd(paste0("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/05_ABBA_BABA/DSuite_perChr/", case, "/"))
    
    # Define the pattern to match the files
    file_pattern <- paste0("Values_Chr", chrID)
    
    # List all files matching the pattern
    file_list <- list.files(pattern = file_pattern)
    
    # Read and combine all files into a single dataframe
    data_list <- lapply(file_list, read_tsv)
    combined_data <- bind_rows(data_list)
    
    combined_data2 <- combined_data %>%
      filter(P2 %in% Cultivar_samples) %>%
      mutate(case=case) %>%
      mutate(Clade="C3")
    
    if (exists("merged_combined_data3")) {
      merged_combined_data3 <- merged_combined_data3 %>%
        rbind(combined_data2)
    } else {
      merged_combined_data3 <- combined_data2
    }
  }
}



setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/05_ABBA_BABA/DSuite_perChr")

write.table(merged_combined_data3, "Values_perChr_FStats_C3.txt", sep = "\t",
            col.names = T, row.names = F, quote = F)

# merged_combined_data3 <- read_table("Values_perChr_FStats_C3.txt", T)

# change Sample labels:

merged_combined_data3 %>%
  mutate(sample = str_replace_all(P2, c("^Sample" = "", "_Hap_" = "_hap"))) %>%
  left_join(label_tableD, by="sample") %>%
  mutate(sig=ifelse(`p-value`< 0.001, "***", ifelse(`p-value`< 0.01,"**", ifelse(`p-value`< 0.05,"*", "")))) %>%
  # mutate(label_or=factor(labelS, levels = H.fit$labels[H.fit$order])) %>%
  # select(-labelS) %>%
  # filter((ABBA/BBAA)>0.2) %>%
  # filter((`f4-ratio`)>0.3) %>%
  # filter(`Z-score`>3) %>%
  ggplot(aes(x = labelS, y = `f4-ratio`, fill=case)) +
  geom_bar(stat="identity", position="dodge")+
  geom_text(aes(label = sig, y=0.4), 
            vjust = 0.9, 
            colour="gray60",
            position = position_dodge(.9), angle=90)+
  # scale_fill_viridis_c(direction= -1) +
  # scale_x_continuous(breaks=seq(0,100,5), expand = c(0,0), limits = c(0,90))+
  scale_y_continuous(breaks=seq(0,1,0.2), expand = c(0,0), limits = c(0,0.6))+
  labs(x = "Haplotype",
       y = "f4",
       fill = "f4")+
  facet_grid(chr~.)+
  theme_minimal()+
  # theme_classic()+
  theme(axis.text = element_text(size=12), 
        axis.text.x = element_text(angle = 90),
        axis.title = element_text(size=18),
        axis.line.x = element_blank())
# legend.position = "none")
  
ggsave(filename=paste0("./plots/02_F4Stat_perChr_Cases_C3.png"),
       device = png, units="cm",
       width=30, height=30)
ggsave(filename=paste0("./plots/02_F4Stat_perChr_Cases_C3.pdf"),
       device = pdf, units="cm",
       width=30, height=30)



```





#### Per Window

Script: Dsuite_perWin.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00


#chr=12
#sample1=SampleE_Hap_4
#sample2=SRR15458993
#sample3=SRR15458986
#sample4=SRR15458989
#output_prefix=bre_vem_Cul_cho_E4

chr=$1
sample1=$2
sample2=$3
sample3=$4
sample4=$5
output_prefix=$6
folder=$7

mkdir -p $folder

ls -1 /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/splitVCF/$chr"/"*.vcf.gz > list_vcf_Chr$chr"_"$output_prefix.txt

vcf=$(head -1 list_vcf_Chr$chr"_"$output_prefix.txt)

bcftools query -l $vcf | awk '{print $1"\txxx"}' | \
    sed "s/${sample1}.*xx/${sample1}\t${sample1}/g" | \
    sed "s/${sample2}.*xx/${sample2}\t${sample2}/g" | \
    sed "s/${sample3}.*xx/${sample3}\t${sample3}/g" | \
    sed "s/${sample4}.*xx/${sample4}\tOutgroup/g" > species_sets_${chr}_${output_prefix}.txt

# delete pre existing output file:
#if test -f ./Values_Chr$chr"_"$output_prefix.txt ; then
#  rm Values_Chr$chr"_"$output_prefix.txt
#fi

echo -e "chr\twin\tP1\tP2\tP3\tDstatistic\tZ-score\tp-value\tf4-ratio\tBBAA\tABBA\tBABA" > Values_Chr$chr"_"$output_prefix.txt


for vcf in $( cat list_vcf_Chr$chr"_"$output_prefix.txt ); do echo $vcf
# Extract the filename
filename="${vcf##*/}"
# Extract the number using parameter expansion
number="${filename#split.}"
window="${number%%.*}"
# delete files from previous window
if test -f ./species_sets_${chr}_${output_prefix}_chr$chr"_"$output_prefix"_tree.txt" ; then
  rm species_sets_${chr}_${output_prefix}_chr$chr"_"$output_prefix"_"*.txt
fi
# calculate D values
/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/Dsuite/Build/Dsuite Dtrios -c -n chr$chr"_"$output_prefix -t astral.merged_trees.WGenome.min25.nwk $vcf species_sets_${chr}_${output_prefix}.txt
# add values to output file:
grep "" species_sets_${chr}_${output_prefix}_chr$chr"_"$output_prefix"_tree.txt" | grep -v "ABBA" | awk -v win="$window" -v chrN="$chr" 'BEGIN {OFS="\t"} {print chrN, win, $0}' >> Values_Chr$chr"_"$output_prefix.txt
done

rm species_sets_${chr}_${output_prefix}.txt list_vcf_Chr$chr"_"$output_prefix.txt

mv Values_Chr$chr"_"$output_prefix.txt $folder"/"


```







```bash

cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/05_ABBA_BABA/

mamba activate env_others

vcf=/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/merged_vcf/chr/chr12_filterQUAL30.vcf.gz

bcftools query -l $vcf | grep Sample > list_Cul_samples.txt

for chr in $(seq -w 1 12)
do echo $chr
for cul in $( cat list_Cul_samples.txt ); do echo $cul
shortCul=$(echo $cul | sed 's/Sample//g' | sed 's/_Hap_//g')
sample2=SRR15458993
sample3=SRR15458986
sample4=SRR15458989
folder=bre_ver_Cul_cho
output_prefix=$folder"_"$shortCul
sbatch ./Dsuite_perWin.sh $chr $cul $sample2 $sample3 $sample4 $output_prefix $folder
done
done
# 24.05.2024: 2547962..2548441, 2587152..2587633

for chr in $(seq -w 1 12)
do echo $chr
for cul in $( cat list_Cul_samples.txt ); do echo $cul
shortCul=$(echo $cul | sed 's/Sample//g' | sed 's/_Hap_//g')
sample2=SRR15458993
sample3=SRR15458982
sample4=SRR15458989
folder=bre_bol_Cul_cho
output_prefix=$folder"_"$shortCul
sbatch ./Dsuite_perWin.sh $chr $cul $sample2 $sample3 $sample4 $output_prefix $folder
done
done
# 24.05.2024: 2548442..2548921

#for i in {4105986..4106225}; do stopserial $i ; done
for i in {2587635..2588120}; do stop $i ; done


for chr in $(seq -w 1 12)
do echo $chr
for cul in $( cat list_Cul_samples.txt ); do echo $cul
shortCul=$(echo $cul | sed 's/Sample//g' | sed 's/_Hap_//g')
sample2=SRR15458993
sample3=SRR15458984
sample4=SRR15458989
folder=bre_com_Cul_cho
output_prefix=$folder"_"$shortCul
sbatch ./Dsuite_perWin.sh $chr $cul $sample2 $sample3 $sample4 $output_prefix $folder
done
done
# 24.05.2024: 2548922..2549401, 2587635..2588120

for chr in $(seq -w 1 12)
do echo $chr
for cul in $( cat list_Cul_samples.txt ); do echo $cul
shortCul=$(echo $cul | sed 's/Sample//g' | sed 's/_Hap_//g')
sample2=SRR15458984
sample3=SRR15458982
sample4=SRR15458989
folder=com_bol_Cul_cho
output_prefix=$folder"_"$shortCul
sbatch ./Dsuite_perWin.sh $chr $cul $sample2 $sample3 $sample4 $output_prefix $folder
done
done
# 24.05.2024: 2549402..2549881



for chr in $(seq -w 1 12)
do echo $chr
for cul in $( cat list_Cul_samples.txt ); do echo $cul
shortCul=$(echo $cul | sed 's/Sample//g' | sed 's/_Hap_//g')
sample2=SRR15458993
sample3=SRR15458989
sample4=SRR15458990
folder=Cul_bre_choma_more
output_prefix=$folder"_"$shortCul
sbatch ./Dsuite_perWin.sh $chr $cul $sample2 $sample3 $sample4 $output_prefix $folder
done
done
# 24.05.2024: 2551967..2552446



for chr in $(seq -w 1 12)
do echo $chr
for cul in $( cat list_Cul_samples.txt ); do echo $cul
shortCul=$(echo $cul | sed 's/Sample//g' | sed 's/_Hap_//g')
sample2=SRR15458993
sample3=SRR15458980
sample4=SRR15458990
folder=Cul_bre_caj_more
output_prefix=$folder"_"$shortCul
sbatch ./Dsuite_perWin.sh $chr $cul $sample2 $sample3 $sample4 $output_prefix $folder
done
done
# 24.05.2024: 2552447..2552926



for chr in $(seq -w 1 12)
do echo $chr
for cul in $( cat list_Cul_samples.txt ); do echo $cul
shortCul=$(echo $cul | sed 's/Sample//g' | sed 's/_Hap_//g')
sample2=SRR15458993
sample3=SRR15458959
sample4=SRR15458990
folder=Cul_bre_and_more
output_prefix=$folder"_"$shortCul
sbatch ./Dsuite_perWin.sh $chr $cul $sample2 $sample3 $sample4 $output_prefix $folder
done
done
# 24.05.2024: 2552927..2553406

```



##### Plots:



```R

setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/05_ABBA_BABA/DSuite")

library(ggplot2)
library(tidyverse)

# Whole chromosomes with haplotype colours, but clustering haplotypes:

{winSize="200kb"
  stepSize="50kb"
  minSS="50SS"
  # chrID="12"
  MaxGap="06"
  n=10
  
  singleL <- c("A","B","C","D","E","F","G","H","I","J")
  AcName <- c("WhR", "BdF", "EgH", "PrW", "Flo", "Ack", "Fla", "Kat", "Lum", "EdB")
  
  singleL_list <- c()
  for (i in singleL){
    singleL_list <- c(singleL_list, rep(i, 48))
  }
  
  AcName_list <- c()
  for (i in AcName){
    AcName_list <- c(AcName_list, rep(i, 48))
  }
  
  label_table <- data.frame(hapN=seq(1,48), singleL_list, AcName_list, simpleHapNum=rep(seq(1,4), 120)) %>%
    mutate(sample=paste0(singleL_list, "_hap", hapN)) %>%
    mutate(label=paste0(AcName_list, "_", simpleHapNum)) %>%
    select(sample, label)
  
  
  singleL <- c("A","B","C","D","E","F","G","H","I","J")
  AcName <- c("WhR", "BdF", "EgH", "PrW", "Flo", "Ack", "Fla", "Kat", "Lum", "EdB")
  
  sample <- c()
  for (i in singleL){
    for (hapN in seq(1,4)){
      sample <- c(sample, paste0(i, "_hap", hapN))
    }
  }
  
  labelS <- c()
  for (i in AcName){
    for (hapN in seq(1,4)){
      labelS <- c(labelS, paste0(i, "_", hapN))
    }
  }
  
  label_tableD <- data.frame(sample, labelS)
  
  
  singleL <- c("A","B","C","D","E","F","G","H","I","J")
  Cultivar_samples <- c()
  for (i in singleL){
    for (hapN in seq(1,4)){
      Cultivar_samples <- c(Cultivar_samples, paste0("Sample", i, "_Hap_", hapN))
    }
  }
  
  
  Cultivar_samples
  
  head(label_table)
}



### Introgression from C4_S

case="bre_bol_Cul_cho"

testDone <- c("bre_bol_Cul_cho", 
              "bre_com_Cul_cho", 
              "bre_ver_Cul_cho", 
              "com_bol_Cul_cho")

chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

chrID="12"
# {

for (case in testDone){
  
  for (chrID in chr_list){
    
    setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/06_mergedHaplotypes_PanGenome")
    
    #Read in haplotypeNumberFile, with haplotype numbers reduced by the python script. 
    hapNums.in <- read.table(paste0("clusteringLargeHaplotypes_", winSize, "_step", stepSize,
                                    "_haplotypeGroups_min", minSS, "_chr", chrID, ".txt"), 
                             header=T, sep="\t")
    
    hapNums.in2 <- hapNums.in %>%
      left_join(label_table, by="sample") 
    
    matrix_hapNums.in <- hapNums.in2 %>%
      filter(Mean_NAN_haplo<0.8) %>%
      select(pos, label, haplotype) %>%
      spread(key=pos, value=haplotype, fill = NA) %>%
      select(-label) %>%
      as.matrix()
    
    samples_ID <- hapNums.in2 %>%
      filter(Mean_NAN_haplo<0.8) %>%
      select(pos, label, haplotype) %>%
      spread(key=pos, value=haplotype, fill = NA) %>%
      select(label) %>% unlist() %>% as.vector()
    
    d <- dist(matrix_hapNums.in, method = "euclidean")
    H.fit <- hclust(d, method="ward.D")
    H.fit$labels<-samples_ID
    
    # this line sets the order of samples:
    hapNums.in2$label <- factor(hapNums.in2$label, levels = H.fit$labels[H.fit$order])
    
    # Min_NWin=5
    # Min_SharedHap=3
    # 
    # MinLenHaplotypeID <- hapNums.in2 %>%
    #   filter(Mean_NAN_haplo<0.8) %>%
    #   group_by(sample, haplotype) %>%
    #   summarise(NWin=n()) %>%
    #   filter(NWin>Min_NWin) %>%
    #   ungroup() %>%
    #   pull(haplotype) %>%
    #   unique()
    # 
    # MinSharedLenHaplotypeID <- hapNums.in2 %>%
    #   filter(Mean_NAN_haplo<0.8) %>%
    #   filter(haplotype %in% MinLenHaplotypeID) %>%
    #   group_by(haplotype) %>%
    #   summarise(NSamples=length(unique(sample))) %>%
    #   filter(NSamples>Min_SharedHap) %>%
    #   ungroup() %>%
    #   pull(haplotype) %>%
    #   unique()
    # 
    # 
    # hapNums.in2 %>%
    #   # filter(Mean_NAN_haplo<MaxGap) %>%
    #   # filter(haplotype %in% MinSharedLenHaplotypeID) %>%
    #   ggplot(aes(pos/1000000, label, fill=factor(haplotype) , colour=factor(haplotype))) +
    #   geom_tile()+
    #   scale_color_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
    #                                     hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
    #   scale_fill_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
    #                                    hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
    #   scale_x_continuous(breaks=seq(0,100,5), name = "Position (Mb)", expand = c(0,0), limits = c(0,90))+
    #   scale_y_discrete(name="Haplotype")+
    #   ggtitle(paste0("Chr:", chrID))+
    #   theme_classic()+
    #   theme(axis.text = element_text(size=12),
    #         axis.title = element_text(size=18),
    #         axis.line.x = element_blank(),
    #         legend.position = "none")
    
    
    # Introgression blocks:
    
    setwd(paste0("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/05_ABBA_BABA/DSuite/", case, "/"))
    
    # Define the pattern to match the files
    file_pattern <- paste0("Values_Chr", chrID)
    
    # List all files matching the pattern
    file_list <- list.files(pattern = file_pattern)
    
    # Read and combine all files into a single dataframe
    data_list <- lapply(file_list, read_tsv)
    combined_data <- bind_rows(data_list)
    
    # Convert 'win' to numeric
    combined_data <- combined_data %>%
      mutate(win = as.numeric(win))
    
    # change Sample labels:
    
    combined_data %>%
      mutate(sample = str_replace_all(P3, c("^Sample" = "", "_Hap_" = "_hap"))) %>%
      left_join(label_tableD, by="sample") %>%
      mutate(label_or=factor(labelS, levels = H.fit$labels[H.fit$order])) %>%
      select(-labelS) %>%
      filter(`p-value`<0.01) %>%
      filter((ABBA/BBAA)>0.2) %>%
      filter((`f4-ratio`)>0.1) %>%
      # filter(`Z-score`>3) %>%
      ggplot(aes(x = win*100000/1000000, y = label_or, fill = `f4-ratio`)) +
      geom_tile() +
      scale_fill_viridis_c(direction= -1) +
      scale_x_continuous(breaks=seq(0,100,10), expand = c(0,0), limits = c(0,90))+
      labs(title = paste0("D Test:", case, " f4-ratio Chr",chrID),
           x = "Pos (MB)",
           y = "Haplotype",
           fill = "f4")+
      theme_classic()+
      theme(axis.text = element_text(size=18, colour="black"),
            axis.title = element_text(size=20),
            title = element_text(size=18),
            legend.text = element_text(size=18),
            axis.line.x = element_blank(), 
            axis.title.x = element_text(hjust = 0.3))
    
    # legend.position = "none")
    
    ggsave(filename=paste0("../plots/01_Case_", case, "_Chr", chrID, "_f4_minf4_03_minABBA_BBAA_02_Win100kb.png"),
           device = png, units="cm", 
           width=25, height=23)
    ggsave(filename=paste0("../plots/01_Case_", case, "_Chr", chrID, "_f4_minf4_03_minABBA_BBAA_02_Win100kb.pdf"),
           device = pdf, units="cm", 
           width=25, height=23)
    
    combined_data %>%
      mutate(sample = str_replace_all(P3, c("^Sample" = "", "_Hap_" = "_hap"))) %>%
      left_join(label_tableD, by="sample") %>%
      mutate(label_or=factor(labelS, levels = H.fit$labels[H.fit$order])) %>%
      select(-labelS) %>%
      filter(`p-value`<0.01) %>%
      filter((ABBA/BBAA)>0.2) %>%
      filter((`f4-ratio`)>0.1) %>%
      # filter((`f4-ratio`)>0.3) %>%
      # filter(`Z-score`>3) %>%
      filter(Dstatistic>0.4)%>%
      ggplot(aes(x = win*100000/1000000, y = label_or, fill = Dstatistic)) +
      geom_tile() +
      scale_fill_viridis_c(direction= -1) +
      scale_x_continuous(breaks=seq(0,100,10), expand = c(0,0), limits = c(0,90))+
      labs(title = paste0("Test:", case, " Dstatistic Chr",chrID),
           x = "Pos (MB)",
           y = "Haplotype",
           fill = "D")+
      theme_classic()+
      theme(axis.text = element_text(size=18, colour="black"),
            axis.title = element_text(size=20),
            title = element_text(size=18),
            legend.text = element_text(size=18),
            axis.line.x = element_blank(), 
            axis.title.x = element_text(hjust = 0.3))
    # legend.position = "none")
    
    ggsave(filename=paste0("../plots/02_Case_", case, "_Chr", chrID, "_D_minD_04_minABBA_BBAA_02_Win100kb.png"),
           device = png, units="cm", 
           width=25, height=23)
    ggsave(filename=paste0("../plots/02_Case_", case, "_Chr", chrID, "_D_minD_04_minABBA_BBAA_02_Win100kb.pdf"),
           device = pdf, units="cm", 
           width=25, height=23)
    
  }
}





### Introgression from C3



head(label_table)

case="Cul_bre_and_more"

testDone <- c("Cul_bre_and_more", 
              "Cul_bre_caj_more", 
              "Cul_bre_choma_more")

chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

chrID="12"
# {

for (case in testDone){
  #   
  for (chrID in chr_list){
    
    setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/06_mergedHaplotypes_PanGenome")
    
    #Read in haplotypeNumberFile, with haplotype numbers reduced by the python script. 
    hapNums.in <- read.table(paste0("clusteringLargeHaplotypes_", winSize, "_step", stepSize,
                                    "_haplotypeGroups_min", minSS, "_chr", chrID, ".txt"), 
                             header=T, sep="\t")
    
    hapNums.in2 <- hapNums.in %>%
      left_join(label_table, by="sample") 
    
    matrix_hapNums.in <- hapNums.in2 %>%
      filter(Mean_NAN_haplo<0.8) %>%
      select(pos, label, haplotype) %>%
      spread(key=pos, value=haplotype, fill = NA) %>%
      select(-label) %>%
      as.matrix()
    
    samples_ID <- hapNums.in2 %>%
      filter(Mean_NAN_haplo<0.8) %>%
      select(pos, label, haplotype) %>%
      spread(key=pos, value=haplotype, fill = NA) %>%
      select(label) %>% unlist() %>% as.vector()
    
    d <- dist(matrix_hapNums.in, method = "euclidean")
    H.fit <- hclust(d, method="ward.D")
    H.fit$labels<-samples_ID
    
    # this line sets the order of samples:
    hapNums.in2$label <- factor(hapNums.in2$label, levels = H.fit$labels[H.fit$order])
    
    # Min_NWin=5
    # Min_SharedHap=3
    # 
    # MinLenHaplotypeID <- hapNums.in2 %>%
    #   filter(Mean_NAN_haplo<0.8) %>%
    #   group_by(sample, haplotype) %>%
    #   summarise(NWin=n()) %>%
    #   filter(NWin>Min_NWin) %>%
    #   ungroup() %>%
    #   pull(haplotype) %>%
    #   unique()
    # 
    # MinSharedLenHaplotypeID <- hapNums.in2 %>%
    #   filter(Mean_NAN_haplo<0.8) %>%
    #   filter(haplotype %in% MinLenHaplotypeID) %>%
    #   group_by(haplotype) %>%
    #   summarise(NSamples=length(unique(sample))) %>%
    #   filter(NSamples>Min_SharedHap) %>%
    #   ungroup() %>%
    #   pull(haplotype) %>%
    #   unique()
    # 
    # 
    # hapNums.in2 %>%
    #   # filter(Mean_NAN_haplo<MaxGap) %>%
    #   # filter(haplotype %in% MinSharedLenHaplotypeID) %>%
    #   ggplot(aes(pos/1000000, label, fill=factor(haplotype) , colour=factor(haplotype))) +
    #   geom_tile()+
    #   scale_color_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
    #                                     hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
    #   scale_fill_manual(values = rep(c(hcl.colors(n, palette = "Berlin"), hcl.colors(n, palette = "Fall"), 
    #                                    hcl.colors(n, palette = "Sunset"), hcl.colors(n, palette = "Dynamic")), 800))+
    #   scale_x_continuous(breaks=seq(0,100,5), name = "Position (Mb)", expand = c(0,0), limits = c(0,90))+
    #   scale_y_discrete(name="Haplotype")+
    #   ggtitle(paste0("Chr:", chrID))+
    #   theme_classic()+
    #   theme(axis.text = element_text(size=12),
    #         axis.title = element_text(size=18),
    #         axis.line.x = element_blank(),
    #         legend.position = "none")
    
    
    # Introgression blocks:
    
    setwd(paste0("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/05_ABBA_BABA/DSuite/", case, "/"))
    
    # Define the pattern to match the files
    file_pattern <- paste0("Values_Chr", chrID)
    
    # List all files matching the pattern
    file_list <- list.files(pattern = file_pattern)
    
    # Read and combine all files into a single dataframe
    data_list <- lapply(file_list, read_tsv)
    combined_data <- bind_rows(data_list)
    
    # Convert 'win' to numeric
    combined_data <- combined_data %>%
      mutate(win = as.numeric(win))
    
    # change Sample labels:
    
    combined_data %>%
      filter(P2 %in% Cultivar_samples) %>%
      mutate(sample = str_replace_all(P2, c("^Sample" = "", "_Hap_" = "_hap"))) %>%
      left_join(label_tableD, by="sample") %>%
      mutate(label_or=factor(labelS, levels = H.fit$labels[H.fit$order])) %>%
      select(-labelS) %>%
      filter(`p-value`<0.01) %>%
      filter((ABBA/BBAA)>0.2) %>%
      filter((`f4-ratio`)>0.1) %>%
      # filter((ABBA/BBAA)>0.2) %>%
      # filter((`f4-ratio`)>0.2) %>%
      ggplot(aes(x = win*100000/1000000, y = label_or, fill = `f4-ratio`)) +
      geom_tile() +
      scale_fill_viridis_c(direction= -1) +
      scale_x_continuous(breaks=seq(0,100,10), expand = c(0,0), limits = c(0,90))+
      labs(title = paste0("D Test:", case, " f4-ratio Chr",chrID),
           x = "Pos (MB)",
           y = "Haplotype",
           fill = "f4")+
      theme_classic()+
      theme(axis.text = element_text(size=18, colour="black"),
            axis.title = element_text(size=20),
            title = element_text(size=18),
            legend.text = element_text(size=18),
            axis.line.x = element_blank(), 
            axis.title.x = element_text(hjust = 0.3))
    # legend.position = "none")
    
    ggsave(filename=paste0("../plots/03_Case_", case, "_Chr", chrID, "_f4_minf4_02_minABBA_BBAA_02_Win100kb.png"),
           device = png, units="cm",
           width=25, height=23)
    ggsave(filename=paste0("../plots/03_Case_", case, "_Chr", chrID, "_f4_minf4_02_minABBA_BBAA_02_Win100kb.pdf"),
           device = pdf, units="cm",
           width=25, height=23)
    
    combined_data %>%
      filter(P2 %in% Cultivar_samples) %>%
      mutate(sample = str_replace_all(P2, c("^Sample" = "", "_Hap_" = "_hap"))) %>%
      left_join(label_tableD, by="sample") %>%
      mutate(label_or=factor(labelS, levels = H.fit$labels[H.fit$order])) %>%
      select(-labelS) %>%
      filter(`p-value`<0.01) %>%
      filter((ABBA/BBAA)>0.2) %>%
      filter((`f4-ratio`)>0.1) %>%
      # filter((`f4-ratio`)>0.3) %>%
      # filter(`Z-score`>3) %>%
      filter(Dstatistic>0.4)%>%
      ggplot(aes(x = win*100000/1000000, y = label_or, fill = Dstatistic)) +
      geom_tile() +
      scale_fill_viridis_c(direction= -1) +
      scale_x_continuous(breaks=seq(0,100,10), expand = c(0,0), limits = c(0,90))+
      labs(title = paste0("Test:", case, " Dstatistic Chr",chrID),
           x = "Pos (MB)",
           y = "Haplotype",
           fill = "D")+
      theme_classic()+
      theme(axis.text = element_text(size=18, colour="black"),
            axis.title = element_text(size=20),
            title = element_text(size=18),
            legend.text = element_text(size=18),
            axis.line.x = element_blank(), 
            axis.title.x = element_text(hjust = 0.3))
    # legend.position = "none")
    
    ggsave(filename=paste0("../plots/04_Case_", case, "_Chr", chrID, "_D_minD_04_minABBA_BBAA_02_Win100kb.png"),
           device = png, units="cm",
           width=25, height=23)
    ggsave(filename=paste0("../plots/04_Case_", case, "_Chr", chrID, "_D_minD_04_minABBA_BBAA_02_Win100kb.pdf"),
           device = pdf, units="cm",
           width=25, height=23)
    
  }
}














#  Combine information from multiple test:
setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/05_ABBA_BABA/DSuite/")

chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")


case="bre_bol_Cul_cho"

chrID="12"

testDone_C4_S <- c("bre_bol_Cul_cho", 
                   "bre_com_Cul_cho", 
                   "bre_ver_Cul_cho", 
                   "com_bol_Cul_cho")

case="Cul_bre_and_more"

testDone_C3 <- c("Cul_bre_and_more", 
                 "Cul_bre_caj_more", 
                 "Cul_bre_choma_more")


for (chrID in chr_list){

file_path_C3 <- paste0("merged_filtered_table_C3_chr", chrID, ".csv")
levels_file_path_C3 <- paste0("factor_levels_C3_chr", chrID, ".rds")

# Define the file paths for merged_filtered_table_C4_S
file_path_C4_S <- paste0("merged_filtered_table_C4_S_chr", chrID, ".csv")
levels_file_path_C4_S <- paste0("factor_levels_C4_S_chr", chrID, ".rds")

# Function to save the data frame and factor levels
save_data_with_levels <- function(data, data_path, levels_path) {
  # Save the data frame
  write.csv(data, file = data_path, row.names = FALSE)
  
  # Save the levels of factors
  factor_levels <- lapply(data, function(column) {
    if (is.factor(column)) {
      levels(column)
    } else {
      NULL
    }
  })
  saveRDS(factor_levels, file = levels_path)
}

# Function to restore factor levels
restore_factor_levels <- function(data, levels_path) {
  factor_levels <- readRDS(levels_path)
  for (col_name in names(factor_levels)) {
    if (!is.null(factor_levels[[col_name]])) {
      data[[col_name]] <- factor(data[[col_name]], levels = factor_levels[[col_name]])
    }
  }
  return(data)
}

setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/05_ABBA_BABA/DSuite/")

# Process merged_filtered_table_C4_S if file does not exist
if (!file.exists(file_path_C4_S)) {
  merged_filtered_table_C4_S <- data.frame(test="Notest", chr="00", pos=0, cul_Sample="SampleA_Hap_1", Clade="C4_S")
  
  for (case in testDone_C4_S) {
    setwd(paste0("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/05_ABBA_BABA/DSuite/", case, "/"))
    file_pattern <- paste0("Values_Chr", chrID)
    file_list <- list.files(pattern = file_pattern)
    data_list <- lapply(file_list, read_tsv)
    combined_data <- bind_rows(data_list)
    combined_data <- combined_data %>%
      mutate(win = as.numeric(win))
    
    filtered_table <- combined_data %>%
      filter(`p-value`<0.01) %>%
      filter((ABBA/BBAA)>0.2) %>%
      filter((`f4-ratio`)>0.2) %>%
      mutate(Clade="C4_S") %>%
      mutate(cul_Sample=P3) %>%
      mutate(pos=win*100000) %>%
      mutate(test=case) %>%
      select(test, chr, pos, cul_Sample, Clade)
    
    merged_filtered_table_C4_S <- merged_filtered_table_C4_S %>%
      rbind(filtered_table)
  }
  
  merged_filtered_table_C4_S <- merged_filtered_table_C4_S %>%
    filter(chr != "00")
  
  setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/05_ABBA_BABA/DSuite/")
  
  # Save the data frame and factor levels
  save_data_with_levels(merged_filtered_table_C4_S, file_path_C4_S, levels_file_path_C4_S)
} else {
  # Load existing data and restore factor levels
  merged_filtered_table_C4_S <- read.csv(file_path_C4_S, stringsAsFactors = FALSE)
  merged_filtered_table_C4_S <- restore_factor_levels(merged_filtered_table_C4_S, levels_file_path_C4_S)
}

setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/05_ABBA_BABA/DSuite/")

# Process merged_filtered_table_C3 if file does not exist
if (!file.exists(file_path_C3)) {
merged_filtered_table_C3 <- data.frame(test="Notest", chr="00", pos=0, cul_Sample="SampleA_Hap_1", Clade="C3")

for (case in testDone_C3) {
  setwd(paste0("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/05_ABBA_BABA/DSuite/", case, "/"))
  file_pattern <- paste0("Values_Chr", chrID)
  file_list <- list.files(pattern = file_pattern)
  data_list <- lapply(file_list, read_tsv)
  combined_data <- bind_rows(data_list)
  combined_data <- combined_data %>%
    mutate(win = as.numeric(win))
  
  filtered_table <- combined_data %>%
    filter(`p-value`<0.01) %>%
    filter((ABBA/BBAA)>0.2) %>%
    filter((`f4-ratio`)>0.2) %>%
    filter(P2 %in% Cultivar_samples) %>%
    mutate(Clade="C3") %>%
    mutate(cul_Sample=P2) %>%
    mutate(pos=win*100000) %>%
    mutate(test=case) %>%
    select(test, chr, pos, cul_Sample, Clade)
  
  merged_filtered_table_C3 <- merged_filtered_table_C3 %>%
    rbind(filtered_table)
}

merged_filtered_table_C3 <- merged_filtered_table_C3 %>%
  filter(chr != "00")

setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/05_ABBA_BABA/DSuite/")

# Save the data frame and factor levels
save_data_with_levels(merged_filtered_table_C3, file_path_C3, levels_file_path_C3)
} else {
  # Load existing data and restore factor levels
  merged_filtered_table_C3 <- read.csv(file_path_C3, stringsAsFactors = FALSE)
  merged_filtered_table_C3 <- restore_factor_levels(merged_filtered_table_C3, levels_file_path_C3)
}

# calculate order of haplotypes: 

setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/06_mergedHaplotypes_PanGenome")

#Read in haplotypeNumberFile, with haplotype numbers reduced by the python script. 
hapNums.in <- read.table(paste0("clusteringLargeHaplotypes_", winSize, "_step", stepSize,
                                "_haplotypeGroups_min", minSS, "_chr", chrID, ".txt"), 
                         header=T, sep="\t")

hapNums.in2 <- hapNums.in %>%
  left_join(label_table, by="sample") 

matrix_hapNums.in <- hapNums.in2 %>%
  filter(Mean_NAN_haplo<0.8) %>%
  select(pos, label, haplotype) %>%
  spread(key=pos, value=haplotype, fill = NA) %>%
  select(-label) %>%
  as.matrix()

samples_ID <- hapNums.in2 %>%
  filter(Mean_NAN_haplo<0.8) %>%
  select(pos, label, haplotype) %>%
  spread(key=pos, value=haplotype, fill = NA) %>%
  select(label) %>% unlist() %>% as.vector()

d <- dist(matrix_hapNums.in, method = "euclidean")
H.fit <- hclust(d, method="ward.D")
H.fit$labels<-samples_ID

# this line sets the order of samples:
hapNums.in2$label <- factor(hapNums.in2$label, levels = H.fit$labels[H.fit$order])


# produce plot:

setwd(paste0("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/05_ABBA_BABA/DSuite"))

merged_filtered_table_C4_S %>%
  rbind(merged_filtered_table_C3) %>%
  mutate(sample = str_replace_all(cul_Sample, c("^Sample" = "", "_Hap_" = "_hap"))) %>%
  left_join(label_tableD, by="sample") %>%
  mutate(label_or=factor(labelS, levels = H.fit$labels[H.fit$order])) %>%
  select(-labelS) %>%
  group_by(chr, pos, cul_Sample, Clade, sample, label_or) %>%
  summarise(SupportingTest=n()) %>%
  # filter(SupportingTest>1) %>%
  group_by(chr, pos, cul_Sample, sample, label_or) %>%
  summarise(NClades=n(), 
            CladeSingle=Clade[1]) %>%
  filter(NClades==1) %>%
  mutate(CladeG=factor(CladeSingle, levels=c("C4_S", "C3"))) %>%
  # filter(`Z-score`>3) %>%
  ggplot(aes(x = pos/1000000, y = label_or, fill = CladeG)) +
  geom_tile() +
  # scale_fill_viridis_c(direction= -1) +
  scale_fill_manual(values = c("#72b680","#7e73b6"))+
  scale_x_continuous(breaks=seq(0,100,10), expand = c(0,0), limits = c(0,90))+
  labs(title = paste0("ABBA-BABA support: Chr",chrID),
       x = "Pos (MB)",
       y = "Haplotype",
       fill = "Clade")+
  theme_classic()+
  theme(axis.text = element_text(size=18, colour="black"),
        axis.title = element_text(size=20),
        title = element_text(size=20),
        legend.text = element_text(size=18),
        axis.line.x = element_blank(), 
        axis.title.x = element_text(hjust = 0.3))

ggsave(filename=paste0("./plots/05_ABBA_BABA_support_Combined_Chr", chrID, "_Win100kb.png"),
       device = png, units="cm",
       width=25, height=23)
ggsave(filename=paste0("./plots/05_ABBA_BABA_support_Combined_Chr", chrID, "_Win100kb.pdf"),
       device = pdf, units="cm",
       width=25, height=23)

# Function to calculate the number of lines per clade in overlapping windows
calculate_windows <- function(df, window_size = 100000, step_size = 20000) {
  result <- list()
  
  for (chr in unique(df$chr)) {
    chr_df <- df %>% filter(chr == chr)
    max_pos <- max(chr_df$pos)
    
    for (start in seq(0, max_pos, by = step_size)) {
      end <- start + window_size
      
      window_df <- chr_df %>%
        filter(pos >= start & pos < end) %>%
        group_by(cul_Sample, sample, label_or, CladeSingle) %>%
        summarize(count = n(), .groups = 'drop')
      
      if (nrow(window_df) > 0) {
        window_df <- window_df %>%
          mutate(chr = chr, start = start, end = end)
        
        result <- append(result, list(window_df))
      }
    }
  }
  
  final_result <- bind_rows(result)
  return(final_result)
}

file_path_C4_S_sim <- paste0("merged_filtered_table_C4_S_C3_sim_chr", chrID, ".csv")
levels_file_path_C4_S_sim <- paste0("factor_levels_C4_S_C3_sim_chr", chrID, ".rds")

# Process merged_filtered_table_C4_S_sim if file does not exist
if (!file.exists(file_path_C4_S_sim)) {
winSizeInt <- 1000000
merged_filtered_table_C4_S_sim <- merged_filtered_table_C4_S %>%
  rbind(merged_filtered_table_C3) %>%
  mutate(sample = str_replace_all(cul_Sample, c("^Sample" = "", "_Hap_" = "_hap"))) %>%
  left_join(label_tableD, by="sample") %>%
  mutate(label_or=factor(labelS, levels = H.fit$labels[H.fit$order])) %>%
  select(-labelS) %>%
  group_by(chr, pos, cul_Sample, Clade, sample, label_or) %>%
  summarise(SupportingTest=n()) %>%
  # filter(SupportingTest>1) %>%
  group_by(chr, pos, cul_Sample, sample, label_or) %>%
  summarise(NClades=n(), 
            CladeSingle=Clade[1]) %>%
  filter(NClades==1) %>%
  calculate_windows(window_size = winSizeInt, step_size = 200000) 
# Save the data frame and factor levels
save_data_with_levels(merged_filtered_table_C4_S_sim, file_path_C4_S_sim, levels_file_path_C4_S_sim)
} else {
  # Load existing data and restore factor levels
  merged_filtered_table_C4_S_sim <- read.csv(file_path_C4_S_sim, stringsAsFactors = FALSE)
  merged_filtered_table_C4_S_sim <- restore_factor_levels(merged_filtered_table_C4_S_sim, levels_file_path_C4_S_sim)
}


winSizeInt=1000000
min_block=0.3*(winSizeInt/100000)

# merged_filtered_table_C4_S_sim %>%
#   filter(count>min_block) %>%
#   mutate(pos=start) %>%
#   mutate(CladeG=factor(CladeSingle, levels=c("C4_S", "C3"))) %>%
#   # filter(`Z-score`>3) %>%
#   ggplot(aes(x = pos/1000000, y = label_or, fill = CladeG, colour=CladeG)) +
#   geom_tile() +
#   # scale_fill_viridis_c(direction= -1) +
#   scale_fill_manual(values = c("#72b680","#7e73b6"))+
#   scale_colour_manual(values = c("#72b680","#7e73b6"))+
#   scale_x_continuous(breaks=seq(0,100,10), expand = c(0,0), limits = c(0,90))+
#   labs(title = paste0("ABBA-BABA support: Chr",chrID),
#        x = "Pos (MB)",
#        y = "Haplotype",
#        fill = "Clade")+
#   theme_classic()+
#   theme(axis.text = element_text(size=18, colour="black"),
#         axis.title = element_text(size=20),
#         title = element_text(size=20),
#         legend.text = element_text(size=18),
#         axis.line.x = element_blank(), 
#         axis.title.x = element_text(hjust = 0.3))
# 
# ggsave(filename=paste0("./plots/06_ABBA_BABA_support_Combined_Chr", chrID, "_Win100kb_sim1MbWinStep200kb.png"),
#        device = png, units="cm",
#        width=25, height=23)
# ggsave(filename=paste0("./plots/06_ABBA_BABA_support_Combined_Chr", chrID, "_Win100kb_sim1MbWinStep200kb.pdf"),
#        device = pdf, units="cm",
#        width=25, height=23)

min_block=0.4*(winSizeInt/100000)

merged_filtered_table_C4_S_sim %>%
  filter(count>min_block) %>%
  mutate(pos=start) %>%
  mutate(CladeG=factor(CladeSingle, levels=c("C4_S", "C3"))) %>%
  # filter(`Z-score`>3) %>%
  ggplot(aes(x = pos/1000000, y = label_or, fill = CladeG, colour=CladeG)) +
  geom_tile() +
  # scale_fill_viridis_c(direction= -1) +
  scale_fill_manual(values = c("#72b680","#7e73b6"))+
  scale_colour_manual(values = c("#72b680","#7e73b6"))+
  scale_x_continuous(breaks=seq(0,100,10), expand = c(0,0), limits = c(0,90))+
  labs(title = paste0("Introgression f4 support: Chr",chrID),
       x = "Pos (MB)",
       y = "Haplotype",
       fill = "Clade")+
  theme_classic()+
  theme(axis.text = element_text(size=18, colour="black"),
        axis.title = element_text(size=20),
        title = element_text(size=20),
        legend.text = element_text(size=18),
        axis.line.x = element_blank(),
        axis.title.x = element_text(hjust = 0.3))

ggsave(filename=paste0("./plots/07_ABBA_BABA_support_Combined_Chr", chrID, "_Win100kb_sim1MbWinStep200kb.png"),
       device = png, units="cm",
       width=25, height=23)
ggsave(filename=paste0("./plots/07_ABBA_BABA_support_Combined_Chr", chrID, "_Win100kb_sim1MbWinStep200kb.pdf"),
       device = pdf, units="cm",
       width=25, height=23)


}






# proportions per Cultivar:

chrID="12"
# Define the list of chromosomes
chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

label_table2 <- data.frame(hapN=seq(1,48), singleL_list, AcName_list, simpleHapNum=rep(seq(1,4), 120)) %>%
  mutate(sample=paste0(singleL_list, "_hap", hapN)) %>%
  mutate(label=paste0(AcName_list, "_", simpleHapNum)) %>%
  select(sample, label, AcName_list)

chrID="12"

# Initialize an empty data frame to store the combined results


# Loop through each chromosome and process the data
combined_summary_table_Cul <- data.frame() ; for (chrID in chr_list) {
  
  # Define file paths for merged_filtered_table_C3
  file_path_C3 <- paste0("merged_filtered_table_C3_chr", chrID, ".csv")
  
  # Define file paths for merged_filtered_table_C4_S
  file_path_C4_S <- paste0("merged_filtered_table_C4_S_chr", chrID, ".csv")
  
  file_path_C4_S_sim <- paste0("merged_filtered_table_C4_S_C3_sim_chr", chrID, ".csv")
  
  # Load and restore factor levels for merged_filtered_table_C4_S
  merged_filtered_table_C4_S <- read.csv(file_path_C4_S)
  
  # Load and restore factor levels for merged_filtered_table_C3
  merged_filtered_table_C3 <- read.csv(file_path_C3)
  
  merged_filtered_table_C4_S_sim <- read.csv(file_path_C4_S_sim, stringsAsFactors = FALSE)
  
  # # Process and combine data
  # summary_table <- merged_filtered_table_C4_S %>%
  #   rbind(merged_filtered_table_C3) %>%
  #   mutate(sample = str_replace_all(cul_Sample, c("^Sample" = "", "_Hap_" = "_hap"))) %>%
  #   left_join(label_table2, by="sample") %>%
  #   group_by(chr, pos, cul_Sample, Clade, sample, label, AcName_list) %>%
  #   summarise(SupportingTest=n(), .groups = 'drop') %>%
  #   group_by(chr, pos, cul_Sample, sample, label, AcName_list) %>%
  #   summarise(NClades=n(), CladeSingle=Clade[1], .groups = 'drop') %>%
  #   filter(NClades==1) %>%
  #   group_by(chr, CladeSingle, AcName_list) %>%
  #   summarise(BlockSize=n()*100000, .groups = 'drop')
  
  winSizeInt <- 1000000
  min_block=0.2*(winSizeInt/100000)
  
  summary_table <- merged_filtered_table_C4_S_sim %>%
    left_join(label_table2 %>%
                select(sample, AcName_list), by="sample") %>%
    filter(count>min_block) %>%
    mutate(pos=start) %>%
    group_by(chr, CladeSingle, AcName_list) %>%
    summarise(BlockSize=n()*200000, .groups = 'drop')
  
  # Combine the summary table with the combined summary table
  combined_summary_table_Cul <- bind_rows(combined_summary_table_Cul, summary_table)
}

# combined_summary_table now contains the combined data from all chromosomes


combined_summary_table_Cul %>%
  group_by(AcName_list, CladeSingle) %>%
  summarise(TotalSize=sum(BlockSize)) %>%
  mutate(CladeG=factor(CladeSingle, levels=c("C4_S", "C3"))) %>%
  ggplot(aes(x = AcName_list, y = TotalSize/1000000, fill = CladeG)) +
  # geom_bar(stat = "identity", position = position_dodge()) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("#72b680","#7e73b6"))+
  scale_y_continuous(breaks=seq(0,2000,400), expand = c(0,0))+
  labs(x = "Sample",
       y = "Total Introgressed\nSize (Mb)",
       fill = "Clade") +
  theme_classic()+
  theme(axis.text.y = element_text(size=18, colour="black"),
        axis.text.x = element_text(size=18, colour="black", angle = 45, hjust = 1),
        axis.title = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=20), 
        legend.position = "bottom")

ggsave(filename=paste0("./plots/08_TotalPerSample_Win100kb_sim1MbWinStep200kb.png"),
       device = png, units="cm",
       width=12, height=12)
ggsave(filename=paste0("./plots/08_TotalPerSample_Win100kb_sim1MbWinStep200kb.pdf"),
       device = pdf, units="cm",
       width=12, height=12)




# proportions per Chr:

table_chr_size <- read_table("chr_size.txt", F)
names(table_chr_size) <- c("chrS", "size")
head(table_chr_size)

chrID="12"
# Define the list of chromosomes
chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

label_table2 <- data.frame(hapN=seq(1,48), singleL_list, AcName_list, simpleHapNum=rep(seq(1,4), 120)) %>%
  mutate(sample=paste0(singleL_list, "_hap", hapN)) %>%
  mutate(label=paste0(AcName_list, "_", simpleHapNum)) %>%
  select(sample, label, AcName_list)

chrID="12"

# Loop through each chromosome and process the data
combined_summary_table_Chr <- data.frame() ; for (chrID in chr_list) {
  
  # Define file paths for merged_filtered_table_C3
  file_path_C3 <- paste0("merged_filtered_table_C3_chr", chrID, ".csv")
  
  # Define file paths for merged_filtered_table_C4_S
  file_path_C4_S <- paste0("merged_filtered_table_C4_S_chr", chrID, ".csv")
  
  file_path_C4_S_sim <- paste0("merged_filtered_table_C4_S_C3_sim_chr", chrID, ".csv")
  
  # Load and restore factor levels for merged_filtered_table_C4_S
  merged_filtered_table_C4_S <- read.csv(file_path_C4_S)
  
  # Load and restore factor levels for merged_filtered_table_C3
  merged_filtered_table_C3 <- read.csv(file_path_C3)
  
  merged_filtered_table_C4_S_sim <- read.csv(file_path_C4_S_sim, stringsAsFactors = FALSE)
  
  # # Process and combine data
  # summary_table <- merged_filtered_table_C4_S %>%
  #   rbind(merged_filtered_table_C3) %>%
  #   mutate(sample = str_replace_all(cul_Sample, c("^Sample" = "", "_Hap_" = "_hap"))) %>%
  #   left_join(label_table2, by="sample") %>%
  #   group_by(chr, pos, cul_Sample, Clade, sample, label, AcName_list) %>%
  #   summarise(SupportingTest=n(), .groups = 'drop') %>%
  #   filter(SupportingTest>1) %>%
  #   group_by(chr, pos, cul_Sample, sample, label, AcName_list) %>%
  #   summarise(NClades=n(), CladeSingle=Clade[1], .groups = 'drop') %>%
  #   filter(NClades==1) %>%
  #   # group_by(chr, label) %>%
  #   group_by(chr, CladeSingle, label) %>%
  #   summarise(BlockSize=n()*100000, .groups = 'drop') %>%
  #   mutate(chrS=paste0("chr", chrID))
  # 
  winSizeInt <- 1000000
  min_block=0.3*(winSizeInt/100000)

  summary_table <-  merged_filtered_table_C4_S_sim %>%
    left_join(label_table2 %>%
                select(sample, AcName_list), by="sample") %>%
    filter(count>min_block) %>%
    mutate(pos=start) %>%
    group_by(chr, CladeSingle, label_or) %>%
    # group_by(chr, label_or) %>%
    summarise(BlockSize=n()*200000, .groups = 'drop') %>%
      mutate(chrS=paste0("chr", chrID))

  # Combine the summary table with the combined summary table
  combined_summary_table_Chr <- bind_rows(combined_summary_table_Chr, summary_table)
}

# combined_summary_table now contains the combined data from all chromosomes


combined_summary_table_Chr %>%
  mutate(CladeG=factor(CladeSingle, levels=c("C4_S", "C3"))) %>%
  left_join(table_chr_size, by="chrS") %>%
  mutate(frac_chr=BlockSize/size) %>%
  ggplot(aes(x = factor(chr), y = frac_chr, fill = CladeG)) +
  # ggplot(aes(x = factor(chr), y = BlockSize/1000000, fill = CladeG)) +
  # ggplot(aes(x = factor(chr), y = BlockSize)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#72b680","#7e73b6"))+
  scale_y_continuous(breaks=seq(0,1,0.1))+
  labs(x = "Chromosome",
       y = "Introgressed Fraction",
       fill = "Clade") +
  theme_classic()+
  theme(axis.text.y = element_text(size=18, colour="black"),
        axis.text.x = element_text(size=18, colour="black"),
        axis.title = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=20), 
        legend.position = "top")

ggsave(filename=paste0("./plots/09_MeanPerChr_Win100kb_sim1MbWinStep200kb.png"),
       device = png, units="cm",
       width=15, height=10)
ggsave(filename=paste0("./plots/09_MeanPerChr_Win100kb_sim1MbWinStep200kb.pdf"),
       device = pdf, units="cm",
       width=15, height=10)



combined_summary_table_Chr %>%
  mutate(CladeG=factor(CladeSingle, levels=c("C4_S", "C3"))) %>%
  left_join(table_chr_size, by="chrS") %>%
  mutate(frac_chr=BlockSize/size) %>%
  group_by(chr, label_or) %>%
  summarise(total_frac_chr=sum(frac_chr)) %>%
  ggplot(aes(x = factor(chr), y = total_frac_chr)) +
  # ggplot(aes(x = factor(chr), y = BlockSize/1000000, fill = CladeG)) +
  # ggplot(aes(x = factor(chr), y = BlockSize)) +
  geom_boxplot(color="black", fill="gray40") +
  # scale_fill_manual(values = c("#72b680","#7e73b6"))+
  scale_y_continuous(breaks=seq(0,1,0.1))+
  labs(x = "Chromosome",
       y = "Introgressed Fraction") +
  theme_classic()+
  theme(axis.text.y = element_text(size=18, colour="black"),
        axis.text.x = element_text(size=18, colour="black"),
        axis.title = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=20), 
        legend.position = "top")

ggsave(filename=paste0("./plots/09b_MeanPerChr_Win100kb_sim1MbWinStep200kb.png"),
       device = png, units="cm",
       width=15, height=10)
ggsave(filename=paste0("./plots/09b_MeanPerChr_Win100kb_sim1MbWinStep200kb.pdf"),
       device = pdf, units="cm",
       width=15, height=10)


combined_summary_table_Chr %>%
  mutate(CladeG=factor(CladeSingle, levels=c("C4_S", "C3"))) %>%
  left_join(table_chr_size, by="chrS") %>%
  mutate(frac_chr=BlockSize/size) %>%
  group_by(chr, label_or) %>%
  summarise(total_frac_chr=sum(frac_chr)) %>%
  ggplot(aes(x = factor(chr), y = total_frac_chr)) +
  # ggplot(aes(x = factor(chr), y = BlockSize/1000000, fill = CladeG)) +
  # ggplot(aes(x = factor(chr), y = BlockSize)) +
  geom_boxplot(color="black", fill="gray40") +
  # scale_fill_manual(values = c("#72b680","#7e73b6"))+
  scale_y_continuous(breaks=seq(0.1,1,0.2))+
  labs(x = "Chromosome",
       y = "Introgressed Fraction") +
  theme_classic()+
  theme(axis.text.y = element_text(size=18, colour="black"),
        axis.text.x = element_text(size=18, colour="black"),
        axis.title = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=20),
        legend.position = "top")+
  coord_flip()
ggsave(filename=paste0("./plots/09b_v2_MeanPerChr_Win100kb_sim1MbWinStep200kb.png"),
       device = png, units="cm",
       width=10, height=16)
ggsave(filename=paste0("./plots/09b_v2_MeanPerChr_Win100kb_sim1MbWinStep200kb.pdf"),
       device = pdf, units="cm",
       width=10, height=16)

getwd()

combined_summary_table_Chr %>%
  mutate(CladeG=factor(CladeSingle, levels=c("C4_S", "C3"))) %>%
  left_join(table_chr_size, by="chrS") %>%
  mutate(frac_chr=BlockSize/size) %>%
  group_by(chr, label_or) %>%
  summarise(total_frac_chr=sum(frac_chr)) %>%
  pull(total_frac_chr) %>%
  max()


combined_summary_table_Chr %>%
  pull(label_or) %>%
  unique() %>%
  length()
  head()

combined_summary_table_Chr %>%
  mutate(CladeG=factor(CladeSingle, levels=c("C4_S", "C3"))) %>%
  left_join(table_chr_size, by="chrS") %>%
  mutate(totalSize=sum(table_chr_size$size)*40) %>%
  mutate(frac_chr=BlockSize/totalSize) %>%
  ungroup() %>%
  pull(frac_chr) %>%
  sum()


combined_summary_table_Chr_frac <- combined_summary_table_Chr %>%
  mutate(CladeG=factor(CladeSingle, levels=c("C4_S", "C3"))) %>%
  left_join(table_chr_size, by="chrS") %>%
  mutate(frac_chr=BlockSize/size) %>%
  group_by(chr, label_or) %>%
  summarise(total_frac_chr=sum(frac_chr)) 


model <- lm(total_frac_chr ~ factor(chr), data=combined_summary_table_Chr_frac)
Anova(model, Type="III") ; summary(model)

model <- lmer(total_frac_chr ~ factor(chr) + (1 | label), data=combined_summary_table_Chr_frac)


# Block size distribution:

chrID="12"
# Define the list of chromosomes
chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

# Loop through each chromosome and process the data
introgression_blocks <- data.frame() ; for (chrID in chr_list) {
  
  file_path_C4_S_sim <- paste0("merged_filtered_table_C4_S_C3_sim_chr", chrID, ".csv")
  merged_filtered_table_C4_S_sim <- read.csv(file_path_C4_S_sim, stringsAsFactors = FALSE)
  
  winSizeInt <- 1000000
  min_block=0.3*(winSizeInt/100000)
  
  # Define the window size (200 kb)
  window_size <- 200000
  
  # Example usage with your data frame
  summary_table <- merged_filtered_table_C4_S_sim %>%
    filter(count > min_block) %>%
    mutate(window = start) %>%
    select(label_or, chr, window) %>%
    arrange(label_or, chr, window) %>%
    group_by(label_or, chr) %>%
    mutate(diff = c(0, diff(window))) %>%
    mutate(block = cumsum(diff != window_size)) %>%
    group_by(label_or, chr, block) %>%
    summarise(start_block = min(window),
              end_block = max(window) + window_size - 1, .groups = 'drop') %>%
    select(label_or, chr, start_block, end_block)
  
  # Combine the summary table with the combined introgression_blocks
  introgression_blocks <- bind_rows(introgression_blocks, summary_table)
}



introgression_blocks %>%
  mutate(blocksize=end_block-start_block+1) %>%
  ggplot(aes(x = blocksize/1000000)) +
  geom_histogram(binwidth = 1, color="black", fill="gray40")+
  scale_x_continuous(breaks=seq(0,60,5))+
  scale_y_continuous(breaks=c(1,10,100,1000, 4000), trans="log1p")+
  labs(x = "Introgressed block sizes (Mb)",
       y = "Count")+
  theme_classic()+
  theme(axis.text.y = element_text(size=18, colour="black"),
        axis.text.x = element_text(size=18, colour="black"),
        axis.title = element_text(size=20))

ggsave(filename=paste0("./plots/10_SizeDistribution_Win100kb_sim1MbWinStep200kb.png"),
       device = png, units="cm",
       width=12, height=10)
ggsave(filename=paste0("./plots/10_SizeDistribution_Win100kb_sim1MbWinStep200kb.pdf"),
       device = pdf, units="cm",
       width=12, height=10)




introgression_blocks %>%
  mutate(blocksize=end_block-start_block+1) %>%
  ggplot(aes(x = blocksize/1000000)) +
  geom_histogram(binwidth = 1, color="black", fill="gray40")+
  scale_x_continuous(breaks=seq(0,60,5))+
  scale_y_continuous(breaks=c(1,10,100, 500,1000, 2000, 4000), trans="log1p")+
  # scale_y_break(c(1000, 3000)) + 
  labs(x = "Block Size (Mb)",
       y = "Count")+
  facet_grid(chr ~ .)+
  theme_classic()+
  theme(axis.text.y = element_text(size=18, colour="black"),
        axis.text.x = element_text(size=18, colour="black"),
        axis.title = element_text(size=20), 
        strip.text = element_text(size=20))

ggsave(filename=paste0("./plots/11_SizeDistributionPerChr_Win100kb_sim1MbWinStep200kb.png"),
       device = png, units="cm",
       width=18, height=32)
ggsave(filename=paste0("./plots/11_SizeDistributionPerChr_Win100kb_sim1MbWinStep200kb.pdf"),
       device = pdf, units="cm",
       width=18, height=32)


introgression_blocks %>%
  mutate(blocksize=end_block-start_block+1) %>%
  pull(blocksize) %>%
  max()








# variation of introgression along the genome:

table_chr_size <- read_table("chr_size.txt", F)
names(table_chr_size) <- c("chrS", "size")
head(table_chr_size)

chrID="12"
# Define the list of chromosomes
chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

# label_table2 <- data.frame(hapN=seq(1,48), singleL_list, AcName_list, simpleHapNum=rep(seq(1,4), 120)) %>%
#   mutate(sample=paste0(singleL_list, "_hap", hapN)) %>%
#   mutate(label=paste0(AcName_list, "_", simpleHapNum)) %>%
#   select(sample, label, AcName_list)

chrID="10"

# Loop through each chromosome and process the data
combined_summary_table_NIntHap <- data.frame() ; for (chrID in chr_list) {
  
  # Define file paths for merged_filtered_table_C3
  file_path_C3 <- paste0("merged_filtered_table_C3_chr", chrID, ".csv")
  
  # Define file paths for merged_filtered_table_C4_S
  file_path_C4_S <- paste0("merged_filtered_table_C4_S_chr", chrID, ".csv")
  
  file_path_C4_S_sim <- paste0("merged_filtered_table_C4_S_C3_sim_chr", chrID, ".csv")
  
  # Load and restore factor levels for merged_filtered_table_C4_S
  merged_filtered_table_C4_S <- read.csv(file_path_C4_S)
  
  # Load and restore factor levels for merged_filtered_table_C3
  merged_filtered_table_C3 <- read.csv(file_path_C3)
  
  merged_filtered_table_C4_S_sim <- read.csv(file_path_C4_S_sim, stringsAsFactors = FALSE)
  
  # # # Process and combine data
  # summary_table <-merged_filtered_table_C4_S %>%
  #   rbind(merged_filtered_table_C3) %>%
  #   mutate(sample = str_replace_all(cul_Sample, c("^Sample" = "", "_Hap_" = "_hap"))) %>%
  #   left_join(label_table2, by="sample") %>%
  #   group_by(chr, pos, cul_Sample, Clade, sample, label, AcName_list) %>%
  #   summarise(SupportingTest=n(), .groups = 'drop') %>%
  #   filter(SupportingTest>1) %>%
  #   group_by(chr, pos, cul_Sample, sample, label, AcName_list) %>%
  #   summarise(NClades=n(), CladeSingle=Clade[1], .groups = 'drop') %>%
  #   filter(NClades==1) %>%
  #   group_by(chr, pos) %>%
  #   summarise(NintrogressedHap=n())
  # # 
  #   
  # summary_table<-data.frame(chr=as.integer(chrID),
  #            pos=seq(min(merged_filtered_table_C4_S$pos),
  #                    max(merged_filtered_table_C4_S$pos), 100000)) %>%
  #   left_join(summary_table, by=c("chr", "pos")) %>%
  #   mutate(NintHap=ifelse(is.na(NintrogressedHap), 0, NintrogressedHap)) %>%
  #   select(-NintrogressedHap)
  
  
  winSizeInt <- 1000000
  min_block=0.3*(winSizeInt/100000)
  #
  summary_table <- merged_filtered_table_C4_S_sim %>%
    filter(count>min_block) %>%
    mutate(pos=start) %>%
    group_by(label_or, chr, pos) %>%
    mutate(maxCount=max(count)) %>%
    filter(count==maxCount) %>%
    summarise(IntHap=1, CladeS=CladeSingle[1]) %>%
    group_by(chr, pos) %>%
    summarise(NintrogressedHap=n())
  
  summary_table<-data.frame(chr=as.integer(chrID),
                            pos=seq(min(merged_filtered_table_C4_S_sim$start),
                                    max(merged_filtered_table_C4_S_sim$start), 200000)) %>%
    left_join(summary_table, by=c("chr", "pos")) %>%
    mutate(NintHap=ifelse(is.na(NintrogressedHap), 0, NintrogressedHap)) %>%
    select(-NintrogressedHap)
  
  # Combine the summary table with the combined summary table
  combined_summary_table_NIntHap <- bind_rows(combined_summary_table_NIntHap, summary_table)
}

# combined_summary_table now contains the combined data from all chromosomes

winSize <- 1000000
combined_summary_table_NIntHap %>%
  mutate(win=floor(pos/winSize)*winSize) %>%
  group_by(chr, win) %>%
  summarise(MeanIntHap=mean(NintHap)) %>%
  # filter(MeanNHaplotypes<24) %>%
  ggplot(aes(win/1000000, MeanIntHap))+
  geom_line(alpha=0.8)+
  scale_x_continuous(breaks=seq(0,100,15)) +
  labs(x = "Pos. (Mb)",
       y = "N. Int. Haplo.") +
  facet_grid(. ~ chr, scale="free", space="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12), 
        axis.title = element_text(size=14))














# variation of introgression along the genome:

table_chr_size <- read_table("chr_size.txt", F)
names(table_chr_size) <- c("chrS", "size")
head(table_chr_size)

chrID="12"
# Define the list of chromosomes
chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")


# for i in range(0,len(allele_count)):
#   for j in range(0,len(allele_count)):
#   if j > i:
#   NumSS += allele_count[i]*allele_count[j]
# NCop = (Nsamples*(Nsamples-1))/2
# pi = NumSS/NCop

chrID="10"

# Function to calculate cumulative recombination rate between two positions
get_pi_IntHap <- function(CladeS) {
  NumSS = 0
  intHap = length(CladeS)
  nonintHap = 40 - intHap
  allHap = c(CladeS, rep("NoInt", nonintHap))
  for (i in seq(1, 40)){
    for (j in seq(1, 40)){
      if (j > i){
        NumSS = NumSS + ((allHap[i]==allHap[j])*1)
      }
    }
  }
  Nsamples=40
  NCop = (Nsamples*(Nsamples-1))/2
  pi = NumSS/NCop
  return(pi)
}

DM_ref_chrSizes <- data.frame(chr = c("chr01", "chr02", "chr03", "chr04", "chr05", "chr06", "chr07", "chr08", "chr09", "chr10", "chr11", "chr12"), 
                              totalLen = c(88591686, 46102915, 60707570, 69236331, 55599697, 59091578, 57639317, 59226000, 67600300, 61044151, 46777387, 59670755))

# Loop through each chromosome and process the data
combined_summary_table_PiIntHap <- data.frame() ; for (chrID in chr_list) {
  
  # Define file paths for merged_filtered_table_C3
  file_path_C3 <- paste0("merged_filtered_table_C3_chr", chrID, ".csv")
  
  # Define file paths for merged_filtered_table_C4_S
  file_path_C4_S <- paste0("merged_filtered_table_C4_S_chr", chrID, ".csv")
  
  file_path_C4_S_sim <- paste0("merged_filtered_table_C4_S_C3_sim_chr", chrID, ".csv")
  
  # Load and restore factor levels for merged_filtered_table_C4_S
  merged_filtered_table_C4_S <- read.csv(file_path_C4_S)
  
  # Load and restore factor levels for merged_filtered_table_C3
  merged_filtered_table_C3 <- read.csv(file_path_C3)
  
  merged_filtered_table_C4_S_sim <- read.csv(file_path_C4_S_sim, stringsAsFactors = FALSE)
  
  # # # Process and combine data
  # summary_table <- merged_filtered_table_C4_S %>%
  #   rbind(merged_filtered_table_C3) %>%
  #   mutate(sample = str_replace_all(cul_Sample, c("^Sample" = "", "_Hap_" = "_hap"))) %>%
  #   left_join(label_table2, by="sample") %>%
  #   group_by(chr, pos, cul_Sample, Clade, sample, label, AcName_list) %>%
  #   summarise(SupportingTest=n(), .groups = 'drop') %>%
  #   filter(SupportingTest>1) %>%
  #   group_by(chr, pos, cul_Sample, sample, label, AcName_list) %>%
  #   summarise(NClades=n(), CladeSingle=Clade[1], .groups = 'drop') %>%
  #   filter(NClades==1) %>%
  #   group_by(chr, pos) %>%
  #   summarise(DivInt=get_pi_IntHap(CladeSingle), NSamples=n())
  # 
  # summary_table<-data.frame(chr=as.integer(chrID),
  #            pos=seq(min(merged_filtered_table_C4_S$pos),
  #                    max(merged_filtered_table_C4_S$pos), 100000)) %>%
  #   left_join(summary_table, by=c("chr", "pos")) %>%
  #   mutate(DivIntAll=ifelse(is.na(DivInt), 0, DivInt)) %>%
  #   select(-DivInt)
  
  
  winSizeInt <- 1000000
  min_block=0.5*(winSizeInt/100000)
  #
  summary_table <-  merged_filtered_table_C4_S_sim %>%
    filter(count>min_block) %>%
    mutate(pos=start) %>%
    group_by(label_or, chr, pos) %>%
    mutate(maxCount=max(count)) %>%
    filter(count==maxCount) %>%
    summarise(IntHap=1, CladeS=CladeSingle[1]) %>%
    group_by(chr, pos) %>%
    summarise(DivInt=get_pi_IntHap(CladeS))
  
  max_chrVal <- DM_ref_chrSizes %>%
    filter(chr==paste0("chr", chrID)) %>%
    pull(totalLen)
  
  summary_table<-data.frame(chr=as.integer(chrID),
                            pos=seq(min(merged_filtered_table_C4_S_sim$start),
                                    max_chrVal, 200000)) %>%
    left_join(summary_table, by=c("chr", "pos")) %>%
    mutate(DivIntAll=ifelse(is.na(DivInt), 0, DivInt)) %>%
    select(-DivInt)
  
  # Combine the summary table with the combined summary table
  combined_summary_table_PiIntHap <- bind_rows(combined_summary_table_PiIntHap, summary_table)
}

centromere_table <- read.csv("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/04_LD/centromeres.csv", header = T) %>%
  separate(chr, into=c("other", "chrn"), sep = "r") %>%
  mutate(chr=as.integer(chrn)) %>%
  select(chr, start, end)

head(centromere_table)

# combined_summary_table now contains the combined data from all chromosomes
winSize <- 1000000
combined_summary_table_PiIntHap %>%
  arrange(chr, pos) %>%
  mutate(group = cumsum(c(TRUE, diff(DivIntAll != 0) != 0)),
         zero_count = ave(DivIntAll == 0, group, FUN = cumsum)) %>%
  group_by(group) %>%
  mutate(run_length = max(zero_count * (DivIntAll == 0))) %>%
  ungroup() %>%
  filter(!(DivIntAll == 0 & run_length <= 10)) %>%
  select(-group, -zero_count, -run_length) %>%
  mutate(win=floor(pos/winSize)*winSize) %>%
  group_by(chr, win) %>%
  summarise(MeanDivInt=mean(DivIntAll)) %>%
  # filter(MeanNHaplotypes<24) %>%
  ggplot(aes(win/1000000, MeanDivInt))+
  geom_ribbon(aes(ymin=ifelse(MeanDivInt<0.5, MeanDivInt, 0.5),
                  max=0.5), fill="darkorange") +
  geom_line(alpha=0.8)+
  # geom_point(alpha=0.3)+
  geom_line(data=centromere_table %>%
              gather(key="posB", value = "pos", -chr) %>%
              mutate(posy = -0.1),
            aes(pos/1000000,posy), colour="darkgreen", linewidth=1.5, alpha=0.5) +
  scale_x_continuous(breaks=seq(0,100,15), expand = c(0,0)) +
  labs(x = "Pos. (Mb)",
       y = "Diversity in\nHap. Ancestry") +
  facet_grid(. ~ chr, scale="free", space="free")+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "grey80"), 
        strip.text = element_text(size=12), 
        axis.text = element_text(size=12), 
        axis.title = element_text(size=14))

ggsave("./plots/12_VariationIntrogressionAloneGenome_1MbWin.png", height=2, width=15)
ggsave("./plots/12_VariationIntrogressionAloneGenome_1MbWin.pdf", height=2, width=15)

getwd()


combined_summary_table_PiIntHap %>%
  arrange(chr, pos) %>%
  mutate(group = cumsum(c(TRUE, diff(DivIntAll != 0) != 0)),
         zero_count = ave(DivIntAll == 0, group, FUN = cumsum)) %>%
  group_by(group) %>%
  mutate(run_length = max(zero_count * (DivIntAll == 0))) %>%
  ungroup() %>%
  filter(!(DivIntAll == 0 & run_length <= 2)) %>%
  select(-group, -zero_count, -run_length) %>%
  mutate(win=floor(pos/winSize)*winSize) %>%
  group_by(chr, win) %>%
  summarise(MeanDivInt=mean(DivIntAll)) %>%
  mutate(pos=win, DivHapAnc=MeanDivInt) %>%
  ungroup() %>%
  mutate(numeros_str = sprintf("%02d", as.numeric(as.character(chr))),
         chr = factor(paste0("chr", numeros_str))) %>%
  select(chr, pos, DivHapAnc) %>%
  write.csv(file = "diversity_hap_ancestry_1MbWin.csv", row.names = FALSE)




# correlation introgression diversity and genetic diversity (Pi):

combined_summary_table_PiIntHap %>%
  arrange(chr, pos) %>%
  mutate(group = cumsum(c(TRUE, diff(DivIntAll != 0) != 0)),
         zero_count = ave(DivIntAll == 0, group, FUN = cumsum)) %>%
  group_by(group) %>%
  mutate(run_length = max(zero_count * (DivIntAll == 0))) %>%
  ungroup() %>%
  filter(!(DivIntAll == 0 & run_length <= 2)) %>%
  select(-group, -zero_count, -run_length) %>%
  mutate(win=floor(pos/winSize)*winSize) %>%
  group_by(chr, win) %>%
  summarise(MeanDivInt=mean(DivIntAll)) %>%
  mutate(pos=win, DivHapAnc=MeanDivInt) %>%
  ungroup() %>%
  mutate(numeros_str = sprintf("%02d", as.numeric(as.character(chr))),
         chr = factor(paste0("chr", numeros_str))) %>%
  select(chr, pos, DivHapAnc) %>%
  head()


getwd()
"/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/05_ABBA_BABA/DSuite"




# Loading Pi values:
setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes")

library(ggplot2)
library(tidyverse)

# population parameters along the genome:
# chr     win     sum_mean_pi     sum_mean_pi2    SS      meanNsamples    mean_pi W_theta Taj_D

table_pop <- read_table("AllSamplesMultiAlleles_Win10000_ConsideringmissingData_AllCHR_popParameters.txt", F)
names(table_pop) <- c("chr", "win", "sum_pi", "SS", "NVar", "mean_pi_pw", "W_theta", "Taj_D", "Mean_NSamples")


winSize = 1000000

table_pop %>%
  mutate(pos=floor(win/winSize)*winSize) %>%
  group_by(chr, pos) %>%
  summarise(meanPi=mean(mean_pi_pw, na.rm = T)) %>%
  head()
  

# merding pi and ancestry values:
merged_combined_summary_table_PiIntHap_Pi <- combined_summary_table_PiIntHap %>%
  arrange(chr, pos) %>%
  mutate(group = cumsum(c(TRUE, diff(DivIntAll != 0) != 0)),
         zero_count = ave(DivIntAll == 0, group, FUN = cumsum)) %>%
  group_by(group) %>%
  mutate(run_length = max(zero_count * (DivIntAll == 0))) %>%
  ungroup() %>%
  filter(!(DivIntAll == 0 & run_length <= 2)) %>%
  select(-group, -zero_count, -run_length) %>%
  mutate(win=floor(pos/winSize)*winSize) %>%
  group_by(chr, win) %>%
  summarise(MeanDivInt=mean(DivIntAll)) %>%
  mutate(pos=win, DivHapAnc=MeanDivInt) %>%
  ungroup() %>%
  mutate(numeros_str = sprintf("%02d", as.numeric(as.character(chr))),
         chr = factor(paste0("chr", numeros_str))) %>%
  select(chr, pos, DivHapAnc) %>%
  left_join(table_pop %>%
              mutate(pos=floor(win/winSize)*winSize) %>%
              group_by(chr, pos) %>%
              summarise(meanPi=mean(mean_pi_pw, na.rm = T)), by=c("chr", "pos")) 



model <- lm(meanPi ~ DivHapAnc + factor(chr), data=merged_combined_summary_table_PiIntHap_Pi)
Anova(model, Type="III") ; summary(model)

  
correlation <- cor(merged_combined_summary_table_PiIntHap_Pi$DivHapAnc, merged_combined_summary_table_PiIntHap_Pi$meanPi)
print(paste("Correlation between DivHapAnc and meanPi: ", correlation))


```



#### Whole genome

Concatenate vcf files

script: mergeChr.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00
#SBATCH -J MergeVCF

vcf=$(ls ../chr/*vcf.gz)
bcftools concat --threads 5 $vcf | grep -v "##contig=<ID=scaffold" | awk -F"\t" 'BEGIN{ FS=OFS="\t" }{if($1 ~ /^chr/) $3="." ; print $0}' | gzip -c > WholeGenome_filterQUAL30.vcf.gz


```



Job:

```bash

cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/merged_vcf/wholegenome

mamba activate env_others
sbatch mergeChr.sh
# 07.06.2024: 2575389

```



Script: ABBA_BABA_WG.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=24:00:00

sample1=$1
sample2=$2
sample3=$3
sample4=$4
output_prefix=$5
folder=$6

mkdir -p $folder

vcf=/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/merged_vcf/wholegenome/WholeGenome_filterQUAL30.vcf.gz

bcftools query -l $vcf | awk '{print $1"\txxx"}' | \
    sed "s/${sample1}.*xx/${sample1}\t${sample1}/g" | \
    sed "s/${sample2}.*xx/${sample2}\t${sample2}/g" | \
    sed "s/${sample3}.*xx/${sample3}\t${sample3}/g" | \
    sed "s/${sample4}.*xx/${sample4}\tOutgroup/g" > species_sets_WholeGenome_${output_prefix}.txt

echo -e "P1\tP2\tP3\tDstatistic\tZ-score\tp-value\tf4-ratio\tBBAA\tABBA\tBABA" > Values_WholeGenome_$output_prefix.txt

# delete files from previous window
if test -f ./species_sets_WholeGenome_${output_prefix}_WholeGenome_$output_prefix"_tree.txt" ; then
  rm species_sets_WholeGenome_${output_prefix}_WholeGenome_$output_prefix"_"*.txt
fi

# calculate D values
/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/Dsuite/Build/Dsuite Dtrios -c -n WholeGenome_$output_prefix -t astral.merged_trees.WGenome.min25.nwk $vcf species_sets_WholeGenome_${output_prefix}.txt

# add values to output file:
grep "" species_sets_WholeGenome_${output_prefix}_WholeGenome_$output_prefix"_tree.txt" | grep -v "ABBA" >> Values_WholeGenome_$output_prefix.txt

rm species_sets_WholeGenome_${output_prefix}.txt

mv Values_WholeGenome_$output_prefix.txt $folder"/"

```

Job:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/05_ABBA_BABA/DSuite_WG

mamba activate env_others


vcf=/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/merged_vcf/wholegenome/WholeGenome_filterQUAL30.vcf.gz

bcftools query -l $vcf | grep Sample > list_Cul_samples.txt

for cul in $( cat list_Cul_samples.txt ); do echo $cul
shortCul=$(echo $cul | sed 's/Sample//g' | sed 's/_Hap_//g')
sample2=SRR15458993
sample3=SRR15458986
sample4=SRR15458989
folder=bre_ver_Cul_cho
output_prefix=$folder"_"$shortCul
sbatch ./ABBA_BABA_WG.sh $cul $sample2 $sample3 $sample4 $output_prefix $folder
done
# 07.06.2024: 2575421..2575460

for cul in $( cat list_Cul_samples.txt ); do echo $cul
shortCul=$(echo $cul | sed 's/Sample//g' | sed 's/_Hap_//g')
sample2=SRR15458993
sample3=SRR15458982
sample4=SRR15458989
folder=bre_bol_Cul_cho
output_prefix=$folder"_"$shortCul
sbatch ./ABBA_BABA_WG.sh $cul $sample2 $sample3 $sample4 $output_prefix $folder
done
# 07.06.2024: 2575461..2575500


for cul in $( cat list_Cul_samples.txt ); do echo $cul
shortCul=$(echo $cul | sed 's/Sample//g' | sed 's/_Hap_//g')
sample2=SRR15458993
sample3=SRR15458984
sample4=SRR15458989
folder=bre_com_Cul_cho
output_prefix=$folder"_"$shortCul
sbatch ./ABBA_BABA_WG.sh $cul $sample2 $sample3 $sample4 $output_prefix $folder
done
# 07.06.2024: 2575501..2575540

for cul in $( cat list_Cul_samples.txt ); do echo $cul
shortCul=$(echo $cul | sed 's/Sample//g' | sed 's/_Hap_//g')
sample2=SRR15458984
sample3=SRR15458982
sample4=SRR15458989
folder=com_bol_Cul_cho
output_prefix=$folder"_"$shortCul
sbatch ./ABBA_BABA_WG.sh $cul $sample2 $sample3 $sample4 $output_prefix $folder
done
# 07.06.2024: 2575541..2575580



for cul in $( cat list_Cul_samples.txt ); do echo $cul
shortCul=$(echo $cul | sed 's/Sample//g' | sed 's/_Hap_//g')
sample2=SRR15458993
sample3=SRR15458989
sample4=SRR15458990
folder=Cul_bre_choma_more
output_prefix=$folder"_"$shortCul
sbatch ./ABBA_BABA_WG.sh $cul $sample2 $sample3 $sample4 $output_prefix $folder
done
# 07.06.2024



for cul in $( cat list_Cul_samples.txt ); do echo $cul
shortCul=$(echo $cul | sed 's/Sample//g' | sed 's/_Hap_//g')
sample2=SRR15458993
sample3=SRR15458980
sample4=SRR15458990
folder=Cul_bre_caj_more
output_prefix=$folder"_"$shortCul
sbatch ./ABBA_BABA_WG.sh $cul $sample2 $sample3 $sample4 $output_prefix $folder
done
# 07.06.2024



for cul in $( cat list_Cul_samples.txt ); do echo $cul
shortCul=$(echo $cul | sed 's/Sample//g' | sed 's/_Hap_//g')
sample2=SRR15458993
sample3=SRR15458959
sample4=SRR15458990
folder=Cul_bre_and_more
output_prefix=$folder"_"$shortCul
sbatch ./ABBA_BABA_WG.sh $cul $sample2 $sample3 $sample4 $output_prefix $folder
done
# 07.06.2024

```





## Pairwise distance to wild species



Script: CalculateDis.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00


chr=$1

vcf=/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/merged_vcf/chr/chr$chr""_filterQUAL30.vcf.gz

#vcf=/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/splitVCF/01/split.000200.vcf.gz

# this is martin version but it is quite slow...
python3 /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/genomics_general/VCF_processing/parseVCF.py -i $vcf -o chr$chr""_filterQUAL30.geno.gz

python /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/genomics_general/distMat.py --windType coordinate -w 100000 -m 200  \
--addWindowID \
--windowDataOutFile chr$chr""_filterQUAL30.winData \
-f phased \
--ploidy 2 \
-g chr$chr""_filterQUAL30.geno.gz \
-o chr$chr""_filterQUAL30.dis \
-T 5
# 


```



Job:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/08_distanceMatrix_perWin


mamba activate mypython3


for chr in $(seq -w 1 12)
do echo $chr
sbatch ./CalculateDis.sh $chr 
done
# 2553496..2553507


```



script: ExtractLowestDis.py

```python
import numpy as np
import sys

def filter_distance_matrix(header, matrix, cultivar_samples, wild_samples):
    """
    Filters the distance matrix to include only cultivar samples in rows and wild samples in columns.

    Parameters:
    header (list): A list of sample names corresponding to the rows and columns of the matrix.
    matrix (numpy.ndarray): A 2D numpy array representing the distance matrix.
    cultivar_samples (list): A list of sample names corresponding to cultivar samples.
    wild_samples (list): A list of sample names corresponding to wild samples.

    Returns:
    tuple: A tuple containing the list of cultivar samples in the correct order, the list of wild samples in the correct order, and the filtered distance matrix.
    """
    
    # Get the indices of the cultivar samples
    cultivar_indices = [header.index(sample) for sample in cultivar_samples if sample in header]
    ordered_cultivar_samples = [header[i] for i in cultivar_indices]
    
    # Get the indices of the wild samples
    wild_indices = [header.index(sample) for sample in wild_samples if sample in header]
    ordered_wild_samples = [header[i] for i in wild_indices]
    
    # Filter the matrix
    filtered_matrix = matrix[np.ix_(cultivar_indices, wild_indices)]
    
    return ordered_cultivar_samples, ordered_wild_samples, filtered_matrix



def find_lowest_distances(distance_matrix, cultivar_samples, wild_samples):
    """
    Finds the lowest distances between cultivar samples and wild samples in the distance matrix.
    Parameters:
    distance_matrix (numpy.ndarray): The distance matrix.
    cultivar_samples (list): List of cultivar samples.
    wild_samples (list): List of wild samples.
    Returns:
    list: A list of tuples containing the lowest distances, cultivar samples, and wild samples.
    """
    lowest_distances = []
    # Iterate over each row in the distance matrix
    for i, cultivar_sample in enumerate(cultivar_samples):
        # Find the indices of the lowest distances in the current row
        sorted_indices = np.argsort(distance_matrix[i])  # Sort indices based on distances
        lowest_indices = sorted_indices[:3]  # Get indices of lowest 3 distance
        # Store the lowest distances and corresponding samples
        lowest_dist_values = [distance_matrix[i, idx] for idx in lowest_indices]
        max_dist = max(lowest_dist_values)
        NTopSamples = sum(distance_matrix[i]<=max_dist)
        # Handle multiple samples with the same top value
        if NTopSamples < 10:
            for idx, dist in enumerate(distance_matrix[i]):
                if dist <= max_dist:
                    lowest_distances.append([cultivar_sample, wild_samples[idx], dist])
        else:
            lowest_distances.append([cultivar_sample, "NA", "NA"])
    return lowest_distances



def read_samples(file_path):
    """
    Reads cultivar samples from a text file and returns a list of these samples.
    Parameters:
    file_path (str): The path to the text file containing the cultivar samples.
    Returns:
    list: A list of cultivar samples.
    """
    with open(file_path, 'r') as file:
        # Read all lines from the file and strip any leading/trailing whitespace characters
        list_samples = [line.strip() for line in file.readlines()]
    return list_samples




# Processing file:

filename="chr01_filterQUAL30.dis"
win_information="ed_chr01_filterQUAL30.winData"
cultivar_samples_file="list_cultivars.txt"
wild_samples_file="list_wild.txt"
output_file_name="test.txt"


filename=str(sys.argv[1])
win_information=str(sys.argv[2])
cultivar_samples_file=str(sys.argv[3])
wild_samples_file=str(sys.argv[4])
output_file_name=str(sys.argv[5])

cultivar_samples = read_samples(cultivar_samples_file)
wild_samples = read_samples(wild_samples_file)
Win_info = read_samples(win_information)
output_file = open(output_file_name, "w")



win=0
for line in open(filename, "r"):
    parts = line.strip().split()
    if len(parts) > 1:
        header.append(parts[0])
        matrix.append([float(x) for x in parts[1:]])
    elif win==0:
        header = []
        matrix = []
        win+=1
    else:
        ordered_cultivar_samples, ordered_wild_samples, filtered_matrix = filter_distance_matrix(header, np.array(matrix), cultivar_samples, wild_samples)
        lowest_distances_array = find_lowest_distances(filtered_matrix, ordered_cultivar_samples, ordered_wild_samples)
        for topvalue in lowest_distances_array:
            if topvalue[2] != "NA":
                values="\t".join(map(str, topvalue))
                output_file.write(f'{Win_info[win-1]}\t{values}\n')
        win+=1
        header = []
        matrix = []

output_file.close()


```



Job:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/08_distanceMatrix_perWin

mamba activate mypython3

head -62 chr01_filterQUAL30.dis | grep -v "^60" | awk '{print $1}' | grep Sample >  list_cultivars.txt
head -62 chr01_filterQUAL30.dis | grep -v "^60" | awk '{print $1}' | grep SRR >  list_wild.txt

# windowID,scaffold,start,end,mid,sites
for chr in $(seq -w 1 12)
do echo $chr
sed 's/sites,/sites\n/g' chr$chr""_filterQUAL30.winData | grep -v window > ed_chr$chr""_filterQUAL30.winData
done


for chr in $(seq -w 1 12)
do echo $chr
python3 ExtractLowestDis.py chr$chr""_filterQUAL30.dis ed_chr$chr""_filterQUAL30.winData list_cultivars.txt list_wild.txt top_chr$chr""_filterQUAL30.txt
done

for chr in $(seq -w 1 12)
do echo $chr
python3 ExtractLowestDis3.py chr$chr""_filterQUAL30.dis ed_chr$chr""_filterQUAL30.winData list_cultivars.txt list_wild.txt top3_chr$chr""_filterQUAL30.txt
done

```

#### Plot:



```R

setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/06_mergedHaplotypes_PanGenome")

# Whole chromosomes with haplotype colours, but clustering haplotypes:
winSize="200kb"
stepSize="50kb"
minSS="50SS"
# chrID="12"
MaxGap="06"
n=10

singleL <- c("A","B","C","D","E","F","G","H","I","J")
AcName <- c("WhR", "BdF", "EgH", "PrW", "Flo", "Ack", "Fla", "Kat", "Lum", "EdB")

singleL_list <- c()
for (i in singleL){
  singleL_list <- c(singleL_list, rep(i, 48))
}

AcName_list <- c()
for (i in AcName){
  AcName_list <- c(AcName_list, rep(i, 48))
}

label_table <- data.frame(hapN=seq(1,48), singleL_list, AcName_list, simpleHapNum=rep(seq(1,4), 120)) %>%
  mutate(sample=paste0(singleL_list, "_hap", hapN)) %>%
  mutate(label=paste0(AcName_list, "_", simpleHapNum)) %>%
  select(sample, label)


singleL <- c("A","B","C","D","E","F","G","H","I","J")
AcName <- c("WhR", "BdF", "EgH", "PrW", "Flo", "Ack", "Fla", "Kat", "Lum", "EdB")

sample <- c()
for (i in singleL){
  for (hapN in seq(1,4)){
    sample <- c(sample, paste0(i, "_hap", hapN))
  }
}

labelS <- c()
for (i in AcName){
  for (hapN in seq(1,4)){
    labelS <- c(labelS, paste0(i, "_", hapN))
  }
}

label_tableD <- data.frame(sample, labelS)


singleL <- c("A","B","C","D","E","F","G","H","I","J")
Cultivar_samples <- c()
for (i in singleL){
  for (hapN in seq(1,4)){
    Cultivar_samples <- c(Cultivar_samples, paste0("Sample", i, "_Hap_", hapN))
  }
}


wildSamplesInfo <- read_tsv("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/08_distanceMatrix_perWin/WildSample_info_sim.txt", T)
wildSamplesInfo$Clade <- factor(wildSamplesInfo$Clade, levels=c("C4_N", "C4_S", "C4_S_C3", "C3_C4N", "C3", "C1_2"))
wildSamplesInfo

Cultivar_samples

head(label_table)

chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

# chrID="10"

for (chrID in chr_list){
    
    # first we load the information of clustering on haplotypes to have the order of samples:
  
  setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/06_mergedHaplotypes_PanGenome")
  
  #Read in haplotypeNumberFile, with haplotype numbers reduced by the python script. 
  hapNums.in <- read.table(paste0("clusteringLargeHaplotypes_", winSize, "_step", stepSize,
                                  "_haplotypeGroups_min", minSS, "_chr", chrID, ".txt"), 
                           header=T, sep="\t")
  
  hapNums.in2 <- hapNums.in %>%
    left_join(label_table, by="sample") 
  
  matrix_hapNums.in <- hapNums.in2 %>%
    filter(Mean_NAN_haplo<0.8) %>%
    select(pos, label, haplotype) %>%
    spread(key=pos, value=haplotype, fill = NA) %>%
    select(-label) %>%
    as.matrix()
  
  samples_ID <- hapNums.in2 %>%
    filter(Mean_NAN_haplo<0.8) %>%
    select(pos, label, haplotype) %>%
    spread(key=pos, value=haplotype, fill = NA) %>%
    select(label) %>% unlist() %>% as.vector()
  
  d <- dist(matrix_hapNums.in, method = "euclidean")
  H.fit <- hclust(d, method="ward.D")
  H.fit$labels<-samples_ID
  
  # this line sets the order of samples:
  hapNums.in2$label <- factor(hapNums.in2$label, levels = H.fit$labels[H.fit$order])
  
  # Sequence similarity blocks:
  
  setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/08_distanceMatrix_perWin")
  
  values_WSp <- read_tsv(paste0("top3_chr", chrID, "_filterQUAL30.txt"), F)
  names(values_WSp) <- c("windowID","chr","win","end","mid","sites", "CulSample", "WildSample", "dis")
  
  # change Sample labels:
  
  values_WSp %>%
    mutate(sample = str_replace_all(CulSample, c("^Sample" = "", "_Hap_" = "_hap"))) %>%
    left_join(label_tableD, by="sample") %>%
    mutate(label_or=factor(labelS, levels = H.fit$labels[H.fit$order])) %>%
    select(-labelS) %>%
    left_join(wildSamplesInfo, by="WildSample") %>%
    filter(dis<0.15) %>%
    group_by(chr, win, CulSample, label_or, Clade) %>%
    summarise(NSamples=n()) %>%
    group_by(chr, win, CulSample, label_or) %>%
    mutate(MaxNSamples=max(NSamples)) %>%
    filter(NSamples==MaxNSamples) %>%
    group_by(chr, win, CulSample, label_or) %>%
    summarise(NClades=length(unique(Clade)), 
              CladeG=Clade[1]) %>%
    filter(NClades==1) %>%
    filter(CladeG!="C4_N") %>%
    ggplot(aes(x = win/1000000, y = label_or, fill = factor(CladeG))) +
    geom_tile() +
    scale_fill_manual(values = c("#72b680","#7e73b6","#ae4f3b"))+
    scale_x_continuous(breaks=seq(0,100,5), expand = c(0,0), limits = c(0,90))+
    labs(title = paste0("Sequence similarity: Chr",chrID),
         x = "Position (Mb)",
         y = "Haplotype",
         fill = "Clade")+
    theme_classic()+
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=18),
          axis.line.x = element_blank(), 
          axis.title.x = element_text(hjust = 0.3))
  
  ggsave(filename=paste0("plots/01_SeqSmilarity_perClade_chr", chrID, "_WinSize100kb_maxDis015.png"), 
         device = png, units="cm", 
         width=30, height=18)
  ggsave(filename=paste0("plots/01_SeqSmilarity_perClade_chr", chrID, "_WinSize100kb_maxDis015.pdf"), 
         device = pdf, units="cm", 
         width=30, height=18)
}


```



### Calculate distance from VCF Cultivar to Wild samples

script: divergence_haplotype_wildSamples.py

```python
import sys
import numpy as np
from cyvcf2 import VCF
from collections import defaultdict

# the script produces a table with SampleCultivarID, SampleWIldID, N variance sites considered, number of matching variants, proportion of variants found (presence absence), proportion of variants found but considering ploidy in wild samples, wild clade, wild species.

# mamba activate mypython3
#genotypes_file = str("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/splitVCF/01/split.000846.vcf.gz")
#output = str("test")

genotypes_file = str(sys.argv[1])
chr = str(sys.argv[2])
output = str(sys.argv[3])

print(f'Processing file: {genotypes_file}')

# processing Genotype table:
vf = VCF(genotypes_file, strict_gt=True)
list_samples = vf.samples
cult_haplotypes_samples = list_samples[20:]
wild_samples = list_samples[:20]

Ncult_haplotypes_samples = len(cult_haplotypes_samples)
Nwild_samples = len(wild_samples)

Nlines = np.zeros(Ncult_haplotypes_samples)
allele_suppost = np.zeros((Ncult_haplotypes_samples, Nwild_samples))
allele_suppost_fq = np.zeros((Ncult_haplotypes_samples, Nwild_samples))

#variant=next(vf)
ini_pos=0
for variant in vf:
    if ini_pos==0:
        ini_pos=variant.POS
    list_gt_total = [item[:-1] for item in variant.genotypes]
    list_gt_cult_haplotypes = list_gt_total[20:]
    list_gt_wild = list_gt_total[:20]
    for index, cul_sample in enumerate(list_gt_cult_haplotypes):
        if cul_sample !=[-1,-1]:
            WSample_SameHaplotype1 = np.array([(sum([(i==cul_sample[0])*1 for i in WSample])>0)*1 for WSample in list_gt_wild])
            WSample_SameHaplotype2 = np.array([(sum([(i==cul_sample[1])*1 for i in WSample])>0)*1 for WSample in list_gt_wild])
            WSample_SameHaplotype = (WSample_SameHaplotype1+WSample_SameHaplotype2)/2
            WSample_SameHaplotypeFq1 = np.array([sum([(i==cul_sample[0])*1 for i in WSample]) for WSample in list_gt_wild])
            WSample_SameHaplotypeFq2 = np.array([sum([(i==cul_sample[1])*1 for i in WSample]) for WSample in list_gt_wild])
            WSample_SameHaplotypeFq = (WSample_SameHaplotypeFq1+WSample_SameHaplotypeFq2)/4
            allele_suppost[index,] += WSample_SameHaplotype
            allele_suppost_fq[index,] += WSample_SameHaplotypeFq
            Nlines[index] += 1

suppost_perWsample_output = open(f'{output}.txt', "w")
suppost_perWsample_output_top = open(f'{output}_top90quantileSamples.txt', "w")

for index, haplotype in enumerate(cult_haplotypes_samples):
  if Nlines[index] != 0:
    pro_allele_suppost = allele_suppost[index,]/Nlines[index]
    pro_allele_suppost_fq = allele_suppost_fq[index,]/(Nlines[index])
    min_Suppost=np.quantile(pro_allele_suppost_fq, 0.90)
    for Windex, WSample in enumerate(wild_samples):
      suppost_perWsample_output_join = '\t'.join(map(str, [chr, ini_pos, haplotype, WSample, Nlines[index], allele_suppost[index,Windex], pro_allele_suppost[Windex], pro_allele_suppost_fq[Windex]]))
      suppost_perWsample_output.write(f'{suppost_perWsample_output_join}\n')
      if pro_allele_suppost_fq[Windex]>=min_Suppost:
        suppost_perWsample_output_top.write(f'{suppost_perWsample_output_join}\n')


suppost_perWsample_output.close()
suppost_perWsample_output_top.close()


```

script: divergence_haplotype_wildSamples.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00

chr=$1
output_prefix=$2

mkdir -p $output_prefix
mkdir -p $output_prefix/chr$chr

ls -1 /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/splitVCF/$chr"/"*.vcf.gz > list_vcf_Chr$chr.txt

for vcf in $( cat list_vcf_Chr$chr.txt ); do echo $vcf
# Extract the filename
filename="${vcf##*/}"
# Extract the number using parameter expansion
number="${filename#split.}"
window="${number%%.*}"
python3 divergence_haplotype_wildSamples.py $vcf $chr $output_prefix/chr$chr"/"dis_chr$chr"_"Win$window
done

```

Job:

```bash

cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/08_distanceMatrix_perWin

mamba activate mypython3

for chr in $(seq -w 1 12)
do echo $chr
sbatch ./divergence_haplotype_wildSamples.sh $chr distance_values
done
# 28.05.2024:2554624..2554635

for chr in $(seq -w 1 12)
do echo $chr
cat ./distance_values/chr$chr"/"dis_*top90quantileSamples.txt > ./distance_values/top90quantileSamples_dis_chr$chr.txt
done
# 


```

#### Plot:



```R

setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/06_mergedHaplotypes_PanGenome")

# Whole chromosomes with haplotype colours, but clustering haplotypes:
winSize="200kb"
stepSize="50kb"
minSS="50SS"
# chrID="12"
MaxGap="06"
n=10

singleL <- c("A","B","C","D","E","F","G","H","I","J")
AcName <- c("WhR", "BdF", "EgH", "PrW", "Flo", "Ack", "Fla", "Kat", "Lum", "EdB")

singleL_list <- c()
for (i in singleL){
  singleL_list <- c(singleL_list, rep(i, 48))
}

AcName_list <- c()
for (i in AcName){
  AcName_list <- c(AcName_list, rep(i, 48))
}

label_table <- data.frame(hapN=seq(1,48), singleL_list, AcName_list, simpleHapNum=rep(seq(1,4), 120)) %>%
  mutate(sample=paste0(singleL_list, "_hap", hapN)) %>%
  mutate(label=paste0(AcName_list, "_", simpleHapNum)) %>%
  select(sample, label)


singleL <- c("A","B","C","D","E","F","G","H","I","J")
AcName <- c("WhR", "BdF", "EgH", "PrW", "Flo", "Ack", "Fla", "Kat", "Lum", "EdB")

sample <- c()
for (i in singleL){
  for (hapN in seq(1,4)){
    sample <- c(sample, paste0(i, "_hap", hapN))
  }
}

labelS <- c()
for (i in AcName){
  for (hapN in seq(1,4)){
    labelS <- c(labelS, paste0(i, "_", hapN))
  }
}

label_tableD <- data.frame(sample, labelS)


singleL <- c("A","B","C","D","E","F","G","H","I","J")
Cultivar_samples <- c()
for (i in singleL){
  for (hapN in seq(1,4)){
    Cultivar_samples <- c(Cultivar_samples, paste0("Sample", i, "_Hap_", hapN))
  }
}


wildSamplesInfo <- read_tsv("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/08_distanceMatrix_perWin/WildSample_info_sim.txt", T)
wildSamplesInfo$Clade <- factor(wildSamplesInfo$Clade, levels=c("C4_N", "C4_S", "C4_S_C3", "C3_C4N", "C3", "C1_2"))
wildSamplesInfo

Cultivar_samples

head(label_table)

chr_list <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

chrID="12"

for (chrID in chr_list){
  
  # first we load the information of clustering on haplotypes to have the order of samples:
  
  setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/06_mergedHaplotypes_PanGenome")
  
  #Read in haplotypeNumberFile, with haplotype numbers reduced by the python script. 
  hapNums.in <- read.table(paste0("clusteringLargeHaplotypes_", winSize, "_step", stepSize,
                                  "_haplotypeGroups_min", minSS, "_chr", chrID, ".txt"), 
                           header=T, sep="\t")
  
  hapNums.in2 <- hapNums.in %>%
    left_join(label_table, by="sample") 
  
  matrix_hapNums.in <- hapNums.in2 %>%
    filter(Mean_NAN_haplo<0.8) %>%
    select(pos, label, haplotype) %>%
    spread(key=pos, value=haplotype, fill = NA) %>%
    select(-label) %>%
    as.matrix()
  
  samples_ID <- hapNums.in2 %>%
    filter(Mean_NAN_haplo<0.8) %>%
    select(pos, label, haplotype) %>%
    spread(key=pos, value=haplotype, fill = NA) %>%
    select(label) %>% unlist() %>% as.vector()
  
  d <- dist(matrix_hapNums.in, method = "euclidean")
  H.fit <- hclust(d, method="ward.D")
  H.fit$labels<-samples_ID
  
  # this line sets the order of samples:
  hapNums.in2$label <- factor(hapNums.in2$label, levels = H.fit$labels[H.fit$order])
  
  # Sequence similarity blocks:
  
  setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/08_distanceMatrix_perWin")
  
  values_WSp <- read_tsv(paste0("./distance_values/top90quantileSamples_dis_chr", chrID, ".txt"), F)
  names(values_WSp) <- c("chr", "win", "CulSample", "WildSample", "NTotalVar", "NVar", "support", "supportfq")
  # names(values_WSp) <- c("windowID","chr","win","end","mid","sites", "CulSample", "WildSample", "dis")
  
  # change Sample labels:
  
  # values_WSp %>%
  #   mutate(sample = str_replace_all(CulSample, c("^Sample" = "", "_Hap_" = "_hap"))) %>%
  #   left_join(label_tableD, by="sample") %>%
  #   mutate(label_or=factor(labelS, levels = H.fit$labels[H.fit$order])) %>%
  #   select(-labelS) %>%
  #   left_join(wildSamplesInfo, by="WildSample") %>%
  #   filter(supportfq>0.80) %>%
  #   group_by(chr, win, CulSample, label_or, Clade) %>%
  #   summarise(NSamples=n()) %>%
  #   group_by(chr, win, CulSample, label_or) %>%
  #   mutate(MaxNSamples=max(NSamples)) %>%
  #   filter(NSamples==MaxNSamples) %>%
  #   group_by(chr, win, CulSample, label_or) %>%
  #   summarise(NClades=length(unique(Clade)), 
  #             CladeG=Clade[1]) %>%
  #   filter(NClades==1) %>%
  #   # filter(CladeG!="C4_N") %>%
  #   ggplot(aes(x = win/1000000, y = label_or, fill = factor(CladeG))) +
  #   geom_tile() +
  #   # scale_fill_manual(values = c("#72b680","#7e73b6","#ae4f3b"))+
  #   scale_x_continuous(breaks=seq(0,100,5), expand = c(0,0), limits = c(0,90))+
  #   labs(title = paste0("Sequence similarity: Chr",chrID),
  #        x = "Position (Mb)",
  #        y = "Haplotype",
  #        fill = "Clade")+
  #   theme_classic()+
  #   theme(axis.text = element_text(size=12),
  #         axis.title = element_text(size=18),
  #         axis.line.x = element_blank(), 
  #         axis.title.x = element_text(hjust = 0.3))
  
  # ggsave(filename=paste0("plots/02_fromVCF_SeqSmilarity_perClade_chr", chrID, "_WinSize100kb_maxDis085.png"), 
  #        device = png, units="cm", 
  #        width=30, height=18)
  # ggsave(filename=paste0("plots/02_fromVCF_SeqSmilarity_perClade_chr", chrID, "_WinSize100kb_maxDis085.pdf"), 
  #        device = pdf, units="cm", 
  #        width=30, height=18)
  
  NSp_Cldes_table <- values_WSp %>%
    mutate(sample = str_replace_all(CulSample, c("^Sample" = "", "_Hap_" = "_hap"))) %>%
    left_join(label_tableD, by="sample") %>%
    mutate(label_or=factor(labelS, levels = H.fit$labels[H.fit$order])) %>%
    select(-labelS) %>%
    left_join(wildSamplesInfo, by="WildSample") %>%
    filter(supportfq>0.80) %>%
    group_by(chr, win, CulSample, label_or) %>%
    summarise(NSp=length(unique(Species)), NClades=length(unique(Clade))) 
  
  values_WSp %>%
    mutate(sample = str_replace_all(CulSample, c("^Sample" = "", "_Hap_" = "_hap"))) %>%
    left_join(label_tableD, by="sample") %>%
    mutate(label_or=factor(labelS, levels = H.fit$labels[H.fit$order])) %>%
    select(-labelS) %>%
    left_join(wildSamplesInfo, by="WildSample") %>%
    filter(supportfq>0.80) %>%
    left_join(NSp_Cldes_table, by=c("chr", "win", "CulSample", "label_or")) %>%
    filter(NSp<5) %>%
    filter(NClades<2) %>%
    group_by(chr, win, CulSample, label_or, Clade) %>%
    summarise(NSamples=n()) %>%
    group_by(chr, win, CulSample, label_or) %>%
    mutate(MaxNSamples=max(NSamples)) %>%
    filter(NSamples==MaxNSamples) %>%
    group_by(chr, win, CulSample, label_or) %>%
    summarise(NClades=length(unique(Clade)), 
              CladeG=Clade[1]) %>%
    filter(NClades==1) %>%
    # filter(CladeG!="C4_N") %>%
    ggplot(aes(x = win/1000000, y = label_or, fill = factor(CladeG))) +
    geom_tile() +
    scale_fill_manual(values = c("#92e7a4", "#72b680","#7e73b6","#ae4f3b"))+
    scale_x_continuous(breaks=seq(0,100,10), expand = c(0,0), limits = c(0,90))+
    labs(title = paste0("Sequence similarity: Chr",chrID),
         x = "Position (Mb)",
         y = "Haplotype",
         fill = "Clade")+
    theme_classic()+
    theme(axis.text = element_text(size=18, colour="black"),
          axis.title = element_text(size=20),
          title = element_text(size=20),
          legend.text = element_text(size=18),
          axis.line.x = element_blank(), 
          axis.title.x = element_text(hjust = 0.3))
  
  ggsave(filename=paste0("plots/02b_fromVCF_SeqSmilarity_perClade_chr", chrID, "_WinSize100kb_maxDis085.png"), 
         device = png, units="cm", 
         width=25, height=23)
  ggsave(filename=paste0("plots/02b_fromVCF_SeqSmilarity_perClade_chr", chrID, "_WinSize100kb_maxDis085.pdf"), 
         device = pdf, units="cm", 
         width=25, height=23)
  
  
  # Function to calculate the number of lines per clade in overlapping windows
  calculate_windows <- function(df, window_size = 1000000, step_size = 200000) {
    
    NsubWin = floor((window_size/100000)*min_block)
    
    result <- list()
    
    for (chr in unique(df$chr)) {
      chr_df <- df %>% filter(chr == chr)
      max_pos <- max(chr_df$win)
      
      for (start in seq(0, max_pos, by = step_size)) {
        end <- start + window_size
        
        window_df <- chr_df %>%
          filter(win >= start & win < end) %>%
          group_by(CulSample, label_or, CladeG) %>%
          summarize(count = n(), .groups = 'drop') %>%
          group_by(CulSample, label_or) %>%
          mutate(total_count = sum(count), 
                 max_count=(max(count)), 
                 NClades=length(unique(CladeG)))
        
        if (nrow(window_df) > 0) {
          window_df <- window_df %>%
            mutate(chr = chr, start = start, end = end)
          
          result <- append(result, list(window_df))
        }
      }
    }
    
    final_result <- bind_rows(result)
    return(final_result)
  }
  
  winSizeInt=1000000
  
  data_file_path <- paste0("values_WSp_WinSlides_chr", chrID, "_WinSize100kb_maxDis085_sim1MbWin.csv")
  levels_file_path <- paste0("factor_levels_chr", chrID, ".rds")
  
  # Function to save the data frame and factor levels
  save_data_with_levels <- function(data, data_path, levels_path) {
    # Save the data frame
    write.csv(data, file = data_path, row.names = FALSE)
    
    # Save the levels of factors
    factor_levels <- lapply(data, function(column) {
      if (is.factor(column)) {
        levels(column)
      } else {
        NULL
      }
    })
    saveRDS(factor_levels, file = levels_path)
  }
  
  # Check if the data file exists
  if (file.exists(data_file_path)) {
    # If the file exists, read it into a data frame
    values_WSp_WinSlides <- read.csv(data_file_path, stringsAsFactors = FALSE)
    
    # Load the factor levels
    factor_levels <- readRDS(levels_file_path)
    
    # Restore the factor levels
    for (col_name in names(factor_levels)) {
      if (!is.null(factor_levels[[col_name]])) {
        values_WSp_WinSlides[[col_name]] <- factor(values_WSp_WinSlides[[col_name]], levels = factor_levels[[col_name]])
      }
    }
  } else {
    # If the file does not exist, create the data frame
    values_WSp_WinSlides <- values_WSp %>%
      mutate(sample = str_replace_all(CulSample, c("^Sample" = "", "_Hap_" = "_hap"))) %>%
      left_join(label_tableD, by = "sample") %>%
      mutate(label_or = factor(labelS, levels = H.fit$labels[H.fit$order])) %>%
      select(-labelS) %>%
      left_join(wildSamplesInfo, by = "WildSample") %>%
      filter(supportfq > 0.80) %>%
      left_join(NSp_Cldes_table, by = c("chr", "win", "CulSample", "label_or")) %>%
      filter(NSp < 5) %>%
      filter(NClades < 2) %>%
      group_by(chr, win, CulSample, label_or, Clade) %>%
      summarise(NSamples = n()) %>%
      group_by(chr, win, CulSample, label_or) %>%
      mutate(MaxNSamples = max(NSamples)) %>%
      filter(NSamples == MaxNSamples) %>%
      group_by(chr, win, CulSample, label_or) %>%
      summarise(NClades = length(unique(Clade)), 
                CladeG = Clade[1]) %>%
      filter(NClades == 1) %>%
      calculate_windows(winSizeInt, 200000)
    
    # Save the data frame and factor levels
    save_data_with_levels(values_WSp_WinSlides, data_file_path, levels_file_path)
  }
  
  values_WSp_WinSlides %>%
    filter(total_count>3) %>%
    mutate(FqClade=count/total_count) %>%
    filter(FqClade>0.6) %>%
    ggplot(aes(x = start/1000000, y = label_or, fill = factor(CladeG), colour = factor(CladeG))) +
    geom_tile() +
    scale_fill_manual(values = c("#92e7a4", "#72b680","#7e73b6","#ae4f3b"))+
    scale_colour_manual(values = c("#92e7a4", "#72b680","#7e73b6","#ae4f3b"))+
    scale_x_continuous(breaks=seq(0,100,10), expand = c(0,0), limits = c(0,90))+
    labs(title = paste0("Sequence similarity: Chr",chrID),
         x = "Position (Mb)",
         y = "Haplotype",
         fill = "Clade")+
    theme_classic()+
    theme(axis.text = element_text(size=18, colour="black"),
          axis.title = element_text(size=20),
          title = element_text(size=20),
          legend.text = element_text(size=18),
          axis.line.x = element_blank(), 
          axis.title.x = element_text(hjust = 0.3))
  
  ggsave(filename=paste0("plots/03_fromVCF_SeqSmilarity_perClade_chr", chrID, "_WinSize100kb_maxDis085_sim1MbWin.png"), 
         device = png, units="cm", 
         width=25, height=23)
  ggsave(filename=paste0("plots/03_fromVCF_SeqSmilarity_perClade_chr", chrID, "_WinSize100kb_maxDis085_sim1MbWin.pdf"), 
         device = pdf, units="cm", 
         width=25, height=23)
  
}


```



## Heterozygosity per sample



Script: Heterozygosity.py

```python
import sys
import gzip
from collections import defaultdict

def calculate_heterozygosity(vcf_file, output_file, exclude_samples):
    # Parse the exclude samples list
    exclude_samples_set = set(exclude_samples.split(','))

    # Initialize data structures
    windows = defaultdict(lambda: defaultdict(list))
    samples = []
    window_size = 1_000_000  # Using underscore for better readability

    # Read the VCF file
    with gzip.open(vcf_file, 'rt') if vcf_file.endswith('.gz') else open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#'):
                headers = line.strip().split('\t')
                samples = headers[9:]
                continue
            
            cols = line.strip().split('\t')
            chrom = cols[0]
            pos = int(cols[1])
            genotypes = cols[9:]
            
            window_start = (pos // window_size) * window_size
            
            for i, genotype in enumerate(genotypes):
                sample = samples[i]
                if sample in exclude_samples_set:
                    continue
                gt = genotype.split(':')[0]
                if (gt == '0/1') or (gt == '1/0'):
                    windows[(chrom, window_start)][sample].append(1)
                elif gt == './.':
                    windows[(chrom, window_start)][sample].append(None)
                else:
                    windows[(chrom, window_start)][sample].append(0)

    # Calculate heterozygosity and write to output file
    with open(output_file, 'w') as out:
        out.write('Sample,Chromosome,Window_Start,Mean_Heterozygosity\n')
        for (chrom, window_start), sample_data in windows.items():
            for sample, genotypes in sample_data.items():
                total_variants = len(genotypes)
                missing_data = genotypes.count(None)
                heterozygotes = sum(1 for gt in genotypes if gt == 1)
                if total_variants - missing_data > 0:
                    mean_heterozygosity = heterozygotes / (total_variants - missing_data)
                else:
                    mean_heterozygosity = 0
                out.write(f'{sample},{chrom},{window_start},{mean_heterozygosity}\n')

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python Heterozygosity.py <input_vcf> <output_csv> <exclude_samples>")
        print("Exclude samples should be a comma-separated list of sample names to exclude")
        sys.exit(1)
    
    input_vcf = sys.argv[1]
    output_csv = sys.argv[2]
    exclude_samples = sys.argv[3]
    calculate_heterozygosity(input_vcf, output_csv, exclude_samples)

```

bash script: Heterozygosity.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=12:00:00
#SBATCH -J He

vcf=$1
output=$2
samplesExcluded=$3

python Heterozygosity.py $vcf $output $samplesExcluded


```

job:

```bash

cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/11_He_perSample

mamba activate mypython3

vcf=/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/merged_vcf/chr/chr10_filterQUAL30.vcf.gz

python Heterozygosity.py $vcf chr10_filterQUAL30.csv

for chr in $(seq -w 1 12); do echo $chr
vcf=/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/merged_vcf/chr/chr$chr"_"filterQUAL30.vcf.gz

sbatch Heterozygosity.sh $vcf chr$chr"_"filterQUAL30.csv SampleA_Hap_1,SampleA_Hap_2,SampleA_Hap_3,SampleA_Hap_4,SampleB_Hap_1,SampleB_Hap_2,SampleB_Hap_3,SampleB_Hap_4,SampleC_Hap_1,SampleC_Hap_2,SampleC_Hap_3,SampleC_Hap_4,SampleD_Hap_1,SampleD_Hap_2,SampleD_Hap_3,SampleD_Hap_4,SampleE_Hap_1,SampleE_Hap_2,SampleE_Hap_3,SampleE_Hap_4,SampleF_Hap_1,SampleF_Hap_2,SampleF_Hap_3,SampleF_Hap_4,SampleG_Hap_1,SampleG_Hap_2,SampleG_Hap_3,SampleG_Hap_4,SampleH_Hap_1,SampleH_Hap_2,SampleH_Hap_3,SampleH_Hap_4,SampleI_Hap_1,SampleI_Hap_2,SampleI_Hap_3,SampleI_Hap_4,SampleJ_Hap_1,SampleJ_Hap_2,SampleJ_Hap_3,SampleJ_Hap_4 
done
# 30.06.2024: 2611007..2611018



```

ArtificialHeterozygosity.py

```python
import sys
import gzip
import random
from collections import defaultdict

def calculate_artificial_heterozygosity(vcf_file, output_file, target_samples):
    # Parse the target samples list
    target_samples_list = target_samples.split(',')

    # Initialize data structures
    windows = defaultdict(lambda: defaultdict(list))
    window_size = 1_000_000  # Using underscore for better readability

    # Create 20 random pairs from the target samples
    random.shuffle(target_samples_list)
    paired_samples = [target_samples_list[i:i + 2] for i in range(0, len(target_samples_list), 2)]

    if len(paired_samples) != 20:
        print("Error: Number of pairs is not equal to 20. Check the input sample list.")
        sys.exit(1)

    # Read the VCF file
    with gzip.open(vcf_file, 'rt') if vcf_file.endswith('.gz') else open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#'):
                headers = line.strip().split('\t')
                samples = headers[9:]
                sample_indices = {sample: idx for idx, sample in enumerate(samples) if sample in target_samples_list}
                continue
            
            cols = line.strip().split('\t')
            chrom = cols[0]
            pos = int(cols[1])
            genotypes = cols[9:]
            
            window_start = (pos // window_size) * window_size

            # Create artificial diploids for each pair
            for pair_idx, (sample1, sample2) in enumerate(paired_samples):
                idx1 = sample_indices[sample1]
                idx2 = sample_indices[sample2]
                
                gt1 = genotypes[idx1].split(':')[0]
                gt2 = genotypes[idx2].split(':')[0]

                if gt1 == './.' or gt2 == './.':
                    artificial_gt = None
                else:
                    alleles1 = gt1.replace('|', '/').split('/')
                    alleles2 = gt2.replace('|', '/').split('/')
                    allele1 = random.choice(alleles1)
                    allele2 = random.choice(alleles2)
                    artificial_gt = f"{allele1}/{allele2}"

                if artificial_gt == '0/1' or artificial_gt == '1/0':
                    windows[(chrom, window_start)][f"Artificial_{pair_idx+1}"].append(1)
                elif artificial_gt is None:
                    windows[(chrom, window_start)][f"Artificial_{pair_idx+1}"].append(None)
                else:
                    windows[(chrom, window_start)][f"Artificial_{pair_idx+1}"].append(0)

    # Calculate heterozygosity and write to output file
    with open(output_file, 'w') as out:
        out.write('Sample,Chromosome,Window_Start,Mean_Heterozygosity\n')
        for (chrom, window_start), sample_data in windows.items():
            for sample, genotypes in sample_data.items():
                total_variants = len(genotypes)
                missing_data = genotypes.count(None)
                heterozygotes = sum(1 for gt in genotypes if gt == 1)
                if total_variants - missing_data > 0:
                    mean_heterozygosity = heterozygotes / (total_variants - missing_data)
                else:
                    mean_heterozygosity = 0
                out.write(f'{sample},{chrom},{window_start},{mean_heterozygosity}\n')

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python Heterozygosity.py <input_vcf> <output_csv> <target_samples>")
        print("Target samples should be a comma-separated list of sample names to include")
        sys.exit(1)
    
    input_vcf = sys.argv[1]
    output_csv = sys.argv[2]
    target_samples = sys.argv[3]
    calculate_artificial_heterozygosity(input_vcf, output_csv, target_samples)

```

ArtificialHeterozygosity.sh

```bash
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=12:00:00
#SBATCH -J HeAr

vcf=$1
output=$2
samples=$3

python ArtificialHeterozygosity.py $vcf $output $samples


```



job:

```bash
cd /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/11_He_perSample

mamba activate mypython3

vcf=/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/merged_vcf/chr/chr10_filterQUAL30.vcf.gz


for chr in $(seq -w 1 12); do echo $chr
vcf=/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/merged_vcf/chr/chr$chr"_"filterQUAL30.vcf.gz

sbatch ArtificialHeterozygosity.sh $vcf CulArt_chr$chr"_"filterQUAL30.csv SampleA_Hap_1,SampleA_Hap_2,SampleA_Hap_3,SampleA_Hap_4,SampleB_Hap_1,SampleB_Hap_2,SampleB_Hap_3,SampleB_Hap_4,SampleC_Hap_1,SampleC_Hap_2,SampleC_Hap_3,SampleC_Hap_4,SampleD_Hap_1,SampleD_Hap_2,SampleD_Hap_3,SampleD_Hap_4,SampleE_Hap_1,SampleE_Hap_2,SampleE_Hap_3,SampleE_Hap_4,SampleF_Hap_1,SampleF_Hap_2,SampleF_Hap_3,SampleF_Hap_4,SampleG_Hap_1,SampleG_Hap_2,SampleG_Hap_3,SampleG_Hap_4,SampleH_Hap_1,SampleH_Hap_2,SampleH_Hap_3,SampleH_Hap_4,SampleI_Hap_1,SampleI_Hap_2,SampleI_Hap_3,SampleI_Hap_4,SampleJ_Hap_1,SampleJ_Hap_2,SampleJ_Hap_3,SampleJ_Hap_4 
done
# 30.06.2024: 2611019..2611030




```

Plots:

```R
setwd("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/11_He_perSample")

# Load necessary libraries
library(ggplot2)
library(dplyr)

# List of chromosome files for wild samples
wild_chrom_files <- list.files(pattern = "^chr[0-9]+_filterQUAL30.csv")
wild_chrom_files

# Read and combine all files into a single data frame for wild samples
wild_data <- lapply(wild_chrom_files, read.csv) %>% bind_rows()

# Add group information for wild samples
wild_data <- wild_data %>%
  mutate(Group = "Wild")

# Calculate the mean heterozygosity per sample across the whole genome for wild samples
mean_heterozygosity_per_sample_wild <- wild_data %>%
  group_by(Sample) %>%
  summarize(mean_heterozygosity = mean(Mean_Heterozygosity, na.rm = TRUE))

# Normalize heterozygosity values for wild samples
wild_data <- wild_data %>%
  left_join(mean_heterozygosity_per_sample_wild, by = "Sample") %>%
  mutate(Normalized_Heterozygosity = Mean_Heterozygosity / mean_heterozygosity)

# List of chromosome files for cultivar samples
cultivar_chrom_files <- list.files(pattern = "CulArt_chr[0-9]+_filterQUAL30.csv")

# Read and combine all files into a single data frame for cultivar samples
cultivar_data <- lapply(cultivar_chrom_files, read.csv) %>% bind_rows()

# Add group information for cultivar samples
cultivar_data <- cultivar_data %>%
  mutate(Group = "Cultivar")

# Calculate the mean heterozygosity per artificial sample across the whole genome for cultivar samples
mean_heterozygosity_per_sample_cultivar <- cultivar_data %>%
  group_by(Sample) %>%
  summarize(mean_heterozygosity = mean(Mean_Heterozygosity, na.rm = TRUE))

# Normalize heterozygosity values for cultivar samples
cultivar_data <- cultivar_data %>%
  left_join(mean_heterozygosity_per_sample_cultivar, by = "Sample") %>%
  mutate(Normalized_Heterozygosity = Mean_Heterozygosity / mean_heterozygosity)

# Combine wild and cultivar data into one data frame
all_data <- bind_rows(wild_data, cultivar_data)

# Extract and sort unique chromosome names from the data
chromosome_levels <- all_data %>%
  pull(Chromosome) %>%
  unique() %>%
  sort()

# Ensure consistent formatting of chromosome names
all_data <- all_data %>%
  mutate(Chromosome = factor(Chromosome, levels = chromosome_levels)) %>%
  arrange(Chromosome, Window_Start)

head(all_data)
# Create a plot

all_data %>%
  #filter(Normalized_Heterozygosity<4) %>%
  ggplot(aes(x = Window_Start, y = Normalized_Heterozygosity, group = Sample)) +
  geom_line(alpha=0.1) +
  facet_grid(Group ~ Chromosome, scales = "free_x", space = "free") +
  # labs(title = "Normalized Mean Heterozygosity along the Genome",
  #      x = "Genomic Position (Window Start)",
  #      y = "Normalized Mean Heterozygosity") +
  theme_minimal() +
  theme(legend.position = "none")

all_data %>%
  #filter((Chromosome!="chr10") | (Chromosome=="chr10" & Group=="Cultivar") | (Chromosome=="chr10" & (Window_Start/1000000)<60)) %>%
  group_by(Group, Chromosome, Window_Start) %>%
  summarise(median_Nor_heterozygosity=median(Normalized_Heterozygosity)) %>%
  #filter(median_Nor_heterozygosity<3) %>%
  ggplot(aes(x = Window_Start/1000000, y = median_Nor_heterozygosity, colour=Group)) +
  geom_line() +
  facet_grid(. ~ Chromosome, scale = "free", space = "free") +
  scale_x_continuous(breaks=seq(0,100,15))+
  labs(title = "Normilized mean heterozygosity",
       x = "Pos (MB)",
       y = "Mean He",
       Colour = "")+
  theme_classic()+
  # theme_minimal()+
  theme(axis.text = element_text(size=12, colour="black"),
        axis.title = element_text(size=14, colour="black"),
        title = element_text(size=20),
        strip.text = element_text(size=12), 
        panel.grid.major = element_line(colour="gray90"),
        # panel.grid = element_line(colour="gray50"),
        legend.position = "bottom", 
        legend.text = element_text(size=14), 
        legend.title = element_blank())


ggsave("01_meanHe_mean1Mb.png", height=3.5, width=15)
ggsave("01_meanHe_mean1Mb.pdf", height=3.5, width=15)


```




