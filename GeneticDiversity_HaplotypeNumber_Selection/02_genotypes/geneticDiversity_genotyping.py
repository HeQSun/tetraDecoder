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
	pi2 = sum(genotypes) * (sampleN-sum(genotypes))
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



