import numpy as np
from collections import defaultdict
import math
import sys
import os

#file="AllSamplesPopPar_Win10000_chr01_missing_genotypeTable.txt"
#chr_used="chr01"
#Nsamples=int(48)
#winSize=int(10000)
#prefix_output="test"

file = str(sys.argv[1])
chr_used = str(sys.argv[2])
Nsamples=int(sys.argv[3])
winSize=int(sys.argv[4])
output=str(sys.argv[5])


def round_up_to_nearest_win(num, win_size):
	return math.ceil(int(num) / win_size) * win_size

# output file 1: table qieh allele fq:
output_AF = open(f'{output}_AF.txt', "w")
AF_dis = defaultdict(int)

#output file 2: population parameters per Window
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
total_e1_pop1 = c1_pop1/a1_pop1
total_e2_pop1 = c2_pop1/((a1_pop1**2)+a2_pop1)


for line in open(file, "r"):
	#chr, start, end, ref, alt, haplotypes  = line.strip().split("\t")
	if line[0:3] == "CHR":
		output_popParameters.write(f'chr\twin\tsum_mean_pi\tsum_mean_pi2\tSS\tmeanNsamples\tmean_pi\tW_theta\tTaj_D\n')
	else:
		line_values  = line.strip().split("\t")
		chr, start, end, ref, alt = line_values[0:5]
		win = round_up_to_nearest_win(start, winSize)
		#genotypes = np.array(list(map(int, line_values[5:])))
		genotypes = [int(x) for x in line_values[5:] if x!='nan']
		Nsamples = len(genotypes)
		p = sum(genotypes)/Nsamples
		q = 1- p
		# calculate allele frequency. values are stored in a dictionary (AF_dis)
		AF = round((p)*100)
		AF_dis[AF/100] += 1
		# calculate Pi per base.
		if Nsamples > 1: 
			pi = (Nsamples/(Nsamples-1))*2*p*q
			pi2 = sum(genotypes) * (Nsamples-sum(genotypes))
			pi2_unco = pi2/(((Nsamples-1)*Nsamples)/2)
			total_sum_harmonicMean = 0
			total_sum_harmonicMean2 = 0
			for sampleN in range(1,Nsamples):
				total_sum_harmonicMean += 1 / sampleN
				total_sum_harmonicMean2 += 1 / (sampleN**2)
			a1_pop1 = float(total_sum_harmonicMean)
			# add value per window:
			if win!=current_win:
				if current_win == 0:
					current_win = win
					mean_pi = pi
					mean_pi2 = pi2_unco
					#total_NVar = 1
					if (sum(genotypes) != Nsamples and sum(genotypes) != 0 ):
						NVar = 1
						meanNsamples = Nsamples
						W_theta = 1 / total_sum_harmonicMean
						W_theta_unco = 1 / a1_pop1
				else:
					# Watterson Theta:
					total_W_theta = W_theta/winSize
					# Tajima's D:
					Taj_D = (mean_pi2 - W_theta_unco)/(math.sqrt((total_e1_pop1*(NVar))+(total_e2_pop1*(NVar)*(NVar-1))))
					# produce output table
					output_popParameters.write(f'{chr}\t{current_win}\t{mean_pi}\t{mean_pi2}\t{NVar}\t{meanNsamples/NVar}\t{mean_pi/winSize}\t{total_W_theta}\t{Taj_D}\n')
					current_win = win
					mean_pi = pi
					mean_pi2 = pi2_unco
					#total_NVar = 1
					if (sum(genotypes) != Nsamples and sum(genotypes) != 0 ):
						meanNsamples = Nsamples
						NVar = 1
						W_theta = 1 / total_sum_harmonicMean
						W_theta_unco = 1 / a1_pop1
			else:
				mean_pi += pi
				mean_pi2 += pi2_unco
				#total_NVar += 1
				if (sum(genotypes) != Nsamples):
					meanNsamples += Nsamples
					NVar += 1
					W_theta += 1 / total_sum_harmonicMean
					W_theta_unco += 1 / a1_pop1

# produce allele fq table output:
for fq, totalVar in AF_dis.items():
	output_AF.write(f'{fq}\t{totalVar}\n')

output_AF.close()
output_popParameters.close()


