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





