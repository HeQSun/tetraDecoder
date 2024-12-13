
import numpy as np
#from collections import defaultdict
#import math
import sys
import os
#from collections import Counter

#file="AllSamplesPopPar_Win10000_chr01_missing_genotypeTable.txt"
file="test_AllSamplesMultiAlleles_chr01_missing_genotypeTable.txt"
#chr_used="chr01"
#Nsamples=int(48)
#winSize=int(10000)
output="LDtest2"
min_dis = int(500)
start_current = int(0)
end_pos = int(100000000)

file = str(sys.argv[1])
output=str(sys.argv[2])
min_dis = int(sys.argv[3])
start_current = int(sys.argv[4])
end_pos = int(sys.argv[5])



def most_frequent(List):
    unique, counts = np.unique(List, return_counts=True)
    index = np.argmax(counts)
    return unique[index]

#output file 1: LD output:
output_popParameters = open(f'{output}.NoSingletons.txt', "w")

for line in open(file, "r"):
	if line[0:3] == "CHR":
		#output_popParameters.write(f'chr\tstart\tstart2\ttotal\tp1\tq1\tp1q1\tD\tDp\tr2\n')
		output_popParameters.write(f'chr\tstart\tstart2\ttotal\tp1\tq1\tp1q1\tD\tr2\n')
	else:
		line_values  = line.strip().split("\t")
		chr, start, end, ref, alt = line_values[0:5]
		if int(start) > end_pos:
			break
		if int(start) - start_current > min_dis:
			genotypes = [x for x in line_values[5:]]
			Non_nan_genotypes = [x != 'nan' for x in genotypes] 
			index_Non_nan_genotypes = [x for x in range(len(Non_nan_genotypes)) if Non_nan_genotypes[x]] 
			if len(index_Non_nan_genotypes) > 10:
				# check start pos is polimorphic:
				tem_filterd_genotypes = [genotypes[i] for i in index_Non_nan_genotypes] 
				tem_allele_1_genotypes = most_frequent(tem_filterd_genotypes)
				tem_total_allele1 = sum([(x == tem_allele_1_genotypes)*1 for x in tem_filterd_genotypes])
				# this position is only used if it is polimorphic after filtering missing data:
				if tem_total_allele1 < (len(tem_filterd_genotypes)-1):
					start_current = int(start)
					start2_current = 0
					for line2 in open(file, "r"):
						if line2[0:3] == "CHR":
							continue
						else:
							line_values  = line2.strip().split("\t")
							chr2, start2, end2, ref2, alt2 = line_values[0:5]
							if (int(start2) > int(start)) and (int(start2) - start2_current > min_dis):
								genotypes2 = [x for x in line_values[5:]]
								Non_nan_genotypes2 = [x != 'nan' for x in genotypes2] 
								index_Non_nan_genotypes2 = [x for x in range(len(Non_nan_genotypes2)) if Non_nan_genotypes2[x]] 
								total_non_na = list(set(index_Non_nan_genotypes + index_Non_nan_genotypes2))
								if len(total_non_na) > 10:
									filterd_genotypes = [genotypes[i] for i in total_non_na] 
									filterd_genotypes2 = [genotypes2[i] for i in total_non_na] 
									allele_1_genotypes = most_frequent(filterd_genotypes)
									allele_1_genotypes2 = most_frequent(filterd_genotypes2)
									alleles_genotypes = [(x == allele_1_genotypes)*1 for x in filterd_genotypes]
									alleles_genotypes2 = [(x == allele_1_genotypes2)*1 for x in filterd_genotypes2]
									total = len(alleles_genotypes)
									Np1=sum(alleles_genotypes)
									Nq1=sum(alleles_genotypes2)
									p1=Np1/total
									q1=Nq1/total
									# final check that both varint sites are polimorphic and no singletons:
									if (Np1 < (total-1)) and (Nq1 < (total-1)):
										start2_current = int(start2)
										p1q1 = sum([(sum(i) == 2)*1 for i in zip(alleles_genotypes, alleles_genotypes2)])/total
										D=p1q1-(p1*q1)
										#if D<0:
										#	Dp = D/(max([(-p1*q1),(-(1-p1)*(1-q1))]))
										#else:
										#	Dp =  D/(min([(p1*(1-q1)),((1-p1)*(q1))])) 
										r2=(D**2)/(p1*(1-p1)*q1*(1-q1))
										#output_popParameters.write(f'{chr}\t{start}\t{start2}\t{total}\t{p1}\t{q1}\t{p1q1}\t{D}\t{Dp}\t{r2}\n')
										output_popParameters.write(f'{chr}\t{start}\t{start2}\t{total}\t{p1}\t{q1}\t{p1q1}\t{D}\t{r2}\n')


output_popParameters.close()





