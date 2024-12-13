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




