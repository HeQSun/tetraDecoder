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




