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



