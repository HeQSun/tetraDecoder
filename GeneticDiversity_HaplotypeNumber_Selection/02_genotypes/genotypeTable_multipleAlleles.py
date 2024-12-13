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


