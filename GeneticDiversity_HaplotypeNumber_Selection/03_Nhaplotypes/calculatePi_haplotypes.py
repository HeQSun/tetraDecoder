import os
import sys
from collections import Counter


file = 'haplotypeGroups_min10SS_chr01.txt'
output = 'test'

file = str(sys.argv[1])
output = str(sys.argv[2])
chrID = str(sys.argv[3])

Nsamples = 44
NCop = (Nsamples*(Nsamples-1))/2
current_pos = '0'
haplotypes = []

output_pi = open(f'{output}.txt', "w")

for line in open(file, "r"):
    pos, sample, hap  = line.strip().split("\t")
    if pos == current_pos:
        haplotypes.append(hap)
    else:
        # Calculate pi:
        NumSS = 0
        allele_count = list(Counter(haplotypes).values())
        for h_index1 in range(len(allele_count)):
            for h_index2 in range(len(allele_count)):
                if h_index2 > h_index1:
                    NumSS += allele_count[h_index1]*allele_count[h_index2]
        pi = NumSS/NCop
        output_pi.write(f'{chrID}\t{pos}\t{pi}\n')
        # re set dic:
        haplotypes = []
        haplotypes.append(hap)
        current_pos = pos

NumSS = 0
allele_count = list(Counter(haplotypes).values())
for h_index1 in range(len(allele_count)):
    for h_index2 in range(len(allele_count)):
        NumSS += allele_count[h_index1]*allele_count[h_index2]
pi = NumSS/NCop
output_pi.write(f'{chrID}\t{pos}\t{pi}\n')

output_pi.close()

