import os
import sys
import numpy as np
import pandas as pd

import scipy.cluster.hierarchy as shc
from scipy.spatial.distance import pdist
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





