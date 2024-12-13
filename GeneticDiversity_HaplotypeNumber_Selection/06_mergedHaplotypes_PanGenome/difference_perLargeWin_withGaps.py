import os
import sys
import numpy as np
import pandas as pd

import scipy.cluster.hierarchy as shc
#from scipy.spatial.distance import pdist
from scipy.spatial.distance import cdist

import math

input = 'haplotypeGroups_min50SS_10KbWin_chr01.txt'
winSize=int(1000000)
step_size = int(50000)
output="test"
chr="chr01"

input = str(sys.argv[1])
winSize = int(sys.argv[2])
step_size = int(sys.argv[3])
output = str(sys.argv[4])
chr = str(sys.argv[5])

table = pd.read_table(input, sep="\t", names=['pos', 'sample', 'haplotype', 'missingP'])  
tableSpread = table[['pos', 'sample', 'haplotype']].pivot(index = 'pos', columns = "sample", values = "haplotype").reset_index().rename_axis("", axis = 1)
tableSpread_GAPS = table[['pos', 'sample', 'missingP']].pivot(index = 'pos', columns = "sample", values = "missingP").reset_index().rename_axis("", axis = 1)

output_file = open(f'{output}.txt', "w")
output_file.write(f'chr\tpos\thap1\thap2\tN10kbWin\tSS\n')

end = max(tableSpread['pos'])
list_haplotypes = list(tableSpread.columns)[1:]
NHap = len(list_haplotypes)
#list_gapHap = list(table[(table['missingP'] >= 0.9)]['haplotype'].unique())


def clean_dis(df1, df2):
	df1na = np.isnan(df1)
	df1clean = df1[~df1na]
	df2clean = df2[~df1na]
	df2na = np.isnan(df2clean)
	df1clean = df1clean[~df2na]
	df2clean = df2clean[~df2na]
	distance = (df1clean != df2clean).sum()
	return distance

for start in range(0,end,step_size):
	tableSpreadWin = tableSpread[(tableSpread['pos'] >= start) & (tableSpread['pos'] < start+winSize)]
	tableSpreadWin_GAPS = tableSpread_GAPS[(tableSpread_GAPS['pos'] >= start) & (tableSpread_GAPS['pos'] < start+winSize)]
	#gapsTable = np.isin(tableSpreadWin[list(tableSpreadWin.columns)[1:]].to_numpy(), list_gapHap)
	haplotypeTable = tableSpreadWin[list(tableSpreadWin.columns)[1:]].to_numpy().astype(float)
	haplotypeTable_GAPS = (tableSpreadWin_GAPS[list(tableSpreadWin_GAPS.columns)[1:]].to_numpy()) > 0.5
	haplotypeTable[haplotypeTable_GAPS] = np.nan
	#distance_haplotypeTable_tem = (cdist(haplotypeTable.T, haplotypeTable.T, 'hamming') - cdist(gapsTable.T, gapsTable.T, 'hamming')) * haplotypeTable.shape[0]
	#distance_haplotypeTable = (cdist(haplotypeTable.T, haplotypeTable.T, 'hamming') - cdist(haplotypeTable_GAPS.T, haplotypeTable_GAPS.T, 'hamming')) * haplotypeTable.shape[0]
	distance_haplotypeTable = (cdist(haplotypeTable.T, haplotypeTable.T, lambda u, v: clean_dis(u, v)))
	for i in range(NHap):
		for j in range(NHap):
			if j > i:
				output_file.write(f'{chr}\t{start}\t{list_haplotypes[i]}\t{list_haplotypes[j]}\t{haplotypeTable.shape[0]}\t{distance_haplotypeTable[i,j]}\n')



