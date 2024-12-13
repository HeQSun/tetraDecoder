import os
import sys
import numpy as np
import pandas as pd

import scipy.cluster.hierarchy as shc
from scipy.spatial.distance import pdist, squareform
import math

# input = 'haplotypeGroups_min50SS_10KbWin_chr07.txt'
# winSize=int(1000000)
# step_size = int(50000)
# output="test"
# chr="chr07"
# pro_max_d = float(0.25)

input = str(sys.argv[1])
winSize = int(sys.argv[2])
step_size = int(sys.argv[3])
output = str(sys.argv[4])
chr = str(sys.argv[5])
pro_max_d = float(sys.argv[6])

table = pd.read_table(input, sep="\t", names=['pos', 'sample', 'haplotype', 'missingP'])  
tableSpread = table[['pos', 'sample', 'haplotype']].pivot(index = 'pos', columns = "sample", values = "haplotype").reset_index().rename_axis("", axis = 1)
tableSpread_missingP = table[['pos', 'sample', 'missingP']].pivot(index = 'pos', columns = "sample", values = "missingP").reset_index().rename_axis("", axis = 1)

output_file = open(f'{output}.txt', "w")
output_file.write(f'chr\tpos\tsample\thaplotype\tMean_NAN_haplo\n')
output_file.close()

end = max(tableSpread['pos'])
list_haplotypes = list(tableSpread.columns)[1:]
NHap = len(list_haplotypes)

max_d = int((winSize/10000)*pro_max_d)
table_hap_groups_pre = pd.DataFrame()


for start in range(0,end,step_size):
  tableSpreadWin = tableSpread[(tableSpread['pos'] >= start) & (tableSpread['pos'] < start+winSize)]
  tableSpreadWin_missingP = tableSpread_missingP[(tableSpread_missingP['pos'] >= start) & (tableSpread_missingP['pos'] < start+winSize)]
  if tableSpreadWin.shape[0] > 0:
    total_missingP_win = [x.mean() for x in tableSpreadWin_missingP[list(tableSpreadWin_missingP.columns)[1:]].to_numpy().T]
    SS_array = np.zeros((NHap, NHap))
    tableSpreadWin_missingP_np = tableSpreadWin_missingP[list(tableSpreadWin.columns)[1:]].to_numpy()>0.9
    pairDistance = [row[:,None] != row for row in tableSpreadWin[list(tableSpreadWin.columns)[1:]].to_numpy()]
    for i in range(len(pairDistance)):
      new_dis = pairDistance[i]
      new_dis[tableSpreadWin_missingP_np[i], :] = False
      new_dis[:, tableSpreadWin_missingP_np[i]] = False
      SS_array += new_dis
    condensed_dist = squareform(SS_array)
    linkresult = shc.linkage(condensed_dist, method='complete')
    clusters = shc.fcluster(linkresult, max_d, criterion='distance')
    if not table_hap_groups_pre.empty:
      max_cluster_pre = table_hap_groups_pre['haplotype_group_pre'].max()
      table_hap_groups = pd.DataFrame(data={'labs': list(list_haplotypes), 'haplotype_group': list(clusters+max_cluster_pre), 'Mean_NAN_haplo': list(total_missingP_win)})
      comparisonDF = table_hap_groups.merge(table_hap_groups_pre, on='labs', how='left').dropna().groupby(['haplotype_group', 'haplotype_group_pre']).agg({'count'})['labs'].reset_index()
      idx = comparisonDF.groupby('haplotype_group')['count'].transform(max) == comparisonDF['count']
      idx2 = comparisonDF[idx].groupby('haplotype_group_pre')['count'].transform(max) == comparisonDF[idx]['count']
      previous_hap_groups = comparisonDF[idx][idx2].groupby('haplotype_group_pre').first().reset_index().groupby('haplotype_group').first().reset_index()
      table_hap_groups_merge = table_hap_groups.merge( previous_hap_groups[['haplotype_group', 'haplotype_group_pre']], on='haplotype_group', how='left')
      table_hap_groups_merge['haplotype_group_ed'] = [int(table_hap_groups_merge['haplotype_group'][x]) if math.isnan(table_hap_groups_merge['haplotype_group_pre'][x]) else int(table_hap_groups_merge['haplotype_group_pre'][x]) for x in range(0, table_hap_groups_merge.shape[0])]
      table_hap_groups_merge['chr'] = chr
      table_hap_groups_merge['pos'] = start
      table_hap_groups_merge[['chr', 'pos', 'labs', 'haplotype_group_ed', 'Mean_NAN_haplo']].to_csv(f"{output}.txt", sep="\t", index=False, header=False, mode='a')
      table_hap_groups_pre = table_hap_groups_merge[['labs', 'haplotype_group_ed']].rename(columns={"labs": "labs", "haplotype_group_ed": "haplotype_group_pre"})
    else:
      table_hap_groups = pd.DataFrame(data={'labs': list(list_haplotypes), 'haplotype_group': list(clusters), 'Mean_NAN_haplo': list(total_missingP_win)})
      table_hap_groups['chr'] = chr
      table_hap_groups['pos'] = start
      table_hap_groups[['chr', 'pos', 'labs', 'haplotype_group', 'Mean_NAN_haplo']].to_csv(f"{output}.txt", sep="\t", index=False, header=False, mode='a')
      table_hap_groups_pre = table_hap_groups[['labs', 'haplotype_group']].rename(columns={"labs": "labs", "haplotype_group": "haplotype_group_pre"})
  else:
    table_hap_groups['chr'] = chr
    table_hap_groups['pos'] = start
    table_hap_groups[['chr', 'pos', 'labs', 'haplotype_group', 'Mean_NAN_haplo']].to_csv(f"{output}.txt", sep="\t", index=False, header=False, mode='a')






