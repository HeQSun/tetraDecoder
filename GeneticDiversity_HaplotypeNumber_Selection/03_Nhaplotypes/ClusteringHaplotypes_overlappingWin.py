import os
import sys
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as shc
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
import math

# genotype_file = str("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr01_missing_genotypeTable.txt")
# win_size = int(10000)
# output_name = str("test")
# max_d = int(10)

# this script takes a genotype table, a window size and and output name and: produce a distance matrix and cluster haplotypes with a maximum of around 10SNPs (max_d parameter).
# the output is a table with window start position, sample and haplotype ID after clustering.
# The scripts uses the haplotype groups produced from the previous window to allocate haplotypes ID in the current window. 

# genotype_file = str("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr01_missing_genotypeTable.txt")
# win_size = int(1000000)
# step_size= int(50000)
# max_d = int(2000)
# output_name = str("test")


genotype_file = str(sys.argv[1])
win_size = int(sys.argv[2])
step_size=int(sys.argv[3])
max_d = int(sys.argv[4])
output_name = str(sys.argv[5])

def str_to_int(x):
    if x!='nan':
        return int(x)
    else:
        return int(-1)


dic_win_SS_array = {}
dic_win_NA_genotypes_mean = {}
dic_win_NVar = {}

for line in open(genotype_file, "r"):
    if line[0:3] == "CHR":
        haplotypes = line.strip().split("\t")[5:]
        list_win_pre = [0] 
        dic_win_SS_array[0] = np.zeros((len(haplotypes), len(haplotypes)))
        dic_win_NA_genotypes_mean[0] = np.zeros((len(haplotypes)))
        dic_win_NVar[0] = 0
        table_hap_groups_pre = pd.DataFrame()
        if os.path.exists(f"{output_name}.txt"):
            os.remove(f"{output_name}.txt")
    else:
        line_values = line.strip().split("\t")
        pos = int(line_values[1])
        list_win = [i for i in range(int((pos-win_size)/step_size),int(pos/step_size)+1) if i>=0]
        list_win_withPre = set(list_win_pre + list_win)
        genotypes = np.asarray([str_to_int(x) for x in line_values[5:]])
        SS_var = genotypes[:,None] != genotypes
        NA_genotypes = (genotypes==(-1))*1
        NA_genotypes_boolean = (genotypes==(-1))
        SS_var[NA_genotypes_boolean, :] = False
        SS_var[:, NA_genotypes_boolean] = False
        #SS_var[NA_genotypes_boolean, :][:, NA_genotypes_boolean] = False
        #SS_var[NA_genotypes_boolean,] = False
        #SS_var[:,NA_genotypes_boolean] = False
        for win_current in list_win_withPre:
            if not win_current in dic_win_SS_array:
                dic_win_SS_array[win_current] = np.zeros((len(haplotypes), len(haplotypes)))
                dic_win_NA_genotypes_mean[win_current] = np.zeros((len(haplotypes)))
                dic_win_NVar[win_current] = 0
            if win_current in list_win:
                dic_win_SS_array[win_current] += SS_var
                dic_win_NA_genotypes_mean[win_current] += NA_genotypes
                dic_win_NVar[win_current] += 1
            else:
                condensed_dist = squareform(dic_win_SS_array[win_current])
                linkresult = shc.linkage(condensed_dist, method='complete')
                clusters = shc.fcluster(linkresult, max_d, criterion='distance')
                #Z = shc.ward(pdist(dic_win_SS_array[win_current]))
                #clusters = shc.fcluster(Z, max_d, criterion='distance')
                if dic_win_NVar[win_current] > 0:
                    frac_NA_genotypes_mean = dic_win_NA_genotypes_mean[win_current]/dic_win_NVar[win_current]
                else:
                    frac_NA_genotypes_mean = dic_win_NA_genotypes_mean[win_current]
                if not table_hap_groups_pre.empty:
                    max_cluster_pre = table_hap_groups_pre['haplotype_group_pre'].max()
                    table_hap_groups = pd.DataFrame(data={'labs': list(haplotypes), 'haplotype_group': list(clusters+max_cluster_pre), 'Mean_NAN_haplo': list(frac_NA_genotypes_mean)})
                    comparisonDF = table_hap_groups.merge(table_hap_groups_pre, on='labs', how='left').dropna().groupby(['haplotype_group', 'haplotype_group_pre']).agg({'count'})['labs'].reset_index()
                    idx = comparisonDF.groupby('haplotype_group')['count'].transform(max) == comparisonDF['count']
                    idx2 = comparisonDF[idx].groupby('haplotype_group_pre')['count'].transform(max) == comparisonDF[idx]['count']
                    previous_hap_groups = comparisonDF[idx][idx2].groupby('haplotype_group_pre').first().reset_index().groupby('haplotype_group').first().reset_index()
                    table_hap_groups_merge = table_hap_groups.merge( previous_hap_groups[['haplotype_group', 'haplotype_group_pre']], on='haplotype_group', how='left')
                    table_hap_groups_merge['haplotype_group_ed'] = [int(table_hap_groups_merge['haplotype_group'][x]) if math.isnan(table_hap_groups_merge['haplotype_group_pre'][x]) else int(table_hap_groups_merge['haplotype_group_pre'][x]) for x in range(0, table_hap_groups_merge.shape[0])]
                    table_hap_groups_merge['win'] = win_current*step_size
                    table_hap_groups_merge[['win', 'labs', 'haplotype_group_ed', 'Mean_NAN_haplo']].to_csv(f"{output_name}.txt", sep="\t", index=False, header=False, mode='a')
                    table_hap_groups_pre = table_hap_groups_merge[['labs', 'haplotype_group_ed']].rename(columns={"labs": "labs", "haplotype_group_ed": "haplotype_group_pre"})
                    del dic_win_SS_array[win_current]
                    del dic_win_NA_genotypes_mean[win_current]
                    del dic_win_NVar[win_current] 
                else:
                    table_hap_groups = pd.DataFrame(data={'labs': list(haplotypes), 'haplotype_group': list(clusters), 'Mean_NAN_haplo': list(frac_NA_genotypes_mean)})
                    table_hap_groups['win'] = win_current*step_size
                    table_hap_groups[['win', 'labs', 'haplotype_group', 'Mean_NAN_haplo']].to_csv(f"{output_name}.txt", sep="\t", index=False, header=False, mode='a')
                    table_hap_groups_pre = table_hap_groups[['labs', 'haplotype_group']].rename(columns={"labs": "labs", "haplotype_group": "haplotype_group_pre"})
                    del dic_win_SS_array[win_current]
                    del dic_win_NA_genotypes_mean[win_current]
                    del dic_win_NVar[win_current]
        list_win_pre = list_win


