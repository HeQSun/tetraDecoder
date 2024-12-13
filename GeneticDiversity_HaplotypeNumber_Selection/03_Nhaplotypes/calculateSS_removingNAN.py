import os
import sys
import numpy as np
import math

#genotype_file = str("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesMultiAlleles_chr01_missing_genotypeTable.txt")
#win_size = int(10000)
#output_name = str("test")

genotype_file = str(sys.argv[1])
win_size = int(sys.argv[2])
output_name = str(sys.argv[3])

output_file = open(f'{output_name}.txt', "w")

def str_to_int(x):
        if x!='nan':
                return int(x)
        else:
                return int(-1)

for line in open(genotype_file, "r"):
        if line[0:3] == "CHR":
                haplotypes = line.strip().split("\t")[5:]
                SS_array = np.zeros((len(haplotypes), len(haplotypes)))
                NVar = 0
                win_current = 0
        else:
                line_values = line.strip().split("\t")
                pos = int(line_values[1])
                win = math.ceil(int(pos) / win_size) * win_size
                genotypes = np.asarray([str_to_int(x) for x in line_values[5:]])
                SS_var = genotypes[:,None] != genotypes
                NA_genotypes_boolean = (genotypes==(-1))
                SS_var[NA_genotypes_boolean, :] = False
                SS_var[:, NA_genotypes_boolean] = False
                if win == win_current:
                        SS_array += SS_var
                        NVar += 1
                else:
                        for i,x in enumerate(haplotypes):
                                for j,y in enumerate(haplotypes):
                                        if i > j:
                                                output_file.write(f'{win_current}\t{x}\t{y}\t{int(SS_array[i,j])}\n')
                        SS_array = SS_var*1
                        NVar = 0
                        win_current = win

output_file.close()


