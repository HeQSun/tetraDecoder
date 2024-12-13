import sys
import numpy as np
from cyvcf2 import VCF
from collections import defaultdict

# the script produces a table with SampleCultivarID, SampleWIldID, N variance sites considered, number of matching variants, proportion of variants found (presence absence), proportion of variants found but considering ploidy in wild samples, wild clade, wild species.

# mamba activate mypython3
#genotypes_file = str("/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/splitVCF/01/split.000846.vcf.gz")
#output = str("test")

genotypes_file = str(sys.argv[1])
chr = str(sys.argv[2])
output = str(sys.argv[3])

print(f'Processing file: {genotypes_file}')

# processing Genotype table:
vf = VCF(genotypes_file, strict_gt=True)
list_samples = vf.samples
cult_haplotypes_samples = list_samples[20:]
wild_samples = list_samples[:20]

Ncult_haplotypes_samples = len(cult_haplotypes_samples)
Nwild_samples = len(wild_samples)

Nlines = np.zeros(Ncult_haplotypes_samples)
allele_suppost = np.zeros((Ncult_haplotypes_samples, Nwild_samples))
allele_suppost_fq = np.zeros((Ncult_haplotypes_samples, Nwild_samples))

#variant=next(vf)
ini_pos=0
for variant in vf:
    if ini_pos==0:
        ini_pos=variant.POS
    list_gt_total = [item[:-1] for item in variant.genotypes]
    list_gt_cult_haplotypes = list_gt_total[20:]
    list_gt_wild = list_gt_total[:20]
    for index, cul_sample in enumerate(list_gt_cult_haplotypes):
        if cul_sample !=[-1,-1]:
            WSample_SameHaplotype1 = np.array([(sum([(i==cul_sample[0])*1 for i in WSample])>0)*1 for WSample in list_gt_wild])
            WSample_SameHaplotype2 = np.array([(sum([(i==cul_sample[1])*1 for i in WSample])>0)*1 for WSample in list_gt_wild])
            WSample_SameHaplotype = (WSample_SameHaplotype1+WSample_SameHaplotype2)/2
            WSample_SameHaplotypeFq1 = np.array([sum([(i==cul_sample[0])*1 for i in WSample]) for WSample in list_gt_wild])
            WSample_SameHaplotypeFq2 = np.array([sum([(i==cul_sample[1])*1 for i in WSample]) for WSample in list_gt_wild])
            WSample_SameHaplotypeFq = (WSample_SameHaplotypeFq1+WSample_SameHaplotypeFq2)/4
            allele_suppost[index,] += WSample_SameHaplotype
            allele_suppost_fq[index,] += WSample_SameHaplotypeFq
            Nlines[index] += 1

suppost_perWsample_output = open(f'{output}.txt', "w")
suppost_perWsample_output_top = open(f'{output}_top90quantileSamples.txt', "w")

for index, haplotype in enumerate(cult_haplotypes_samples):
  if Nlines[index] != 0:
    pro_allele_suppost = allele_suppost[index,]/Nlines[index]
    pro_allele_suppost_fq = allele_suppost_fq[index,]/(Nlines[index])
    min_Suppost=np.quantile(pro_allele_suppost_fq, 0.90)
    for Windex, WSample in enumerate(wild_samples):
      suppost_perWsample_output_join = '\t'.join(map(str, [chr, str(ini_pos), haplotype, WSample, Nlines[index], allele_suppost[index,Windex], pro_allele_suppost[Windex], pro_allele_suppost_fq[Windex]]))
      suppost_perWsample_output.write(f'{suppost_perWsample_output_join}\n')
      if pro_allele_suppost_fq[Windex]>=min_Suppost:
        suppost_perWsample_output_top.write(f'{suppost_perWsample_output_join}\n')


suppost_perWsample_output.close()
suppost_perWsample_output_top.close()


