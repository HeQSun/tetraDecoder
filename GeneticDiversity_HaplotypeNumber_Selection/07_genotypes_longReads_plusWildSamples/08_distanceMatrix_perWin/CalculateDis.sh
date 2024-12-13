#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00


chr=$1

vcf=/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/merged_vcf/chr/chr$chr""_filterQUAL30.vcf.gz

#vcf=/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/splitVCF/01/split.000200.vcf.gz

# this is martin version but it is quite slow...
python3 /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/genomics_general/VCF_processing/parseVCF.py -i $vcf -o chr$chr""_filterQUAL30.geno.gz


python /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/genomics_general/distMat.py --windType coordinate -w 100000 -m 200  \
--addWindowID \
--windowDataOutFile chr$chr""_filterQUAL30.winData \
-f phased \
--ploidy 2 \
-g chr$chr""_filterQUAL30.geno.gz \
-o chr$chr""_filterQUAL30.dis \
-T 5
# 


