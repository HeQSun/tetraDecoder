#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=4-00:00:00
#SBATCH -J IQTree_Chr


chr=$1

mkdir -p IQTree
mkdir -p IQTree/$chr

ls -1 /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/07_phylogeny_withCulWildSamples/splitVCF_phy/$chr"/"*.phy > list_phy_$chr.txt

for phy in $(cat list_phy_$chr.txt); do echo $phy
phy_file=${phy##*/}
phy_file_base=${phy_file%.phy}
iqtree -s $phy -m GTR -B 1000 -T 4 --prefix ./IQTree/$chr"/"$phy_file_base
done

"done finished"


