#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00
#SBATCH -J clusterHaplotypes

genotype_file=$1
win_size=$2
output_name=$3
maxdif=$4

# python3 /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/02_genotypes/AllSamplesPopPar_Win10000_chr01_missing_genotypeTable.txt 10000 Cluster_hapIDs_chr01
python3 ClusteringHaplotypes.py $genotype_file $win_size $output_name $maxdif



