#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=12:00:00
#SBATCH -J Pop_Par

genotype_table=$1
chr_used=$2
Nsamples=$3
win_size=$4
prefix_output=$5

#genotype_table="AllSamplesPopPar_Win10000_chr01_missing_genotypeTable.txt"
#chr_used="chr01"
#Nsamples=int(48)
#win_size=int(10000)
#prefix_output="test"

# calculate population parameters:
python3 geneticDiversityCal_addmissingData.py $genotype_table $chr_used $Nsamples $win_size $prefix_output



