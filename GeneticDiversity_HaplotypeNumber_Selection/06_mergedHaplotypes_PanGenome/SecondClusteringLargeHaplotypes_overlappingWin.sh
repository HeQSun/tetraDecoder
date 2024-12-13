#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=12:00:00
#SBATCH -J Second_Clustering

# input = 'haplotypeGroups_min50SS_10KbWin_chr06.txt'
# winSize=int(1000000)
# step_size = int(50000)
# output="test"
# chr="chr06"
# pro_max_d = float(0.25)

input=$1
winSize=$2
step_size=$3
output=$4
chr=$5
pro_max_d=$6

# calculate pairwise difference between large windows:
python3 SecondClusteringLargeHaplotypes_overlappingWin.py $input $winSize $step_size $output $chr $pro_max_d




