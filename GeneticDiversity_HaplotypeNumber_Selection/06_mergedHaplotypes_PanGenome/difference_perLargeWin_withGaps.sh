#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=12:00:00
#SBATCH -J diffLargeWin


# input='/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/03_Nhaplotypes/haplotypeGroups_min10SS_chr06.txt'
# winSize=1000000
# step_size=50000
# output=test
# chr=chr06

input=$1
winSize=$2
step_size=$3
output=$4
chr=$5

# calculate pairwise difference between large windows:
python3 difference_perLargeWin_withGaps.py $input $winSize $step_size $output $chr



