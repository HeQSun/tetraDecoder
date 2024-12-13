#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=8-00:00:00
#SBATCH -J LD_cal

file=$1
output=$2
min_dis=$3
start_current=$4
end_pos=$5

python3 calculate_LD.py $file $output $min_dis $start_current $end_pos


