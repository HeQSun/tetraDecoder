#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=10-00:00:00
#SBATCH -J LD_cal

file=$1
output=$2
start_current=$3
end_pos=$4
max_distance=$5

python3 calculate_LD_withoutSingletons.py $file $output $start_current $end_pos $max_distance


