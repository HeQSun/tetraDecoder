#!/bin/bash
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2:00:00
#SBATCH -J msyd


samples=$1
output=$2

msyd call -i $samples -o $output


