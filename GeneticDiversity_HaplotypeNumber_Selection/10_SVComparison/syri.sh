#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=1:00:00
#SBATCH -J syri


ref=$1
cul=$2
aligment=$3
sample=$4

mkdir -p $sample

syri --nc 3 -F B --cigar --dir $sample -c $aligment -r $ref -q $cul -k --lf $sample.syri.log



