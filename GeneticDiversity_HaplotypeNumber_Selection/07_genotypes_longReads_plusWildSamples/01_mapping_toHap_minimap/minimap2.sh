#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00
#SBATCH -J Minimap2

Reference=$1
sample=$2
Read1=$3

minimap2 -ax map-hifi -t 20 -N 1 --secondary=no $Reference $Read1 | samtools view -@ 20 -bS - | samtools sort -@ 20 -o $sample.bam -

samtools index -@ 20 $sample.bam


