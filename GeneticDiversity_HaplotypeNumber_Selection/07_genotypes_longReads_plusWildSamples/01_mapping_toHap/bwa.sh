#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=6-00:00:00
#SBATCH -J BWA_Mapping

Threads=20

Reference=$1
sample=$2
Read1=$3

ReadGroup="@RG\tID:wHAIPI016410-33\tPL:PacBio\tPU:150211_I191_FCC6L7EANXX_L5_wHAIPI016410-33\tLB:wHAIPI016410-33\tSM:"$sample"\tCN:BGI"

echo "the files are:"
echo $Read1
echo $Reference

bwa mem -x pacbio -t $Threads -R $ReadGroup -M $Reference $Read1 | samtools view -@ 20 -bS - | samtools sort -@ 20 -o $sample.bam -

