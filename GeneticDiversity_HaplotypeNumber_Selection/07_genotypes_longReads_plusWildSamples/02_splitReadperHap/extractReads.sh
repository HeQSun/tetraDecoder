#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00
#SBATCH -J ExtractReads


bam_file=$1
list_chr_haplotype=$2
haplotype_name=$3

#samtools view -@ 4 -bL $list_chr_haplotype -o $haplotype_name.bam $bam_file
#samtools sort -@ 4 -n -o $haplotype_name.sort.bam $haplotype_name.bam
#samtools fastq -@ 4 $haplotype_name.sort.bam > $haplotype_name.fastq


samtools view -@ 5 -bL $list_chr_haplotype $bam_file | samtools sort -@ 5 - | samtools fastq -@ 5 - > $haplotype_name.fastq

gzip $haplotype_name.fastq


