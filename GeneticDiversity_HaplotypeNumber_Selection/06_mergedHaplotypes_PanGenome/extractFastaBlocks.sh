#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=48:00:00
#SBATCH -J PFastaBlocks

chr_used=$1
outputfolder=$2

if test -f "allRefSeq_$chr_used"_"winSize200kb_stepSize50kb.fasta"; then
    rm allRefSeq_$chr_used"_"winSize200kb_stepSize50kb.fasta
fi

while IFS= read -r line
do
  echo "$line"
  seq_name=$(echo $line | sed 's/ /_/g')
  sample=$(echo $line | awk '{print $6}')  
  sampleHap=$(echo $line | awk '{print $7}')
  start=$(echo $line | awk '{print $8}')
  end=$(echo $line | awk '{print $9}')
  #haplotype_name=$line
  #sbatch tem_haplotype_blocks_coordenates.sh $chr_used $haplotype_name
  python3 extract_region.py ./ref_haplotypes/"ref_"$chr_used"_"$sample"_"$sampleHap.fasta $chr_used"_"$sampleHap $start $end ./fileTemporal/seq_$seq_name.fasta $seq_name
  cat ./fileTemporal/seq_$seq_name.fasta >> allRefSeq_$chr_used"_"winSize10kb_stepSize5kb.fasta
  rm ./fileTemporal/seq_$seq_name.fasta
done < filtered_merged_allNewCoor_$chr_used"_"winSize10kb_stepSize5kb.txt

mv allRefSeq_$chr_used"_"winSize10kb_stepSize5kb.fasta $outputfolder"/"



