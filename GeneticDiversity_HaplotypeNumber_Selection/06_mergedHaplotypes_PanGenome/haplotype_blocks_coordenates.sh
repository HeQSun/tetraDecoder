#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=48:00:00
#SBATCH -J Fasta_blocks

#chr_used=chr06
#haplotype_name=J_hap24

chr_used=$1
haplotype_name=$2

sample_name=$(echo $haplotype_name | sed 's@_.*@@g' )
single_haplotype_name=$(echo $haplotype_name | sed 's@.*_@@g' )

echo $haplotype_name
echo $chr_used
echo $sample_name
echo $single_haplotype_name

grep -w $haplotype_name Unique_haplotypeGroups_$chr_used"_"winSize200kb_stepSize50kb.txt > Unique_haplotypeGroups_$chr_used"_"winSize200kb_stepSize50kb_$haplotype_name.txt

python3 extract_chr.py  ../../../02_reference/Haplotype_resolved_potato_assemblies_v1_20230327/full48_$sample_name"/"$sample_name"_"full_genome.fasta $chr_used"_"$single_haplotype_name ref_$chr_used"_"$haplotype_name.fasta

#python3 extract_chr.py  DM_1-3_516_R44_potato_genome_assembly.v6.1.fasta $chr_used DMref_$chr_used"_"$haplotype_name.fasta 

minimap2 -t 2 -ax asm5 --eqx DMref_$chr_used.fasta  ref_$chr_used"_"$haplotype_name.fasta | samtools view -@ 2 -bS - | samtools sort -@ 2 - > $haplotype_name"_"$chr_used.output.bam

samtools index $haplotype_name"_"$chr_used.output.bam

if test -f "allNewCoor_$chr_used"__"$haplotype_name.txt"; then
    rm allNewCoor_$chr_used"__"$haplotype_name.txt 
fi

while IFS= read -r line
do
  echo "$line"
  startP=$(echo $line | awk '{print $2}')
  endP=$(echo $line | awk '{print $3}')
  hapID=$(echo $line | awk '{print $4}')  
  segment=$(echo $line | awk '{print $5}')  
  
  startP_S=$(echo $startP | awk '{pos=$1-1000 ; if(pos<1) {print "1"} else {print pos}}')
  hometools mapbp $chr_used":"$startP_S"-"$startP $haplotype_name"_"$chr_used.output.bam | awk -F"-" '{print $2}' > Start_$chr_used"_"$startP"_"$endP"_"$hapID"_"$segment"__"$haplotype_name.txt
  
  endP_S=$(echo $endP | awk '{pos=$1-1000 ; if(pos<1) {print "1"} else {print pos}}')
  hometools mapbp $chr_used":"$endP_S"-"$endP $haplotype_name"_"$chr_used.output.bam | awk -F"-" '{print $2}' > End_$chr_used"_"$startP"_"$endP"_"$hapID"_"$segment"__"$haplotype_name.txt
  
  grep "" *_$chr_used"_"$startP"_"$endP"_"$hapID"_"$segment"__"$haplotype_name.txt | sed 's/.txt:/\t/g' | sed 's/__/\t/g' | sed 's/_/\t/g'  >> allNewCoor_$chr_used"__"$haplotype_name.txt 
  
  rm Start_$chr_used"_"$startP"_"$endP"_"$hapID"_"$segment"__"$haplotype_name.txt End_$chr_used"_"$startP"_"$endP"_"$hapID"_"$segment"__"$haplotype_name.txt 
    
done < Unique_haplotypeGroups_$chr_used"_"winSize200kb_stepSize50kb_$haplotype_name.txt

#rm ref_$chr_used"_"$haplotype_name.fasta DMref_$chr_used"_"$haplotype_name.fasta 




