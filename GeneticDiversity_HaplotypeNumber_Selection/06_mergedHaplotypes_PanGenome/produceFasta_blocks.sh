#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=48:00:00
#SBATCH -J Fasta_blocks

#chr_used=chr06
#haplotype_name=J_hap24
#outputfolder=test

chr_used=$1
haplotype_name=$2
outputfolder=$3

#haplotype_name=$(echo $haplotype_file | sed 's@.*/01_nucmer_syri_vs_DM/nucmer_syri_@@g' | sed 's@vs_DM/refDM_chr.*_vs_qry_@@g' | sed 's@/syri.out@@g' )

sample_name=$(echo $haplotype_name | sed 's@_.*@@g' )
single_haplotype_name=$(echo $haplotype_name | sed 's@.*_@@g' )

echo $haplotype_name
echo $chr_used
echo $sample_name
echo $single_haplotype_name

grep -w $haplotype_name Unique_haplotypeGroups_$chr_used"_"winSize200kb_stepSize50kb.txt > Unique_haplotypeGroups_$chr_used"_"winSize200kb_stepSize50kb_$haplotype_name.txt

python extract_chr.py  ../../../02_reference/Haplotype_resolved_potato_assemblies_v1_20230327/full48_$sample_name"/"$sample_name"_"full_genome.fasta $chr_used"_"$single_haplotype_name ref_$chr_used"_"$haplotype_name.fasta

while IFS= read -r line
do
  echo "$line"
  startP=$(echo $line | awk '{print $2}')
  endP=$(echo $line | awk '{print $3}')
  
  python3 extract_region_for_rev.py DM_1-3_516_R44_potato_genome_assembly.v6.1.fasta $chr_used $startP $endP $haplotype_name"_"$chr_used"_"$startP"_"$endP.output.fasta $haplotype_name"_"$chr_used"_"$startP"_"$endP.output_rev.fasta $chr_used"_"$single_haplotype_name
  
  minimap2 -ax asm5 --eqx ref_$chr_used"_"$haplotype_name.fasta $haplotype_name"_"$chr_used"_"$startP"_"$endP.output.fasta > $haplotype_name"_"$chr_used"_"$startP"_"$endP.output.sam
  
  python3 /dss/dsshome1/lxc03/di36guz2/miniconda3/envs/syri_env/bin/syri -c $haplotype_name"_"$chr_used"_"$startP"_"$endP.output.sam -r ref_$chr_used"_"$haplotype_name.fasta -q $haplotype_name"_"$chr_used"_"$startP"_"$endP.output.fasta -k -F S --prefix $haplotype_name"_"$chr_used"_"$startP"_"$endP"_F_"
  
  minimap2 -ax asm5 --eqx ref_$chr_used"_"$haplotype_name.fasta $haplotype_name"_"$chr_used"_"$startP"_"$endP.output_rev.fasta > $haplotype_name"_"$chr_used"_"$startP"_"$endP.output.sam
  
  python3 /dss/dsshome1/lxc03/di36guz2/miniconda3/envs/syri_env/bin/syri -c $haplotype_name"_"$chr_used"_"$startP"_"$endP.output.sam -r ref_$chr_used"_"$haplotype_name.fasta -q $haplotype_name"_"$chr_used"_"$startP"_"$endP.output_rev.fasta -k -F S --prefix $haplotype_name"_"$chr_used"_"$startP"_"$endP"_R_"
  
  cat $haplotype_name"_"$chr_used"_"$startP"_"$endP"_"[F,R]_syri.out | awk '{if($8-$7>50 && $3-$2>50 && $3-$2<300000) print $1"\t"$2"\t"$3"\t"$6"\t"$7"\t"$8"\t"$3-$2"\t"$8-$7}' | sort -k8,8n | tail -n 1 | sed 's/\t/ /g' > $haplotype_name"_"$chr_used"_"$startP"_"$endP"_"NewCoordinates.txt

  newCoord=$(echo $line $(cat $haplotype_name"_"$chr_used"_"$startP"_"$endP"_"NewCoordinates.txt))
  
  echo $newCoord >> NewCoor_Unique_haplotypeGroups_$chr_used"_"winSize200kb_stepSize50kb_$haplotype_name.txt
  
  nameNSeq=$(echo $newCoord | awk -F" " '{print $1"__"$2"__"$3"__"$4"__"$5"__"$6}')
  Nstart=$(echo $newCoord | awk -F" " '{if($8<$9) {print $8} else {print $9}}')
  Nend=$(echo $newCoord | awk -F" " '{if($8<$9) {print $9} else {print $8}}')
  
  python3 extract_region.py ref_$chr_used"_"$haplotype_name.fasta $chr_used"_"$single_haplotype_name $Nstart $Nend $outputfolder"/"NSeq_$haplotype_name"_"$chr_used"_"$startP"_"$endP.fasta $nameNSeq
  
  rm $haplotype_name"_"$chr_used"_"$startP"_"$endP""*  

done < Unique_haplotypeGroups_$chr_used"_"winSize200kb_stepSize50kb_$haplotype_name.txt




