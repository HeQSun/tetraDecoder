# Here we align initial assembly of hifiasm v0p7 to DM to find 12 linkage groups (with four haplotypes still mixed)
#    DM sequence: /netscratch/dep_mercier/grp_schneeberger/projects/Potato_multipleCultivars/s2_10_Cultivars_PacBio_HiFi/zaux_reference_dm/reference_sequence/DM_1-3_516_R44_potato_genome_assembly.v6.1.fa
#    caution: paths and/or file names of DM and cultivar assemblies should be updated.
#    customized tool: ref_linkage_grouper

# step 1. align full (purged) assembly to reference dm

wd=/your/working/path/
cd $wd
mkdir a4_alignment_based_linkage_grouping

# sample O is otava
#
# 1.1 alignment using DM as reference 
#  
for sample in A B C D E F G H I J O; do 
   cd ${wd}/a4_alignment_based_linkage_grouping/
   mkdir zhifiasm_v0p7_dm_${sample}_asm
   cd zhifiasm_v0p7_dm_${sample}_asm
   dm_ref=../../zaux_reference_dm/reference_sequence/DM_1-3_516_R44_potato_genome_assembly.v6.1.fa
   subset_id=4  # update this 
   assembly=../../a2_initial_assembly/sample_${sample}/hifiasm_asm_v0p7/subset${subset_id}_illu_re_align_purge_ovl/DNB_align/clipped${subset_id}_${sample}_hifiasm.p_utg.gfa.fa
   ll ${dm_ref} ${assembly}
   thread=4
   bsub -q multicore20 -n ${thread} -R "span[hosts=1] rusage[mem=48000]" -M 48000 -o minimap2.log -e minimap2.err "minimap2 -ax asm20 -t ${thread} ${dm_ref} ${assembly} | samtools view -@ ${thread} -bS - | samtools sort -@ ${thread} -o ${sample}_against_dm.bam -"
done

# 1.2 alignment using Otava Stie1 as reference

for sample in A B C D E F G H I J O; do 
   cd ${wd}/a4_alignment_based_linkage_grouping/
   mkdir zhifiasm_v0p7_stie1_${sample}_asm
   cd zhifiasm_v0p7_stie1_${sample}_asm
   st1_ref=../../zaux_reference_stie1/reference_genome/Stie1_genome.fa
   subset_id=4
   assembly=../../a2_initial_assembly/sample_${sample}/hifiasm_asm_v0p7/subset${subset_id}_illu_re_align_purge_ovl/DNB_align/clipped${subset_id}_${sample}_hifiasm.p_utg.gfa.fa
   ll ${st1_ref} ${assembly}
   bsub -q multicore20 -n 4 -R "span[hosts=1] rusage[mem=24000]" -M 24000 -o minimap2.log -e minimap2.err "minimap2 -cx asm20 -t 2 ${st1_ref} ${assembly} > ${sample}_against_stie1.paf"
done

# step 2. separate alignment to linkage group (=> 12 groups)

cd ${wd}/a4_alignment_based_linkage_grouping/
for sample in A B C D E F G H I J O; do 
   cd ${wd}/a4_alignment_based_linkage_grouping/
   # against dm
   cd zhifiasm_v0p7_dm_${sample}_asm
   ref_linkage_grouper ${sample}_against_dm.paf ../../zaux_reference_dm/reference_sequence/DM_1-3_516_R44_potato_genome_assembly.v6.1_main12.chrids &> DM_based.log&
   cd ${wd}/a4_alignment_based_linkage_grouping/
   #
done

# step 3. get grouping details and summary

cd ${wd}/a4_alignment_based_linkage_grouping/
for sample in A B C D E F G H I J O; do 
#for sample in C; do
   # against dm
   cd ${wd}/a4_alignment_based_linkage_grouping/   
   cd zhifiasm_v0p7_dm_${sample}_asm
   grep 'res' DM_based.log   | sed 's/   //g' > dm_res_grouping_details_${sample}.txt
   grep 'final' DM_based.log | sed 's/   //g' > dm_res_total_group_size_${sample}.txt
   #
done

# step 4. visualize - see the example pdf; update the script if you would run it: 1. you do not have grouping using stie1 genome; 2. you have only two assemblies; 3. related paths

cd ${wd}/a4_alignment_based_linkage_grouping/
Rscript viz_s2_a4_ref_based_grouping_res_hifiasm_v0p7.R


# step 5. get chr-wise group of contigs for omni-c read extraction - this is for omnic_read_extracter

cd ${wd}/a4_alignment_based_linkage_grouping/
for sample in A B C D E F G H I J O; do 
#for sample in C; do
   # against dm
   cd ${wd}/a4_alignment_based_linkage_grouping/   
   cd zhifiasm_v0p7_dm_${sample}_asm
   #
   for lg in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12; do 
       awk -v chr="$lg" '$7==chr ' dm_res_grouping_details_${sample}.txt | cut -d' ' -f2,4,7,11 | sed 's/ /\t/g' > dm_res_grouping_details_${sample}_${lg}.txt
   done
   #
   ll  dm_res_grouping_details_${sample}_*.txt   
   cat dm_res_grouping_details_${sample}_*.txt | cut -d' ' -f1 | sort | uniq -c | wc -l 
   cat dm_res_grouping_details_${sample}_*.txt | wc -l
   #
done









