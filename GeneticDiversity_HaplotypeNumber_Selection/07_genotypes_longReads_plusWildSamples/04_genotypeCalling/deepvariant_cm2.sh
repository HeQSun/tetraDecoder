#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH --cpus-per-task=8
#SBATCH --mem=15000mb
#SBATCH --time=20:00:00
#SBATCH -J deepvariant

#echo $TMPDIR

export TMPDIR=/var/tmp/

module use /dss/dsslegfs02/pn73se/pn73se-dss-0000/spack/modules/$LRZ_INSTRSET/linux* ; module load user_spack
module load charliecloud

export OMP_NUM_THREADS=8


sample=$1
run=$2
region=$3
#regions.bed

mkdir -p ./vcf
mkdir -p ./vcf/cm_$run

#export TMPDIR=$PWD/tmpD/$sample.$run
#echo $TMPDIR

ch-run --bind $PWD:/bedfiles/ --bind $PWD/vcf/cm_$run:/output/ --bind /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/02_reference/DM_1-3_516_R44_potato_genome_assembly.v6.1:/reference_file --bind /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/03_mapping_ReferenceDM:/bam_files -w /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/deepvariant/deepvariant_1_4_0_charliecloud_container -- /opt/deepvariant/bin/run_deepvariant \
    --model_type PACBIO \
    --ref /reference_file/DM_1-3_516_R44_potato_genome_assembly.v6.1.sm.fa \
    --reads /bam_files/$sample.bam \
    --output_vcf /output/$sample.$run.vcf.gz \
    --output_gvcf /output/$sample.$run.gvcf.gz \
    --regions /bedfiles/$region \
    --num_shards 8 
    #--intermediate_results_dir /tmp


