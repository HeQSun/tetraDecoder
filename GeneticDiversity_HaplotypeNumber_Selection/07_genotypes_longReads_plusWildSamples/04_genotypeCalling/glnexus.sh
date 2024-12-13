#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --cpus-per-task=12
#SBATCH --mem=40000mb
#SBATCH --time=3-00:00:00
#SBATCH -J glnexus


module use /dss/dsslegfs02/pn73se/pn73se-dss-0000/spack/modules/$LRZ_INSTRSET/linux* ; module load user_spack
module load charliecloud

export TMPDIR=/var/tmp/
export OMP_NUM_THREADS=12


run=$1

mkdir -p ./merged_vcf
mkdir -p ./merged_vcf/$run


ch-run --bind $SCRATCH:/temporal --bind $PWD/vcf/$run:/gvcf/ --bind $PWD/merged_vcf/$run:/output_mergedvcf/ -w /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_raul/glnexus/glnexus_image \
    -- bash -c "/usr/local/bin/glnexus_cli \
    --config DeepVariant \
    --mem-gbytes 40 \
    --threads 12 \
    --dir /temporal/$run \
    /gvcf/*.g.vcf.gz > /output_mergedvcf/TChrGroup_$run.bcf" 

# | bcftools view -Oz - > /output_mergedvcf/ChrGroup_$run.vcf.gz && \
#    tabix -p vcf /output_mergedvcf/ChrGroup_$run.vcf.gz"



