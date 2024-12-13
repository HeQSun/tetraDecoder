#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00
#SBATCH -J vcf2phy_chr


vcf=$1

mkdir -p chr_phy

python /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/vcf2phylip/vcf2phylip.py --input $vcf --fasta --nexus --nexus-binary --min-samples-locus 10 --output-folder ./chr_phy

"done finished"

