#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000mb
#SBATCH --time=12:00:00
#SBATCH -J Example_script

vcf=$1

mkdir -p chr_phy_serial

python /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/vcf2phylip/vcf2phylip.py --input $vcf --fasta --nexus --nexus-binary --min-samples-locus 10 --output-folder ./chr_phy_serial

"done finished"

