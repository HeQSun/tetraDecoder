#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00
#SBATCH -J IQTree_Chr

phy=$1

mkdir -p IQTree_Chr

phy_file=${phy##*/}
phy_file_base=${phy_file%.phy}
iqtree -s $phy -m GTR -B 1000 -T 4 --prefix ./IQTree_Chr/$phy_file_base

echo "done finished"



