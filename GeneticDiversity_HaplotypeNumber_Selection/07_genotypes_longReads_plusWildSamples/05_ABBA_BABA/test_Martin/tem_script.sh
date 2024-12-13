#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00
#SBATCH -J Dis_M



python ../genomics_general-master/distMat.py --windType coordinate -w 50000 -m 250  \
--addWindowID \
--windowDataOutFile test_chr12_filterQUAL30.winData \
-f phased \
--ploidy 2 \
-g test_chr12_filterQUAL30.geno.gz \
-o test_chr12_filterQUAL30.dis \
-T 5



