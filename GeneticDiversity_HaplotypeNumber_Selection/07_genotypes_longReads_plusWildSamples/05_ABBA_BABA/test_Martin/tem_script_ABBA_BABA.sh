#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00
#SBATCH -J ABBA_BABA



python ../genomics_general-master/ABBABABAwindows.py \
-g test_chr12_filterQUAL30.geno.gz -f phased \
-o test_chr12_filterQUAL30.ABBABABA_test.csv.gz \
-P1 SRR15458982 -P2 SRR15458984 -P3 SampleE_Hap_1 -O SRR15458990 \
--popsFile chr12_filterQUAL30.pop.txt -w 100000 -m 250 --T 5





