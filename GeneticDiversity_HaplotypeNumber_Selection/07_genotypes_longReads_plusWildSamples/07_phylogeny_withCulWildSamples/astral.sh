#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00
#SBATCH -J Astral


genetrees=$1.treefile
spptree=astral.$1.treefile
output1=astral.$1.scored.treefile
output2=astral.$1.poly.treefile
log1=astral.$1.log
log2=astral.$1.scored.log
log3=astral.$1.poly.log


# Astral
java -Xmx15G -D"java.library.path=/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/SITG_packages/miniconda3/envs/ete3/share/astral-tree-5.7.8-1/lib/" -jar /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/SITG_packages/miniconda3/envs/ete3/share/astral-tree-5.7.8-1/astral.5.7.8.jar -i $genetrees -o $spptree 2> $log1

# scoring_astral:
java -Xmx15G -D"java.library.path=/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/SITG_packages/miniconda3/envs/ete3/share/astral-tree-5.7.8-1/lib/" -jar /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/SITG_packages/miniconda3/envs/ete3/share/astral-tree-5.7.8-1/astral.5.7.8.jar -q $spptree -i $genetrees -o $output1 -t 8 2> $log2

# polytomy_test_astral:
java -Xmx15G -D"java.library.path=/dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/SITG_packages/miniconda3/envs/ete3/share/astral-tree-5.7.8-1/lib/" -jar /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/SITG_packages/miniconda3/envs/ete3/share/astral-tree-5.7.8-1/astral.5.7.8.jar -q $spptree -i $genetrees -o $output2 -t 10 2> $log3




