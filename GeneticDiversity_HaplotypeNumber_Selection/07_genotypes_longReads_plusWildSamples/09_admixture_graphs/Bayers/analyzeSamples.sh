#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=3:00:00

input=$1
output=$2

python /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/AdmixtureBayes/admixturebayes/analyzeSamples.py --mcmc_results $input --result_file $output.csv 

echo "finished analyzeSamples"


python /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/AdmixtureBayes/admixturebayes/makePlots.py --plot estimates --posterior $output.csv --write_rankings rank.$output.txt --output_prefix Estimate_$output

echo "finished Make Plots"



