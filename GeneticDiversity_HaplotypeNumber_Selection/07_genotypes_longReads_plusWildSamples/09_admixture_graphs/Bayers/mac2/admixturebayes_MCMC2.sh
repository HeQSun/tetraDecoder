#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --cpus-per-task=8
#SBATCH --mem=15000mb
#SBATCH --time=3-00:00:00
#SBATCH -J admixturebayes


input=$1
output=$2

python /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/AdmixtureBayes/admixturebayes/runMCMC.py --input_file $input --outgroup C12 --n 450000 --result_file chain1.$output --MCMC_chains 8

echo "finished chain 1"

python /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/AdmixtureBayes/admixturebayes/runMCMC.py --input_file $input --outgroup C12 --n 450000 --result_file chain2.$output --MCMC_chains 8

echo "finished chain 2"


python /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/AdmixtureBayes/admixturebayes/runMCMC.py --input_file $input --outgroup C12 --n 450000 --result_file chain3.$output --MCMC_chains 8

echo "finished chain 3"

