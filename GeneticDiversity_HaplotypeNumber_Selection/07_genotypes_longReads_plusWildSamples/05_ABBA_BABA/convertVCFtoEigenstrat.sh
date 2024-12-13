#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=2-00:00:00
#SBATCH -J ConvertVCFtoEigen


# Script to convert vcf to eigenstrat format for ADMIXTOOLS
# Written by Joana Meier
# It takes a single argument: the vcf file (can be gzipped)

# It requires vcftools and admixtools

# for some clusters, it is needed to load these modules:
# module load bioinfo-tools gcc/4.8.2 vcftools openblas/0.2.13_seq perl/5.18.4 admixtools
# module load bioinfo-tools gcc/4.8.2 vcftools openblas/0.2.20 perl/5.26.2 AdmixTools/5.0-20171024

# if [ $# -ne 1 ]
#  then
#  echo "Please provide the vcf file to parse"
#  exit 1
# fi

file=$1
file=${file%.gz}
file=${file%.vcf}
file_base=${file##*/}
file_path=${file%/*}

# # Removing Samples:
# remove_fam_file=$2
# # if the vcf file is gzipped:
# if [ -s $file.vcf.gz ]
# then
#  # Get a .map and .ped file
#  vcftools --gzvcf $file".vcf.gz" \
#          --plink --mac 1.0 --remove-indels --max-alleles 2 \
#          --remove $remove_fam_file \
#          --out $file_base
# else
# vcftools --vcf $file".vcf" \
#          --plink --mac 1.0 --remove-indels --max-alleles 2 \
#          --remove $remove_fam_file \
#          --out $file_base
# fi

# Without removing Samples:
# if the vcf file is gzipped:
if [ -s $file.vcf.gz ]
then
 # Get a .map and .ped file
 vcftools --gzvcf $file".vcf.gz" \
         --plink --mac 1.0 --remove-indels --max-alleles 2 \
         --out $file_base
else
vcftools --vcf $file".vcf" \
         --plink --mac 1.0 --remove-indels --max-alleles 2 \
         --out $file_base
fi


# Change the .map file to match the requirements of ADMIXTOOLS
awk -F"\t" 'BEGIN{scaff="";add=0}{
        split($2,newScaff,":")
        if(!match(newScaff[1],scaff)){
                scaff=newScaff[1]
                add=lastPos
        }
        count+=0.0001
        pos=add+$4
        print newScaff[1]"\t"$2"\t"count"\t"pos
        lastPos=pos
}' ${file_base}.map  | sed 's/^chr0//' | sed 's/^chr//' > better.map
sed 's/chr0//' better.map | sed 's/chr//' > better2.map
cp ${file_base}.map tem.${file_base}.map
mv better2.map ${file_base}.map

# Change the .ped file to match the ADMIXTOOLS requirements
awk 'BEGIN{ind=1}{printf ind"\t"$2"\t0\t0\t0\t1\t";
 for(i=7;i<=NF;++i) printf $i"\t";ind++;printf "\n"}' ${file_base}.ped > tmp.ped
cp ${file_base}.ped tem.${file_base}.ped
mv tmp.ped ${file_base}.ped


# create an inputfile for convertf
echo "genotypename:    ${file_base}.ped" > par.PED.EIGENSTRAT.${file_base}
echo "snpname:         ${file_base}.map" >> par.PED.EIGENSTRAT.${file_base}
echo "indivname:       ${file_base}.ped" >> par.PED.EIGENSTRAT.${file_base}
echo "outputformat:    EIGENSTRAT" >> par.PED.EIGENSTRAT.${file_base}
echo "genotypeoutname: ${file_base}.eigenstratgeno" >> par.PED.EIGENSTRAT.${file_base}
echo "snpoutname:      ${file_base}.snp" >> par.PED.EIGENSTRAT.${file_base}
echo "indivoutname:    ${file_base}.ind" >> par.PED.EIGENSTRAT.${file_base}
echo "familynames:     NO" >> par.PED.EIGENSTRAT.${file_base}


# Use CONVERTF to parse PED to eigenstrat
convertf -p par.PED.EIGENSTRAT.${file_base}

# change the snp file for ADMIXTOOLS:
awk 'BEGIN{i=0}{i=i+1; print $1"\t"$2"\t"$3"\t"i"\t"$5"\t"$6}' $file_base.snp > $file_base.snp.tmp
cp $file_base.snp tem.$file_base.snp
mv $file_base.snp.tmp $file_base.snp




