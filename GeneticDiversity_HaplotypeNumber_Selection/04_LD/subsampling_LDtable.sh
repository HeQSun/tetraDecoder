#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4700mb
#SBATCH --time=8-00:00:00
#SBATCH -J SubSample_LD


max_distance=$1
max_distance_name=$2
sampling_lines=$3
chr="chr"$4

#for i in $(seq -w 1 12)
#i=01
#chr="chr"$i
#chrSize=$(grep $chr chr_size.txt | awk '{print $2}') ; echo $chrSize
#ProSampling=$(echo $chrSize | awk -v SSize=$sampling_lines -v maxD=/$max_distance '{print SSize/(((($1/500)-1)*($1/500))/2)}') ; echo $ProSampling
# cat ./LD_NoSingletons_$chr"_maxDis"$max_distance_name"/LD_"$chr"_"* | awk -v PROB=$ProSampling -v CHR=$chr -v distanceN=$max_distance_name 'BEGIN  {srand()} !/^$/ { if (rand() <= PROB || FNR==1) print > "zLD_"CHR"_NoSingletons_"distanceN".txt"}' 
cat ./LD_NoSingletons_$chr"_maxDis"$max_distance_name"/LD_"$chr"_"* > total_LD_"$chr"_"maxDis"$max_distance_name.subSample$sampling_lines.txt

python3 sampleLines.py -n $sampling_lines --file total_LD_"$chr"_"maxDis"$max_distance_name.subSample$sampling_lines.txt > all_LD_"$chr"_"maxDis"$max_distance_name.subSample$sampling_lines.txt

awk '{print $1"\t"$3-$2"\t"$9}' all_LD_"$chr"_"maxDis"$max_distance_name.subSample$sampling_lines.txt > dis_all_LD_"$chr"_"maxDis"$max_distance_name.subSample$sampling_lines.txt

rm total_LD_"$chr"_"maxDis"$max_distance_name.subSample$sampling_lines.txt
#done



