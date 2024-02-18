# Here we check quality of hifi reads.

wd=/your/working/path/
# e.g.,
# wd=/netscratch/dep_mercier/grp_schneeberger/projects/Potato_multipleCultivars/s2_10_Cultivars_PacBio_HiFi/a1_HiFi_stats/

# done 
for sample in A B C D E F G H I J; do 
    cd ${wd}
    mkdir sample_${sample}
    cd sample_${sample}
    ccs=/biodata/dep_mercier/grp_schneeberger/reads/Potato/multi_cultivars_2021/3_10_Cultivars_PacBio_HiFi/sample_${sample}/*/*.gz
    bsub -q normal -o zcat.log -e zcat.err  "echo 0; zcat ${ccs} > sample_${sample}_zcat_HiFi.fastq"
done

# done
for sample in A B C D E F G H I J; do 
    #
    cd ${wd}
    #
    # quality checking: python3 necessay!
    bsub -q ioheavy -n 20 -R "span[hosts=1] rusage[mem=180000]" -M 180000 -o longQC.log -e longQC.err "/netscratch/dep_mercier/grp_schneeberger/bin/bin/anaconda/install/bin/python3.7 /netscratch/dep_mercier/grp_schneeberger/bin/LongQC/longQC.py sampleqc -x pb-rs2 -c sample_${sample}_zcat_HiFi_trim.fastq -o qc_hifi -p 20 sample_${sample}_zcat_HiFi.fastq; cat sample_${sample}_zcat_HiFi_trim.fastq | awk '{if(NR%4==2) print length(\$1)}' > sample_${sample}_zcat_HiFi_length_4cells.txt"
    #
    # length checking: raw data 
    ccs=sample_${sample}_zcat_HiFi.fastq
    bsub -q normal -R "span[hosts=1] rusage[mem=100]" -M 100 -o hifi_length.log -e hifi_length.err "echo 0; cat ${ccs} | awk '{if(NR%4==2) print length(\$1)}' > sample_${sample}_zcat_HiFi_length_4cells_raw.txt"
    #
done

# compress

for sample in A B C D E F G H I J; do 
    cd ${wd}
    cd sample_${sample}
    bsub -q normal -o gzip.log -e gzip.err  "gzip sample_${sample}_zcat_HiFi_trim.fastq"
    # nohup gzip sample_${sample}_zcat_HiFi_trim.fastq &>gzip.log&
done


# visualize read length: need manually update path etc in the R script 

Rscript sample_x_10ptt_HiFi_readsizes.R

