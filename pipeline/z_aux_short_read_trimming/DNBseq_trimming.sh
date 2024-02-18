# Here we clean DNB reads


for sample in R T; do
#
cd /your/path/to/collect/clean/reads/
mkdir clean_${sample}
cd clean_${sample}
#
readpath=/path/to/BGI_reads/sample_${sample}    
ln -s ${readpath}/${sample}_BGI_DNB_R1.fq.gz ln_${sample}_BGI_DNB_R1.fq.gz
ln -s ${readpath}/${sample}_BGI_DNB_R2.fq.gz ln_${sample}_BGI_DNB_R2.fq.gz
R1=ln_${sample}_BGI_DNB_R1.fq.gz
R2=ln_${sample}_BGI_DNB_R2.fq.gz
# tool here: /netscratch/dep_mercier/grp_schneeberger/bin/SOAPnuke/SOAPnuke
bsub -q ioheavy -R "rusage[mem=1000]" -M 2000 -n 4 -o ${sample}_SOAPnuke.log -e ${sample}_SOAPnuke.err "SOAPnuke filter -1 ${R1} -2 ${R2} -T 4 -l 20 -q 0.2 -n 0.001 -4 100 -f AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA -r AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG -C clean_${sample}_BGI_DNB_R1.fq.gz -D clean_${sample}_BGI_DNB_R2.fq.gz -o ./"
#
done
