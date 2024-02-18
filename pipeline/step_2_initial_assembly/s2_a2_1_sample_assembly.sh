# Here we assembly the genomes: result of hifiasm v0.7 used.

wd=/working/path/
cd ${wd}

########################################################################################################################
#
# b.1.1 Assembly with hifiasm
#       currently no purging
#       homozygous coverage: coverage threshold for homozygous reads. Hifiasm prints it as: [M::purge_dups] homozygous read coverage threshold: X. If it is not around homozygous coverage, the final assembly might be either too large or too small. To fix this issue, please set --hom-cov to homozygous coverage.
#
#
wd=/working/path/
cd ${wd}
#
for sample in A B C D E F G H I J; do 
    wd=/working/path/
    cd ${wd}
    mkdir sample_${sample}
    #
    cd ${wd}/sample_${sample}
    mkdir hifiasm_asm_v0p7
    cd hifiasm_asm_v0p7
    #
    ll ../../../a1_HiFi_stats/sample_${sample}/sample_${sample}_zcat_HiFi_trim.fastq.gz
    ln -s ../../../a1_HiFi_stats/sample_${sample}/sample_${sample}_zcat_HiFi_trim.fastq.gz sample_${sample}.fastq.gz
    # this higher version of hifiasm has many redundant contigs
    # ln -s /netscratch/dep_mercier/grp_schneeberger/bin/hifiasm/hifiasm-0.16/hifiasm hifiasm
    #
    hifi=sample_${sample}.fastq.gz
    # discarded - as v0.16 gave many redundant unitigs and low N50!
    # hifiasm v0p7, with fewer redundant unitigs and higher N50 (than version 0.16) - finally used 20220921
    bsub -q ioheavy -R "rusage[mem=200000]" -M 200000 -n 40 -o ${sample}_hifiasm.log -e ${sample}_hifiasm.err "./hifiasm -t 40 -o ${sample}_hifiasm ${hifi}"    
    # 
done

# b.1.2 Get assembly stats
#
wd=/working/path/
cd ${wd}
for sample in A B D E F G H I J; do 
    cd ${wd}/sample_${sample}
    cd hifiasm_asm_v0p7
    ####################################################################################################################
    # outprefix=${sample}_hifiasm.bp.p_utg # for hifiasm v0.16
    outprefix=${sample}_hifiasm.p_utg    # for hifiasm v0.7
    # utg-level
    cat ${outprefix}.gfa | grep '^S' | cut -f2,3 | awk '{print ">"$1"\n"$2}' > ${outprefix}.gfa.fasta
    # total:  bp
    fasta_length ${outprefix}.gfa.fasta | grep '>' | sed 's/>//g' | sort -k2,2n > ${outprefix}.gfa.ctgsizes
    awk '{s+=$2} END {print s}' ${outprefix}.gfa.ctgsizes
    # links between contigs and their overlaps
    cat ${outprefix}.gfa | grep '^L' > ${outprefix}.gfa.links
    #
    # N50
    /path/to/GitHub-schneeberger-group-jac/toolbox/Analysis/Assembly/calc_CN50.pl ${outprefix}.gfa.fasta 844000000 1 &> ${outprefix}.gfa.N50.calc.result.txt&
    #
    ####################################################################################################################
#
done
#
########################################################################################################################
#
# b.2.1 (optional) Assembly with canu 2.2: corOutCoverage=200 errorRate=0.013 "batOptions=-dg 3 -db 3 -dr 1 -ca 500 -cp 50"
#       /netscratch/dep_mercier/grp_schneeberger/bin/canu/canu-2.2/bin/canu
#
wd=/working/path/
cd ${wd}
sample="C"
tool="canu"
cd ${wd}/sample_${sample}
mkdir ${tool}_asm
cd ${tool}_asm
#
ln -s ../../../a1_HiFi_stats/sample_${sample}/sample_${sample}_zcat_HiFi_trim.fastq sample_${sample}.fastq
hifi=sample_${sample}.fastq
bsub -q ioheavy -R "rusage[mem=120000]" -M 120000 -n 20 -o ${sample}_${tool}.log -e ${sample}_${tool}.err "${tool} minInputCoverage=2 -p ${sample}_${tool} -d ${sample}_${tool} -pacbio-hifi ${hifi} useGrid=false genomeSize=840m corMinCoverage=1 corOutCoverage=30 minOverlapLength=1000 correctedErrorRate=0.03 executiveThreads=20 > ${tool}_${sample}.log"
#
/path/to/GitHub-schneeberger-group-jac/toolbox/Analysis/Assembly/calc_CN50.pl C_canu.contigs.fasta 844000000 1 &> C_canu.bp.p_utg.gfa.N50.calc.result.txt&
