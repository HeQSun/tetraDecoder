#!/bin/bash

#usage
if [ "$#" -le 3 ]
  then
echo "gmap.sh <genome> <ref_cds> <chr_gff> <nthread>"
exit 1
fi

genome=$1
ref_cds=$2
chr_gff=$3
nthread=$4

# step 1: ~5 mins
gmap_build -D . -d DB ${genome} >gmap_build.log
# step 2: note: ~30 mins
#         1) target.genome is the contig level of polyploid genome assembly
#         2) $N could be the ploidy of your target genome, for example $N=4 if it is a tetraploid
#         3) ref_cds is the coding sequences of diploid genome, which could be reference to get allelic table   

N=4
gmap -D . -d DB -t ${nthread} -f 2 -n $N ${ref_cds} > gmap.gff3
sed 's/No/\nNo/g' gmap.gff3 | grep -v 'No' > gmap2.gff3 # unexpected cases
mv gmap.gff3 gmap_raw.gff3
mv gmap2.gff3 gmap.gff3
sed -i 's/DM_gene_models_cdna_chr/DM_models_cdna_chr/g' gmap.gff3 # flag "gene" is hardcoded in gmap2AlleleTableBED2.pl by default
# step 3: note: 1 second
#         1) ref_gff3 is the gff3 annotation of the diploid genome
#         2) gmap2Alleletable.pl could be get in the following link: https://github.com/tangerzhang/ALLHiC/blob/master/scripts/gmap2AlleleTable.pl   
#         3) format of Allele.ctg.table: chromosomeID position contig1 [contig2 contig 3 contig4]. 
#                                        Prune step will remove the Hi-C linked reads between allelic contigs.
cut -f1,4,5,9 ${chr_gff} | sed 's/=/\t/g' | sed 's/;/\t/g' | cut -f 1,2,3,5 > ref.bed
# caution: I updated gmap2AlleleTableBED2 to get the expected info: chr pos ctg1[ ctg2 ctg3 ctg4]
# gmap.gff3 generated above is a hardcoded input of gmap2AlleleTableBED2.pl
perl /path/to/ALLHiC/scripts/gmap2AlleleTableBED2.pl ref.bed >gmap2AlleleTableBED2.log
#
