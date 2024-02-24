tetraDecoder
=
Here we introduce the pipeline for achieving haplotype-resolved, chromosome-level genome assembly of tetraploids using long read sequencing and HiC sequencing. This is to support the analysis in [Sun_and_Tusso_et_al_2024].

Pipeline of assembly
=

Befrore running the pipeline, install necessary customized tools, by

cd /your/path/to/install/tetraDecoder/

git clone https://github.com/HeQSun/tetraDecoder.git

then

cd src_gb_stage4_short_read_separation/gzlib

make clean

make

cd ../../

cd src_gb_stage4_short_read_separation_chr/gzlib

make clean

make

cd ../../

for si in src*/; do cd $si; make ; cd .. ; done

To start the process, check the [pipeline](https://github.com/HeQSun/tetraDecoder/tree/main/pipeline)

Note, in the pipeline, there are tools like hifiasm, bowtie2 etc required, which are all publicly available. Please install following their instructions.
