tetraDecoder
=
Here we introduce the pipeline for achieving haplotype-resolved, chromosome-level genome assembly of tetraploids using long read sequencing and HiC sequencing. This is to support the analysis in [Sun_and_Tusso_et_al_2024].

Pipeline of assembly
=

Befrore running the pipeline, install necessary customized tools, by

cd /your/path/to/install/tetraDecoder/

git clone https://USER-GITHUB-ACCOUNT:USER-TOKEN@github.com/HeQSun/tetraDecoder

or

git clone https://github.com/HeQSun/tetraDecoder.git

where SER-GITHUB-ACCOUNT and USER-TOKEN need to be updated according to your own settings.

then

for si in src*/; do cd $si; make ; cd .. ; done

To start the process, check the [pipeline](https://github.com/HeQSun/tetraDecoder/tree/main/pipeline)

Note, in the pipeline, there are tools like hifiasm, bowtie2 etc required, which are all publicly available. Please install following their instructions.
