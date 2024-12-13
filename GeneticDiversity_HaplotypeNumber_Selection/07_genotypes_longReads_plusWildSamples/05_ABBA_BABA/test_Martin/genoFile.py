from cyvcf2 import VCF
import sys

# mamba activate mypython3
file="/dss/dsslegfs01/pn29fi/pn29fi-dss-0012/04_analyses/05_haplotypeBased_analyses/07_genotypes_longReads_plusWildSamples/04_genotypeCalling/merged_vcf/chr/chr12_filterQUAL30.vcf.gz"
output="test.geno"

file = sys.argv[1]
output = str(sys.argv[2])

vf = VCF(file, strict_gt=True)

header_samples = '\t'.join(map(str, vf.samples))
header = "#CHROM\tPOS\t" + header_samples

outputFile = open(f'{output}', "w")

outputFile.write(f'{header}\n')

for variant in vf:
    if len(variant.ALT)==1:
        alleles = list(variant.REF) + variant.ALT
        genotypes = []
        for sampleIndex in range(len(vf.samples)):
            gt = variant.genotypes[sampleIndex]
            gt_filtered = [alleles[s] if s != -1 else 'N' for s in gt[:-1]]
            genotypes.append('/'.join(map(str, gt_filtered)))
        line_genotypes = '\t'.join(map(str, genotypes))
        outputFile.write(f'{variant.CHROM}\t{variant.POS}\t{line_genotypes}\n')



outputFile.close()

