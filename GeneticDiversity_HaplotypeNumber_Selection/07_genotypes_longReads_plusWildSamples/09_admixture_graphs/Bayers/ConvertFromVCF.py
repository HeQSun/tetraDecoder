from cyvcf2 import VCF
import sys
import numpy as np
import gzip

#file = "sampleData.vcf"
#output = "adbayesinput.txt"
#sample_table_file = "samples.txt"

# Read the input arguments (uncomment these lines for actual usage)
file = sys.argv[1]
output = sys.argv[2]
sample_table_file = sys.argv[3]

# Load the sample names and populations from the provided file
sample_names = []
populations = []
with open(sample_table_file, 'r') as f:
    for line in f:
        parts = line.strip().split()
        sample_names.append(parts[0])
        populations.append(parts[1])

# Open the VCF file
vf = VCF(file, strict_gt=True)

# Get the indices of the samples to be included in the output
sample_indices = [vf.samples.index(sample) for sample in sample_names if sample in vf.samples]

# Ensure populations order matches the order of sample names in the VCF file
populations = [populations[sample_names.index(name)] for name in vf.samples if name in sample_names]

# Number of haplotypes and unique populations
NumberOfHaplotypes = 2 * len(sample_indices)
HaplotypePopulations = [pop for pop in populations for _ in (0, 1)]
UniquePops = sorted(set(populations))

# Prepare to write output
with open(output, 'w') as outputFile:
    # Write the populations as the header
    outputFile.write(" ".join(UniquePops) + "\n")
    # Process each variant in the VCF file
    for variant in vf:
        if len(variant.ALT) == 1:  # Check there is only one alternative allele
            # Prepare the genotype data matrix for this SNP
            Data = np.full(NumberOfHaplotypes, -1)
            for i, sampleIndex in enumerate(sample_indices):
                gt = variant.genotypes[sampleIndex]
                Data[2*i] = 0 if gt[0] == 0 else 1 if gt[0] == 1 else -1
                Data[2*i + 1] = 0 if gt[1] == 0 else 1 if gt[1] == 1 else -1
            # Match the data to each population
            output_line = []
            total_derived = 0
            total_ancestral = 0
            for pop in UniquePops:
                relevant_indices = [index for index, p in enumerate(HaplotypePopulations) if p == pop]
                nderived = np.sum(Data[relevant_indices] == 1)
                nancestral = np.sum(Data[relevant_indices] == 0)
                output_line.append(f"{nderived},{nancestral}")
                total_derived += nderived
                total_ancestral += nancestral
            #Include only polymorphic sites:
            if total_derived > 1 and total_ancestral > 1:
                # Write the output line for this SNP
                outputFile.write(" ".join(output_line) + "\n")



