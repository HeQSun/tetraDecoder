import sys
import gzip
from collections import defaultdict

def calculate_heterozygosity(vcf_file, output_file, exclude_samples):
    # Parse the exclude samples list
    exclude_samples_set = set(exclude_samples.split(','))

    # Initialize data structures
    windows = defaultdict(lambda: defaultdict(list))
    samples = []
    window_size = 1_000_000  # Using underscore for better readability

    # Read the VCF file
    with gzip.open(vcf_file, 'rt') if vcf_file.endswith('.gz') else open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#'):
                headers = line.strip().split('\t')
                samples = headers[9:]
                continue
            
            cols = line.strip().split('\t')
            chrom = cols[0]
            pos = int(cols[1])
            genotypes = cols[9:]
            
            window_start = (pos // window_size) * window_size
            
            for i, genotype in enumerate(genotypes):
                sample = samples[i]
                if sample in exclude_samples_set:
                    continue
                gt = genotype.split(':')[0]
                gt = gt.replace('|', '/')  # Replace "|" with "/" for uniformity
                if gt in ['0/1', '1/0']:
                    windows[(chrom, window_start)][sample].append(1)
                elif gt == './.':
                    windows[(chrom, window_start)][sample].append(None)
                else:
                    windows[(chrom, window_start)][sample].append(0)

    # Calculate heterozygosity and write to output file
    with open(output_file, 'w') as out:
        out.write('Sample,Chromosome,Window_Start,Mean_Heterozygosity\n')
        for (chrom, window_start), sample_data in windows.items():
            for sample, genotypes in sample_data.items():
                total_variants = len(genotypes)
                missing_data = genotypes.count(None)
                heterozygotes = sum(1 for gt in genotypes if gt == 1)
                if total_variants - missing_data > 0:
                    mean_heterozygosity = heterozygotes / (total_variants - missing_data)
                else:
                    mean_heterozygosity = 0
                out.write(f'{sample},{chrom},{window_start},{mean_heterozygosity}\n')

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python Heterozygosity.py <input_vcf> <output_csv> <exclude_samples>")
        print("Exclude samples should be a comma-separated list of sample names to exclude")
        sys.exit(1)
    
    input_vcf = sys.argv[1]
    output_csv = sys.argv[2]
    exclude_samples = sys.argv[3]
    calculate_heterozygosity(input_vcf, output_csv, exclude_samples)




