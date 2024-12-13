import sys
import gzip
import random
from collections import defaultdict

def calculate_artificial_heterozygosity(vcf_file, output_file, target_samples):
    # Parse the target samples list
    target_samples_list = target_samples.split(',')

    # Initialize data structures
    windows = defaultdict(lambda: defaultdict(list))
    window_size = 1_000_000  # Using underscore for better readability

    # Create 20 random pairs from the target samples
    random.shuffle(target_samples_list)
    paired_samples = [target_samples_list[i:i + 2] for i in range(0, len(target_samples_list), 2)]

    if len(paired_samples) != 20:
        print("Error: Number of pairs is not equal to 20. Check the input sample list.")
        sys.exit(1)

    # Read the VCF file
    with gzip.open(vcf_file, 'rt') if vcf_file.endswith('.gz') else open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#'):
                headers = line.strip().split('\t')
                samples = headers[9:]
                sample_indices = {sample: idx for idx, sample in enumerate(samples) if sample in target_samples_list}
                continue
            
            cols = line.strip().split('\t')
            chrom = cols[0]
            pos = int(cols[1])
            genotypes = cols[9:]
            
            window_start = (pos // window_size) * window_size

            # Create artificial diploids for each pair
            for pair_idx, (sample1, sample2) in enumerate(paired_samples):
                idx1 = sample_indices[sample1]
                idx2 = sample_indices[sample2]
                
                gt1 = genotypes[idx1].split(':')[0]
                gt2 = genotypes[idx2].split(':')[0]

                if gt1 == './.' or gt2 == './.':
                    artificial_gt = None
                else:
                    alleles1 = gt1.replace('|', '/').split('/')
                    alleles2 = gt2.replace('|', '/').split('/')
                    allele1 = random.choice(alleles1)
                    allele2 = random.choice(alleles2)
                    artificial_gt = f"{allele1}/{allele2}"

                if artificial_gt == '0/1' or artificial_gt == '1/0':
                    windows[(chrom, window_start)][f"Artificial_{pair_idx+1}"].append(1)
                elif artificial_gt is None:
                    windows[(chrom, window_start)][f"Artificial_{pair_idx+1}"].append(None)
                else:
                    windows[(chrom, window_start)][f"Artificial_{pair_idx+1}"].append(0)

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
        print("Usage: python Heterozygosity.py <input_vcf> <output_csv> <target_samples>")
        print("Target samples should be a comma-separated list of sample names to include")
        sys.exit(1)
    
    input_vcf = sys.argv[1]
    output_csv = sys.argv[2]
    target_samples = sys.argv[3]
    calculate_artificial_heterozygosity(input_vcf, output_csv, target_samples)




