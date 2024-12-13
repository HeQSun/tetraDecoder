import os
import argparse
import random
import gzip

# Function to process a VCF line and extract relevant information
def process_vcf_line(line, header):
    columns = line.strip().split()
    chr_, pos, _, ref, alt = columns[:5]
    genotype_info = columns[9:]  # Skip the first 9 columns to get genotype information
    return chr_, pos, ref, alt, genotype_info

# Function to randomly sample one allele from a diploid genotype
def sample_allele(genotype):
    alleles = genotype.split('|')
    if len(alleles) == 1:
        alleles = genotype.split('/')
    return random.choice(alleles) if len(alleles) > 1 else alleles[0]

def calculate_allele_counts(alleles):
    allele_counts = {}
    for allele in alleles:
        if allele != '.' and allele != 'nan':
            if allele in allele_counts:
                allele_counts[allele] += 1
            else:
                allele_counts[allele] = 1
    return allele_counts

def get_minor_allele_count(allele_counts):
    if len(allele_counts) < 2:
        return 0
    sorted_counts = sorted(allele_counts.values(), reverse=True)
    return sorted_counts[1] if len(sorted_counts) > 1 else 0

def main(input_file, output_dir, final_output):
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Read VCF header and initialize files
    with open(input_file, 'rt') if input_file.endswith('.vcf') else gzip.open(input_file, 'rt') as f:
        for line in f:
            if line.startswith('#CHROM'):
                header = line.strip().split()[9:]  # Sample names start from the 10th column
                break

    # Create a file for each sample
    files = {sample: open(os.path.join(output_dir, f"{sample}.txt"), 'w') for sample in header}

    map_output = final_output.replace(".hap", ".map")
    with open(map_output, 'w') as map_file:
        # Process the input VCF file line by line
        with open(input_file, 'rt') if input_file.endswith('.vcf') else gzip.open(input_file, 'rt') as f:
            variant_id = 1
            for line in f:
                if line.startswith('#'):
                    continue
                
                chr_, pos, ref, alt, genotype_info = process_vcf_line(line, header)
                
                # Randomly sample one allele per sample
                alleles = [sample_allele(genotype.split(':')[0]) for genotype in genotype_info]

                # Calculate allele counts and MAC
                allele_counts = calculate_allele_counts(alleles)
                mac = get_minor_allele_count(allele_counts)

                # Check if number of unique alleles does not exceed 3 and MAC > 2
                if len(allele_counts) <= 3 and mac > 2:
                    # Convert alleles to numeric representation
                    alt_alleles = alt.split(',')
                    numeric_alt = ','.join(str(i + 1) for i in range(len(alt_alleles)))

                    # Write to map file
                    map_file.write(f"rs{variant_id} {chr_} {pos} 0 {numeric_alt}\n")
                    variant_id += 1

                    # Write to intermediate files
                    for sample, allele in zip(header, alleles):
                        if allele == '.' or allele == 'nan':
                            allele = '.'  # Replace 'nan' with '.'
                        files[sample].write(f"{chr_}_{pos} {allele}\n")

    # Close all intermediate files
    for f in files.values():
        f.close()

    # Combine intermediate files into the final output
    with open(final_output, 'w') as fout:
        for sample in header:
            sample_file = os.path.join(output_dir, f"{sample}.txt")
            with open(sample_file, 'r') as fin:
                lines = fin.readlines()
                fout.write(f"{sample} " + " ".join([line.split()[1] for line in lines]) + "\n")

    # Remove intermediate files
    for sample in header:
        os.remove(os.path.join(output_dir, f"{sample}.txt"))

    print(f"Data saved to {final_output}")
    print(f"Map data saved to {map_output}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process a large VCF file into a specified format.')
    parser.add_argument('--input_file', required=True, help='Path to the input VCF file (can be .vcf or .vcf.gz).')
    parser.add_argument('--output_dir', required=True, help='Directory to store intermediate files.')
    parser.add_argument('--final_output', required=True, help='Path to the final output file.')

    args = parser.parse_args()

    main(args.input_file, args.output_dir, args.final_output)

