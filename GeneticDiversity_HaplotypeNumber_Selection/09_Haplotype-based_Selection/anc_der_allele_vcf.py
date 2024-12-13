import argparse
import random
from collections import defaultdict
from cyvcf2 import VCF, Writer

def count_alleles(variant, ancestral_samples, list_samples):
    allele_counts = defaultdict(lambda: defaultdict(int))
    
    # Count alleles in ancestral samples
    for sample in ancestral_samples:
        sample_index = list_samples.index(sample)
        genotype = variant.genotypes[sample_index][:2]  # Get the alleles
        for allele in genotype:
            if allele != -1:  # Ignore missing alleles
                allele_counts[sample][allele] += 1
    
    return allele_counts

def determine_most_common_allele(allele_counts):
    overall_counts = defaultdict(int)
    # Aggregate counts across all samples
    for sample_counts in allele_counts.values():
        for allele, count in sample_counts.items():
            overall_counts[allele] += count
    if not overall_counts:
        return None
    # Find most common allele
    most_common_allele = max(overall_counts, key=overall_counts.get)
    return most_common_allele


def count_alleles_cul(variant, output_samples, list_samples):
    allele_counts2 = defaultdict(lambda: defaultdict(int))
    
    # Count alleles in ancestral samples
    for sample in output_samples:
        sample_index = list_samples.index(sample)
        genotype = variant.genotypes[sample_index][:2]  # Get the alleles
        for allele in genotype:
            if allele != -1:  # Ignore missing alleles
                allele_counts2[sample][allele] += 1
    
    return allele_counts2


def determine_second_common_allele(allele_counts2, output_samples, most_common_allele):
    output_alleles = defaultdict(int)
    
    # Aggregate counts of alleles across output samples
    for sample in output_samples:
        for allele, count in allele_counts2[sample].items():
            if str(allele) != str(most_common_allele):
                output_alleles[allele] += count
    
    # Sort output alleles by count descending
    sorted_output_alleles = sorted(output_alleles.items(), key=lambda x: x[1], reverse=True)
    
    # Second most common allele is at index 1 (if exists)
    if len(sorted_output_alleles) > 0:
        second_common_allele, second_common_count = sorted_output_alleles[0]
        # Check if MAC is greater than 2
        if second_common_count > 2:
            return second_common_allele
        else:
            return None
    else:
        return None

def update_genotype(genotype, ancestral_allele, allele_map):
    new_genotype = []
    for allele in genotype:
        if str(allele) == str(ancestral_allele):
            new_genotype.append(0)
        else:
            new_genotype.append(allele_map.get(allele, allele))
    return new_genotype

def main(vcf_file, output_file, ancestral_samples, output_samples):
    vcf_reader = VCF(vcf_file)
    vcf_writer = Writer(output_file, vcf_reader)

    list_samples = vcf_reader.samples

    output_sample_indices = [list_samples.index(sample) for sample in output_samples]

    random.seed(42)  # Seed for reproducibility
    
    count = 0
    included_count = 0
    
    for record in vcf_reader:
        count += 1
        
        # Generate random number between 0 and 1
        random_number = random.random()
        
        # Include variant with 10% probability (random number < 0.1)
        if random_number < 0.1:
            included_count += 1
            
            # Count alleles in ancestral and output samples
            allele_counts = count_alleles(record, ancestral_samples, list_samples)

            
            # Determine most common allele among ancestral samples
            ancestral_allele = determine_most_common_allele(allele_counts)
            
            # Determine second most common allele among output samples
            allele_counts2 = count_alleles_cul(record, output_samples, list_samples)
            derived_allele = determine_second_common_allele(allele_counts2, output_samples, ancestral_allele)

            # Update genotypes
            if (ancestral_allele is not None) and (derived_allele is not None):

                # Map alternative alleles to new indices
                all_alleles = [record.REF] + record.ALT
                #ancestral_allele_index = all_alleles.index(sample)
                ancestral_allele_Base = all_alleles[int(ancestral_allele)]
                alt_alleles = list(set(all_alleles) - {ancestral_allele_Base})
                allele_map = {all_alleles.index(allele): i+1 for i, allele in enumerate(alt_alleles)}

                Allsample_indices = [list_samples.index(sample) for sample in list_samples]

                for i in Allsample_indices:
                    genotype = record.genotypes[i][:2]
                    new_genotype = update_genotype(genotype, ancestral_allele, allele_map)
                    record.genotypes[i] = new_genotype + [False]

                record.genotypes = record.genotypes


                record.REF = ancestral_allele_Base
                record.ALT = alt_alleles

                # Write the updated record to output VCF file
                vcf_writer.write_record(record)
    
    vcf_writer.close()
    print(f"Total variants processed: {count}")
    print(f"Variants included in output: {included_count}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process VCF file and modify genotypes.')
    parser.add_argument('vcf_file', metavar='vcf_file', type=str, help='Input VCF file (compressed)')
    parser.add_argument('output_file', metavar='output_file', type=str, help='Output VCF file')
    parser.add_argument('--ancestral_samples', metavar='ancestral_samples', type=str, nargs='+', help='List of ancestral sample names')
    parser.add_argument('--output_samples', metavar='output_samples', type=str, nargs='+', help='List of output sample names')
    
    args = parser.parse_args()
    
    ancestral_samples = args.ancestral_samples
    output_samples = args.output_samples
    
    main(args.vcf_file, args.output_file, ancestral_samples, output_samples)




