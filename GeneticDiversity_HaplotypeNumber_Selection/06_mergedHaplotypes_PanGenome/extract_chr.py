from Bio import SeqIO
import sys

def extract_chromosome(reference_file, chromosome, output_file):
    """
    Extracts the entire sequence of a specified chromosome from the reference genome
    and saves it to a new FASTA file.
    """
    sequences = SeqIO.to_dict(SeqIO.parse(reference_file, "fasta"))

    if chromosome not in sequences:
        print(f"Chromosome {chromosome} not found in the reference file.")
        sys.exit(1)

    chromosome_sequence = sequences[chromosome]

    with open(output_file, "w") as output_handle:
        SeqIO.write(chromosome_sequence, output_handle, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python extract_chromosome.py <reference_file> <chromosome> <output_file>")
        sys.exit(1)

    reference_file = sys.argv[1]
    chromosome = sys.argv[2]
    output_file = sys.argv[3]

    extract_chromosome(reference_file, chromosome, output_file)




