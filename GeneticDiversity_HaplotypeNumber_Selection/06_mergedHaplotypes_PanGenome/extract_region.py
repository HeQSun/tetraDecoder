from Bio import SeqIO
import sys

def extract_region(reference_file, chromosome, start, end, output_file, sequence_name):
    """
    Extracts a region from the reference genome and saves it to a new FASTA file.
    """
    sequences = SeqIO.to_dict(SeqIO.parse(reference_file, "fasta"))

    if chromosome not in sequences:
        print(f"Chromosome {chromosome} not found in the reference file.")
        sys.exit(1)

    reference_sequence = sequences[chromosome]
    region_sequence = reference_sequence[start-1:end]
    region_sequence.id = sequence_name  # Set the sequence name
    region_sequence.description = ""  # Clear the description to have only the sequence name in the header

    with open(output_file, "w") as output_handle:
        SeqIO.write(region_sequence, output_handle, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 7:
        print("Usage: python extract_region.py <reference_file> <chromosome> <start> <end> <output_file> <sequence_name>")
        sys.exit(1)

    reference_file = sys.argv[1]
    chromosome = sys.argv[2]
    start = int(sys.argv[3])
    end = int(sys.argv[4])
    output_file = sys.argv[5]
    sequence_name = sys.argv[6]

    extract_region(reference_file, chromosome, start, end, output_file, sequence_name)




