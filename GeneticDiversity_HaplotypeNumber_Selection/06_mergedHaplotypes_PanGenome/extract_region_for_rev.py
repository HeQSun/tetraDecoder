from Bio import SeqIO
from Bio.Seq import Seq
import sys

def extract_region(reference_file, chromosome, start, end, output_file_forward, output_file_reverse, sequence_name):
    """
    Extracts a region from the reference genome and saves it to two FASTA files (forward and reverse complement).
    """
    sequences = SeqIO.to_dict(SeqIO.parse(reference_file, "fasta"))

    if chromosome not in sequences:
        print(f"Chromosome {chromosome} not found in the reference file.")
        sys.exit(1)

    reference_sequence = sequences[chromosome]
    region_sequence = reference_sequence[start-1:end]
    region_sequence.id = sequence_name  # Set the sequence name
    region_sequence.description = ""  # Clear the description to have only the sequence name in the header

    # Save forward strand to the first output file
    with open(output_file_forward, "w") as output_handle_forward:
        SeqIO.write(region_sequence, output_handle_forward, "fasta")

    # Save reverse complement to the second output file
    reverse_complement_sequence = region_sequence.reverse_complement()
    reverse_complement_sequence.id = sequence_name #+ "_reverse_complement"
    reverse_complement_sequence.description = ""  # Clear the description to have only the sequence name in the header
    with open(output_file_reverse, "w") as output_handle_reverse:
        SeqIO.write(reverse_complement_sequence, output_handle_reverse, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 8:
        print("Usage: python extract_region.py <reference_file> <chromosome> <start> <end> <output_file_forward> <output_file_reverse> <sequence_name>")
        sys.exit(1)

    reference_file = sys.argv[1]
    chromosome = sys.argv[2]
    start = int(sys.argv[3])
    end = int(sys.argv[4])
    output_file_forward = sys.argv[5]
    output_file_reverse = sys.argv[6]
    sequence_name = sys.argv[7]

    extract_region(reference_file, chromosome, start, end, output_file_forward, output_file_reverse, sequence_name)




