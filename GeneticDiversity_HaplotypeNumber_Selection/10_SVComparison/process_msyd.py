import sys

def process_sv_table(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        header = infile.readline()
        outfile.write(header)

        for line in infile:
            columns = line.strip().split('\t')
            new_columns = []

            # Process the columns
            for i, col in enumerate(columns):
                if i < 4:  # Assuming the first four columns are chromosome, start, end, and annotation
                    new_columns.append(col)
                else:
                    if col == "-":
                        new_columns.append('0')
                    else:
                        try:
                            # If it can be converted to a number, keep it as is
                            float(col)
                            new_columns.append(col)
                        except ValueError:
                            # Otherwise, replace with 1
                            new_columns.append('1')
            outfile.write('\t'.join(new_columns) + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python process_sv.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    process_sv_table(input_file, output_file)



