import sys

def process_sv_table(input_file, output_file):
    window_size = 1000000  # 1 Mb

    def calculate_pi(sample_data):
        # Calculate genetic diversity (pi)
        n = len(sample_data)
        if n < 2:
            return 0.0
        
        pairwise_differences = 0
        for i in range(n):
            for j in range(i + 1, n):
                if sample_data[i] != sample_data[j]:
                    pairwise_differences += 1
        
        num_comparisons = n * (n - 1) / 2
        pi = pairwise_differences / num_comparisons if num_comparisons > 0 else 0.0
        return pi

    with open(input_file, 'r') as infile:
        header = infile.readline().strip().split('\t')
        samples = header[4:]  # Get sample names
        rows = [line.strip().split('\t') for line in infile]

    # Convert columns to appropriate types
    processed_rows = []
    for row in rows:
        if len(row) < 4:
            print(f"Skipping line due to insufficient columns: {row}")
            continue
        try:
            chrom = row[0]
            start = int(row[1])
            end = int(row[2])
            ann = row[3]
            sample_data = row[4:]
            processed_rows.append([chrom, start, end, ann] + sample_data)
        except ValueError as e:
            print(f"Skipping line due to conversion error: {row}, error: {e}")
            continue

    # Calculate windows and overlap
    results = []
    for row in processed_rows:
        chrom, start, end, ann, *sample_data = row
        window_start = (start // window_size) * window_size
        window_end = window_start + window_size

        while window_start <= end:
            overlap_start = max(start, window_start)
            overlap_end = min(end, window_end)
            overlap_bp = max(0, overlap_end - overlap_start + 1)

            if overlap_bp > 0:
                num_samples_with_1 = sum(1 for x in sample_data if x == '1')
                num_samples_with_0 = sum(1 for x in sample_data if x == '0')
                pi = calculate_pi(sample_data)
                results.append([chrom, window_start, window_end, start, end, overlap_bp, num_samples_with_1, num_samples_with_0, pi])

            window_start += window_size
            window_end += window_size

    # Write results to output file
    with open(output_file, 'w') as outfile:
        outfile.write('CHR\tWINDOW_START\tWINDOW_END\tSTART\tEND\tOVERLAP_BP\tNUM_SAMPLES_WITH_1\tNUM_SAMPLES_WITH_0\tPI\n')
        for result in results:
            outfile.write('\t'.join(map(str, result)) + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python process_sv_windows.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    process_sv_table(input_file, output_file)



