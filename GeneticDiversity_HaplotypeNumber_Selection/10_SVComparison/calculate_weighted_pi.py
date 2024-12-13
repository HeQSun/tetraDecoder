import sys

def calculate_weighted_pi(input_file, output_file, window_size=1000000):
    windows = {}
    chromosomes = set()
    
    with open(input_file, 'r') as infile:
        header = infile.readline().strip().split('\t')
        for line in infile:
            columns = line.strip().split('\t')
            chrom = columns[0]
            window_start = int(columns[1])
            overlap_bp = int(columns[5])
            pi = float(columns[8])
            
            window_key = (chrom, window_start)
            if window_key not in windows:
                windows[window_key] = {'total_weighted_pi': 0, 'total_bp': 0}
            
            weight = overlap_bp / window_size
            weighted_pi = pi * weight
            
            windows[window_key]['total_weighted_pi'] += weighted_pi
            windows[window_key]['total_bp'] += overlap_bp
            chromosomes.add(chrom)

    # Ensure all windows are present
    for chrom in chromosomes:
        max_end = max(key[1] for key in windows if key[0] == chrom)
        for start in range(0, max_end, window_size):
            window_key = (chrom, start)
            if window_key not in windows:
                windows[window_key] = {'total_weighted_pi': 0, 'total_bp': 0}

    with open(output_file, 'w') as outfile:
        outfile.write('chr\tpos\tpi\n')
        for window_key in sorted(windows.keys()):
            chrom, window_start = window_key
            total_weighted_pi = windows[window_key]['total_weighted_pi']
            total_bp = windows[window_key]['total_bp']
            average_pi = total_weighted_pi if total_bp > 0 else 0
            outfile.write(f'{chrom}\t{window_start}\t{average_pi}\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python calculate_weighted_pi.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    calculate_weighted_pi(input_file, output_file)



