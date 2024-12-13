import pandas as pd
import numpy as np
import sys

file = str(sys.argv[1])
output = str(sys.argv[2])
winSize = int(sys.argv[3])

# Read the data from the text file
df1 = pd.read_csv(file, sep='\t')

# Define window size
window_size = winSize

# Initialize an empty list for recombination rates per window
window_data = []

# Iterate over each chromosome
for chrom, chrom_group in df1.groupby('Chrom'):
    # Calculate recombination rate (centimorgan/Mbp) for each interval within the chromosome
    chrom_group['RecombinationRate'] = chrom_group['MapPosition'].diff() / (chrom_group['Position'].diff() / 1e6)
    # Remove the first row of each chromosome because it will have NaN for recombination rate
    chrom_group = chrom_group.dropna(subset=['RecombinationRate'])
    # Create a function to map position to window
    def position_to_window(position, window_size):
        return int((position // window_size) + 1)
    # Iterate over each row in the chromosome group to fill recombination rates per window
    for i in range(1, len(chrom_group)):
        start_pos = chrom_group.iloc[i-1]['Position']
        end_pos = chrom_group.iloc[i]['Position']
        rate = chrom_group.iloc[i]['RecombinationRate']
        start_window = position_to_window(start_pos, window_size)
        end_window = position_to_window(end_pos, window_size)
        for window in range(start_window, end_window + 1):
            window_start_pos = (window - 1) * window_size
            window_end_pos = window * window_size
            overlap_start = max(start_pos, window_start_pos)
            overlap_end = min(end_pos, window_end_pos)
            # Calculate the proportion of the interval that falls within the current window
            proportion = (overlap_end - overlap_start) / (end_pos - start_pos)
            #window_rate = rate * proportion
            size_bp = (end_pos - start_pos) * proportion
            window_data.append([int(chrom), window, size_bp, rate])

# Convert the list to a dataframe
recombination_per_window = pd.DataFrame(window_data, columns=['Chrom', 'Window', 'size_bp', 'RecombinationRate'])

recombination_per_window['GenDistance'] = (recombination_per_window['size_bp']*recombination_per_window['RecombinationRate'])/winSize

# Group by Chrom and Window and calculate the mean recombination rate per window
recombination_per_window = recombination_per_window[['Chrom', 'Window', 'GenDistance']].groupby(['Chrom', 'Window'], as_index=False).sum()

# Ensure the 'Chrom' column is of integer type
recombination_per_window['Chrom'] = recombination_per_window['Chrom'].astype(int)

# Save the output to a new text file
recombination_per_window.to_csv(output, sep='\t', index=False)



