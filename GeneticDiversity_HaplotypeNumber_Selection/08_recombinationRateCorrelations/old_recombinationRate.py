import pandas as pd
import sys

file = str(sys.argv[1])
output = str(sys.argv[2])
winSize = int(sys.argv[3])

# Read the data from the text file
df1 = pd.read_csv(file, sep='\t')

# Define window size (500 kb)
window_size = winSize

# Add window column based on the position
df1['Start Window'] = (df1['Position'] // window_size) + 1
df1['End Window'] = df1['Start Window'].shift(-1, fill_value=df1['Start Window'].iloc[-1])

# Calculate recombination rate (centimorgan/Mbp) for each interval
df1['RecombinationRate'] = df1['MapPosition'].diff() / (df1['Position'].diff() / 1e6)

# Initialize an empty dataframe for recombination rates per window
recombination_per_window = pd.DataFrame(columns=['Chrom', 'Window', 'RecombinationRate'])

# Iterate over each row in the first dataframe to fill recombination rates per window
for i, row in df1.iterrows():
    if i < len(df1) - 1:  # Skip the last row to prevent index out of range
        start_window = int(row['Start Window'])
        end_window = int(df1.loc[i + 1, 'Start Window'])
        recombination_rate = row['RecombinationRate']
        chrom = row['Chrom']
        for window in range(start_window, end_window):
            recombination_per_window = recombination_per_window.append(
                {'Chrom': chrom, 'Window': window, 'RecombinationRate': recombination_rate}, ignore_index=True
            )

# Exclude NaN values before calculating the mean recombination rate per window
recombination_per_window.dropna(subset=['RecombinationRate'], inplace=True)

# Group by Chrom and Window and calculate the mean recombination rate per window
recombination_per_window = recombination_per_window.groupby(['Chrom', 'Window'], as_index=False).mean()

# Save the output to a new text file
recombination_per_window.to_csv(output, sep='\t', index=False)



