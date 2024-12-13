import numpy as np
import sys

def filter_distance_matrix(header, matrix, cultivar_samples, wild_samples):
    """
    Filters the distance matrix to include only cultivar samples in rows and wild samples in columns.

    Parameters:
    header (list): A list of sample names corresponding to the rows and columns of the matrix.
    matrix (numpy.ndarray): A 2D numpy array representing the distance matrix.
    cultivar_samples (list): A list of sample names corresponding to cultivar samples.
    wild_samples (list): A list of sample names corresponding to wild samples.

    Returns:
    tuple: A tuple containing the list of cultivar samples in the correct order, the list of wild samples in the correct order, and the filtered distance matrix.
    """
    
    # Get the indices of the cultivar samples
    cultivar_indices = [header.index(sample) for sample in cultivar_samples if sample in header]
    ordered_cultivar_samples = [header[i] for i in cultivar_indices]
    
    # Get the indices of the wild samples
    wild_indices = [header.index(sample) for sample in wild_samples if sample in header]
    ordered_wild_samples = [header[i] for i in wild_indices]
    
    # Filter the matrix
    filtered_matrix = matrix[np.ix_(cultivar_indices, wild_indices)]
    
    return ordered_cultivar_samples, ordered_wild_samples, filtered_matrix



def find_lowest_distances(distance_matrix, cultivar_samples, wild_samples):
    """
    Finds the lowest distances between cultivar samples and wild samples in the distance matrix.
    Parameters:
    distance_matrix (numpy.ndarray): The distance matrix.
    cultivar_samples (list): List of cultivar samples.
    wild_samples (list): List of wild samples.
    Returns:
    list: A list of tuples containing the lowest distances, cultivar samples, and wild samples.
    """
    lowest_distances = []
    # Iterate over each row in the distance matrix
    for i, cultivar_sample in enumerate(cultivar_samples):
        # Find the indices of the lowest distances in the current row
        sorted_indices = np.argsort(distance_matrix[i])  # Sort indices based on distances
        lowest_indices = sorted_indices[:3]  # Get indices of lowest 3 distances
        # Store the lowest distances and corresponding samples
        lowest_dist_values = [distance_matrix[i, idx] for idx in lowest_indices]
        max_dist = max(lowest_dist_values)
        NTopSamples = sum(distance_matrix[i]<=max_dist)
        # Handle multiple samples with the same top value
        if NTopSamples < 10:
            for idx, dist in enumerate(distance_matrix[i]):
                if dist <= max_dist:
                    lowest_distances.append([cultivar_sample, wild_samples[idx], dist])
        else:
            lowest_distances.append([cultivar_sample, "NA", "NA"])
    return lowest_distances



def read_samples(file_path):
    """
    Reads cultivar samples from a text file and returns a list of these samples.
    Parameters:
    file_path (str): The path to the text file containing the cultivar samples.
    Returns:
    list: A list of cultivar samples.
    """
    with open(file_path, 'r') as file:
        # Read all lines from the file and strip any leading/trailing whitespace characters
        list_samples = [line.strip() for line in file.readlines()]
    return list_samples




# Processing file:

filename="chr01_filterQUAL30.dis"
win_information="ed_chr01_filterQUAL30.winData"
cultivar_samples_file="list_cultivars.txt"
wild_samples_file="list_wild.txt"
output_file_name="test.txt"


filename=str(sys.argv[1])
win_information=str(sys.argv[2])
cultivar_samples_file=str(sys.argv[3])
wild_samples_file=str(sys.argv[4])
output_file_name=str(sys.argv[5])

cultivar_samples = read_samples(cultivar_samples_file)
wild_samples = read_samples(wild_samples_file)
Win_info = read_samples(win_information)
output_file = open(output_file_name, "w")



win=0
for line in open(filename, "r"):
    parts = line.strip().split()
    if len(parts) > 1:
        header.append(parts[0])
        matrix.append([float(x) for x in parts[1:]])
    elif win==0:
        header = []
        matrix = []
        win+=1
    else:
        ordered_cultivar_samples, ordered_wild_samples, filtered_matrix = filter_distance_matrix(header, np.array(matrix), cultivar_samples, wild_samples)
        lowest_distances_array = find_lowest_distances(filtered_matrix, ordered_cultivar_samples, ordered_wild_samples)
        for topvalue in lowest_distances_array:
            if topvalue[2] != "NA":
                values="\t".join(map(str, topvalue))
                output_file.write(f'{Win_info[win-1]}\t{values}\n')
        win+=1
        header = []
        matrix = []

output_file.close()


