# ==============================================
# Section 1: DNA to Protein Translation
# ==============================================

# Standard genetic code (codon to amino acid mapping)

CODON_TABLE = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W'
}
# Function to translate DNA sequence to protein sequence

def translate_dna_to_protein(dna_sequence):
    # Step 1: Convert DNA to uppercase
    dna_sequence = dna_sequence.upper()
    
    # Step 2: Split into codons
    codons = [dna_sequence[i:i+3] for i in range(0, len(dna_sequence), 3)]
    
    # Step 3: Translate codons to amino acids
    protein_sequence = ''
    for codon in codons:
        if len(codon) == 3:  # Ensure it's a complete codon
            amino_acid = CODON_TABLE.get(codon, '?')  # Use '?' for unknown codons
            if amino_acid == '_':  # Stop codon
                break
            protein_sequence += amino_acid
    
    return protein_sequence

# Example DNA sequence
dna_sequence = "ATGGCCATTGTAATGGGCCGCTGAA"

# Translate to protein
protein_sequence = translate_dna_to_protein(dna_sequence)
print(f"DNA: {dna_sequence}")
print(f"Protein: {protein_sequence}")

# ============================================== 
# Initialiation of required libraries
# ==============================================

import random
import pandas as pd # type: ignore
import matplotlib.pyplot as plt # type: ignore

# ==============================================
# Section 2: Logistic Population Growth
# ==============================================

# Function to simulate logistic growth with randomized lag and exponential phases

def logistic_growth(K, P0, r, total_time, lag_mean, lag_std, exp_mean, exp_std):
    """
    Simulates logistic population growth with randomized lag and exponential phases.

    Parameters:
        K (float): Carrying capacity.
        P0 (float): Initial population size.
        r (float): Growth rate.
        total_time (int): Total simulation time.
        lag_mean (float): Mean duration of the lag phase.
        lag_std (float): Standard deviation of the lag phase.
        exp_mean (float): Mean duration of the exponential phase.
        exp_std (float): Standard deviation of the exponential phase.

    Returns:
        time (list): List of time points.
        population (list): List of population sizes (OD) at each time point.
    """
    # Step 1: Randomize lag and exponential phases
    lag_duration = max(0, int(random.gauss(lag_mean, lag_std)))  # Ensure lag is non-negative
    exp_duration = max(0, int(random.gauss(exp_mean, exp_std)))  # Ensure exp is non-negative

    # Step 2: Initialize variables
    time = list(range(total_time))  # Time points
    population = [P0] * total_time  # Population sizes (OD) at each time point

    # Step 3: Simulate logistic growth
    for t in range(1, total_time):
        if t < lag_duration:
            # Lag phase: slow growth
            population[t] = population[t - 1]
        elif t < lag_duration + exp_duration:
            # Exponential phase: rapid growth
            delta_population = r * population[t - 1] * (1 - population[t - 1] / K)
            population[t] = population[t - 1] + delta_population
        else:
            # Stationary phase: growth levels off
            population[t] = K

    return time, population

# Example usage
K = 1000  # Carrying capacity
P0 = 10   # Initial population size
r = 0.2   # Growth rate
total_time = 100  # Total simulation time
lag_mean = 10     # Mean lag phase duration
lag_std = 2       # Standard deviation of lag phase
exp_mean = 30     # Mean exponential phase duration
exp_std = 5       # Standard deviation of exponential phase

# Simulate growth
time, population = logistic_growth(K, P0, r, total_time, lag_mean, lag_std, exp_mean, exp_std)

# Print the growth curve (OD vs Time)
print("Time\tOD")
for t, od in zip(time, population):
    print(f"{t}\t{od:.2f}")

# ==============================================
# Section 3: Generate a DataFrame with 100 Different Growth Curves
# ==============================================

# Function to generate multiple logistic growth curves and store them in a DataFrame

def generate_growth_curves(num_curves, K, P0, r, total_time, lag_mean, lag_std, exp_mean, exp_std):
    """
    Generates multiple logistic growth curves and stores them in a DataFrame.

    Parameters:
        num_curves (int): Number of growth curves to generate.
        K (float): Carrying capacity.
        P0 (float): Initial population size.
        r (float): Growth rate.
        total_time (int): Total simulation time.
        lag_mean (float): Mean duration of the lag phase.
        lag_std (float): Standard deviation of the lag phase.
        exp_mean (float): Mean duration of the exponential phase.
        exp_std (float): Standard deviation of the exponential phase.

    Returns:
        df (pd.DataFrame): DataFrame containing time and population data for all curves.
    """
    # Initialize an empty DataFrame
    df = pd.DataFrame()

    # Generate growth curves
    for curve_id in range(1, num_curves + 1):
        time, population = logistic_growth(K, P0, r, total_time, lag_mean, lag_std, exp_mean, exp_std)
        
        # Create a temporary DataFrame for the current curve
        temp_df = pd.DataFrame({
            "curve_id": curve_id,
            "time": time,
            "population": population
        })
        
        # Append the temporary DataFrame to the main DataFrame
        df = pd.concat([df, temp_df], ignore_index=True)

    return df

# Step 5: Example usage
K = 1000  # Carrying capacity
P0 = 10   # Initial population size
r = 0.2   # Growth rate
total_time = 100  # Total simulation time
lag_mean = 10     # Mean lag phase duration
lag_std = 2       # Standard deviation of lag phase
exp_mean = 30     # Mean exponential phase duration
exp_std = 5       # Standard deviation of exponential phase

# Generate 100 growth curves
growth_df = generate_growth_curves(100, K, P0, r, total_time, lag_mean, lag_std, exp_mean, exp_std)

# Display the first few rows of the DataFrame
print(growth_df.head())

# Function to determine the time to reach 80% of the carrying capacity

def time_to_80_percent(time, population, K):
    """
    Determines the time to reach 80% of the carrying capacity.

    Parameters:
        time (list): List of time points.
        population (list): List of population sizes (OD) at each time point.
        K (float): Carrying capacity.

    Returns:
        t_80 (int or None): Time to reach 80% of carrying capacity, or None if not reached.
    """
    target_population = 0.8 * K
    for t, pop in zip(time, population):
        if pop >= target_population:
            return t
    return None

# Generate 100 growth curves
growth_df = generate_growth_curves(100, K, P0, r, total_time, lag_mean, lag_std, exp_mean, exp_std)

# Convert the DataFrame to a list of dictionaries for easier processing
growth_curves = growth_df.groupby('curve_id').apply(lambda x: x.to_dict(orient='list')).tolist()

# Determine time to reach 80% of carrying capacity for each curve
times_to_80_percent = []
for curve in growth_curves:
    curve_id = curve["curve_id"][0]
    time = curve["time"]
    population = curve["population"]
    K = 1000  # Carrying capacity

    # Calculate time to reach 80% of K
    t_80 = time_to_80_percent(time, population, K)
    times_to_80_percent.append(t_80)

    # Print the result
    if t_80 is not None:
        print(f"Curve {curve_id}: Time to reach 80% of carrying capacity = {t_80}")
    else:
        print(f"Curve {curve_id}: Did not reach 80% of carrying capacity.")

# Plot the growth curves
plt.figure(figsize=(10, 6))
for curve in growth_curves:
    time = curve["time"]
    population = curve["population"]
    plt.plot(time, population, alpha=0.5)

plt.xlabel('Time')
plt.ylabel('Population (OD)')
plt.title('Logistic Growth Curves')
plt.show()

# Plot histogram of times to reach 80% of carrying capacity for better visualization
plt.figure(figsize=(10, 6))
plt.hist([t for t in times_to_80_percent if t is not None], bins=20, edgecolor='black')
plt.xlabel('Time to reach 80% of carrying capacity')
plt.ylabel('Frequency')
plt.title('Histogram of Times to Reach 80% of Carrying Capacity')
plt.show()

# ==============================================
# Section 4: Calculating The Hamming Distance
# ==============================================

# Function to calculate the Hamming distance between two strings


def hamming_distance(str1, str2):
    # Pad the shorter string with spaces to make both strings of equal length
    max_len = max(len(str1), len(str2))
    str1 = str1.ljust(max_len)
    str2 = str2.ljust(max_len)
    
    # Calculate the Hamming distance
    distance = sum(1 for x, y in zip(str1, str2) if x != y)
    return distance

# Example usage
slack_username = "Adams"
linkedin_handle = "J.Adams Nyaata"
distance = hamming_distance(slack_username, linkedin_handle)
print(f"Hamming Distance: {distance}")

# ==============================================
# Github Link: 
# ==============================================

j_adams_nyaata_github_link = "https://github.com/John-Adams-N/HackBio-Biocoding"

print(f"Github Link: {j_adams_nyaata_github_link}")
