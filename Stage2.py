# ==============================================
# Section 1: Importing Datasets & Ensuring Data Consistency
# ==============================================

# ==============================================
# 1.1: Loading Datasets & Setting Working Directory
# ==============================================

# importing the necessary libraries
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns

# Define the path to your working directory
working_dir = os.getenv("HACKBIO_WORKING_DIR", ".")

# Define full file paths
sift_csv = os.path.join(working_dir, "sift_dataset.csv")
foldx_csv = os.path.join(working_dir, "foldx_dataset.csv")

# Print file paths to debug
print(f"SIFT CSV Path: {sift_csv}")
print(f"FoldX CSV Path: {foldx_csv}")

# Load CSV files into Pandas DataFrames
try:
    sift_df = pd.read_csv(sift_csv)
    foldx_df = pd.read_csv(foldx_csv)
except FileNotFoundError as e:
    print(e)
    print("Loading datasets from URLs instead...")
    sift_url = "Stage 2 Multi BioProjects/Amino Acid Mutation Analysis/sift_dataset.csv"
    foldx_url = "Stage 2 Multi BioProjects/Amino Acid Mutation Analysis/foldx_dataset.csv"
    sift_df = pd.read_csv(sift_url)
    foldx_df = pd.read_csv(foldx_url)

# ==============================================
# 1.2: Ensuring Data Consistency
# ==============================================

# Ensure column names are correctly formatted (strip leading/trailing spaces if needed)
sift_df.columns = sift_df.columns.str.strip()
foldx_df.columns = foldx_df.columns.str.strip()

# ==============================================
# 2: Merging Datasets and Filtering Deleterious Mutations
# ==============================================

# Create 'specific_Protein_aa' column by concatenating 'Protein' and 'Amino_Acid'
sift_df["specific_Protein_aa"] = sift_df["Protein"] + "_" + sift_df["Amino_Acid"]
foldx_df["specific_Protein_aa"] = foldx_df["Protein"] + "_" + foldx_df["Amino_Acid"]

# Merge both datasets on 'specific_Protein_aa'
merged_df = pd.merge(sift_df, foldx_df, on="specific_Protein_aa", suffixes=('_sift', '_foldx'))

# Filter mutations that are deleterious in both function and structure
deleterious_mutations = merged_df[
    (merged_df["sift_Score"] < 0.05) & (merged_df["foldX_Score"] > 2)
]

# ==============================================
# 3: Analyzing Deleterious Mutations
# ==============================================

# Extract the first amino acid from the 'Amino_Acid' column
deleterious_mutations.loc[:, "First_AA"] = deleterious_mutations["Amino_Acid_sift"].str[0]

# Generate frequency table for amino acids
amino_freq = deleterious_mutations["First_AA"].value_counts()

# Display results
print("\nDeleterious Mutations:")
print(deleterious_mutations)

print("\nAmino Acid Frequency Table:")
print(amino_freq)

# Save the table as an image
fig, ax = plt.subplots(figsize=(8, 4))
ax.axis('tight')
ax.axis('off')
table_data = [deleterious_mutations.columns.tolist()] + deleterious_mutations.head(10).values.tolist()
table = ax.table(cellText=table_data, colLabels=None, loc='center', cellLoc='center')
table.auto_set_font_size(False)
table.set_fontsize(8)
table.auto_set_column_width([i for i in range(len(deleterious_mutations.columns))])
plt.savefig("deleterious_mutations_table.png", dpi=300)
plt.show()

# ==============================================
# 4: Visualizing Results
# ==============================================

# Amino Acid Legend
amino_acid_dict = {
    "A": "Alanine", "R": "Arginine", "N": "Asparagine", "D": "Aspartic Acid", "C": "Cysteine",
    "E": "Glutamic Acid", "Q": "Glutamine", "G": "Glycine", "H": "Histidine", "I": "Isoleucine",
    "L": "Leucine", "K": "Lysine", "M": "Methionine", "F": "Phenylalanine", "P": "Proline",
    "S": "Serine", "T": "Threonine", "W": "Tryptophan", "Y": "Tyrosine", "V": "Valine"
}

# Convert dictionary to sorted lists
amino_acids = list(amino_freq.index)
counts = list(amino_freq.values)

# Sort by frequency
amino_acids, counts = zip(*sorted(zip(amino_acids, counts), key=lambda x: x[1], reverse=True))
labels = [f"{aa} - {amino_acid_dict[aa]}" for aa in amino_acids]

# Bar Plot
plt.figure(figsize=(12, 6))
sns.barplot(x=list(amino_acids), y=list(counts), palette="viridis")
plt.xlabel("Amino Acid", fontsize=14)
plt.ylabel("Frequency", fontsize=14)
plt.title("Frequency of Amino Acids in Deleterious Mutations", fontsize=16)
plt.xticks(rotation=45)
plt.legend(labels, title="Amino Acid Key", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.show()

# Pie Chart
plt.figure(figsize=(10, 10))
plt.pie(counts, labels=amino_acids, autopct=None, colors=sns.color_palette("viridis", len(amino_acids)))
plt.title("Proportion of Amino Acids in Deleterious Mutations", fontsize=16)
plt.legend([f"{aa} - {amino_acid_dict[aa]} ({count/sum(counts)*100:.1f}%)" for aa, count in zip(amino_acids, counts)], 
           title="Amino Acid Key", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.show()
