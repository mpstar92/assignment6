import pandas as pd
import glob
import os

# Get all RSEM genes.results files
input_files = glob.glob("results/rsem/*.genes.results")

# Empty list to hold all dataframes
dfs = []

# Process each file
for file in input_files:
    sample = os.path.basename(file).replace(".genes.results", "")
    df = pd.read_csv(file, sep="\t", usecols=["gene_id", "TPM"])
    df = df.rename(columns={"TPM": sample})
    dfs.append(df.set_index("gene_id"))

# Combine all into one matrix
merged = pd.concat(dfs, axis=1)

# Save output
merged.to_csv("results/summary/rsem_tpm_matrix.tsv", sep="\t")