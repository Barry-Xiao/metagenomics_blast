# scripts/merge_outputs.py
#
# This script is executed by Snakemake as the final step in the pipeline.
# It reads all the individual annotated TSV files, adds a 'sample' column
# to each one, and concatenates them into a single, final output file.

import pandas as pd
import sys
import os

# --- Snakemake I/O ---
try:
    input_files = snakemake.input
    output_file = snakemake.output[0]
except NameError:
    # This block allows the script to be run standalone for testing
    print("Usage: python merge_outputs.py <output_file> <input_file1> <input_file2> ...")
    sys.exit(1)


# --- Main Logic ---

# A list to hold the DataFrame for each sample
all_dfs = []

# Iterate over all the input files provided by Snakemake
for f in input_files:
    # Check if the file is not empty to avoid errors
    if os.path.getsize(f) > 0:
        try:
            # Read the annotated TSV into a DataFrame
            df = pd.read_csv(f, sep='\t')

            # Extract the sample name from the file path
            # e.g., 'results/annotated_hits/sample1.annotated.tsv' -> 'sample1'
            sample_name = os.path.basename(f).split('.')[0]

            # Add a new column with the sample name
            df['sample'] = sample_name

            # Add the processed DataFrame to our list
            all_dfs.append(df)
        except pd.errors.EmptyDataError:
            print(f"Warning: File {f} is empty or just a header. Skipping.")
            continue

# Concatenate all the DataFrames in the list into a single one
if all_dfs:
    merged_df = pd.concat(all_dfs, ignore_index=True)

    # Reorder columns to have 'sample' first for better readability
    cols = ['sample'] + [col for col in merged_df.columns if col != 'sample']
    merged_df = merged_df[cols]

    # Write the final merged DataFrame to the output file
    merged_df.to_csv(output_file, sep='\t', index=False)
else:
    # If all input files were empty, create an empty output file
    print("All input files were empty. Creating an empty merged output file.")
    with open(output_file, 'w') as out:
        pass # Creates an empty file
