# scripts/find_best_hit.py
#
# This script is executed by Snakemake.
# It reads a BLAST tabular output file (format 6) and, for each unique query,
# it identifies the single best hit based on the highest bitscore.
# This version uses the pandas library for efficient and readable data manipulation.

import pandas as pd
import sys

# --- Snakemake I/O ---
# The snakemake object is injected by the script directive in the Snakefile.
try:
    input_file = snakemake.input.blast_result
    output_file = snakemake.output[0]
except NameError:
    # This block allows the script to be run standalone for testing
    # e.g., python scripts/find_best_hit.py input.blast.tsv output.best_hit.tsv
    if len(sys.argv) != 3:
        print("Usage: python find_best_hit.py <input_blast_tsv> <output_best_hit_tsv>")
        sys.exit(1)
    input_file = sys.argv[1]
    output_file = sys.argv[2]


# --- Main Logic ---

# Define the column names for BLAST output format 6
blast_columns = [
    'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
    'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
]

try:
    # Read the BLAST output file into a pandas DataFrame.
    # If the file is empty (no hits), this will raise an error which we catch below.
    df = pd.read_csv(
        input_file,
        sep='\t',
        header=None,
        names=blast_columns
    )

    if not df.empty:
        # Sort the entire DataFrame by bitscore in descending order.
        # This brings the best hits for each query to the top of their respective groups.
        df.sort_values('bitscore', ascending=False, inplace=True)

        # Drop duplicate query IDs, keeping only the first occurrence,
        # which is now guaranteed to be the best hit due to the sort.
        best_hits_df = df.drop_duplicates('qseqid', keep='first')

        # Write the resulting DataFrame to the output file.
        best_hits_df.to_csv(output_file, sep='\t', index=False)

    else:
        # If the input file was empty, create an empty output file with just the header.
        pd.DataFrame(columns=blast_columns).to_csv(output_file, sep='\t', index=False)

except pd.errors.EmptyDataError:
    # This handles the case where blastn produced no hits and the input file is empty.
    # We create an empty output file with the correct headers.
    print(f"Input file {input_file} is empty. Creating an empty output file.")
    pd.DataFrame(columns=blast_columns).to_csv(output_file, sep='\t', index=False)


