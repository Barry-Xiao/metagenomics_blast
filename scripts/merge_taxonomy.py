import pandas as pd
import sys

# Log to file
log_file = snakemake.log[0]
sys.stderr = sys.stdout = open(log_file, "w")

try:
    # --- Load Data ---
    # Load the best hits from the find_best_hit rule
    hits_df = pd.read_csv(snakemake.input.hits, sep='\t')

    # Load the taxonomy data fetched from NCBI
    # The xtract output from entrez-direct has no header, so we name the columns
    tax_df = pd.read_csv(
        snakemake.input.tax,
        sep='\t',
        header=None,
        names=['sseqid', 'organism', 'taxid']
    )

    # --- Merge Data ---
    # Perform a left merge. This keeps all original blast hits, even if
    # a taxonomy lookup failed for some reason. Failed lookups will have NaN
    # in the 'organism' and 'taxid' columns.
    merged_df = pd.merge(
        hits_df,
        tax_df,
        on='sseqid',
        how='left'
    )

    # --- Save Output ---
    # Save the final annotated table to the output file
    merged_df.to_csv(snakemake.output[0], sep='\t', index=False)

except Exception as e:
    print(f"Error merging taxonomy data: {e}")
    # Create an empty file to prevent Snakemake from failing on subsequent rules
    # if there was an error (e.g., no BLAST hits to begin with).
    open(snakemake.output[0], 'a').close()
    sys.exit(1)

