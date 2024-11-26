import argparse
import pandas as pd
import numpy as np
import os
import re
import logging

def process_max_similarity_for_reference(input_dir, nbins):
    """
    Processes the reference dataset, computes the maximum similarity per column, 
    and bins the results into histograms. It assumes the reference dataset 
    (diagonal elements are overwritten with -np.inf).

    Parameters:
        input_dir (str): Directory containing similarity matrix files.
        nbins (int): Number of bins for the histogram.

    Returns:
        None
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Reading files from directory: {input_dir}")

    # Get sorted list of files (sorted by the numeric part)
    files = sorted(
        [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith(".csv")],
        key=lambda x: int(x.split('_')[1])  # Sort files by the numeric part (XX)
    )

    # Prepare variables to track the global maximum similarities and column count
    global_max_column_similarities = []
    columns_processed = 0  # Track the number of columns processed globally

    # Loop through each file one at a time
    for file in files:
        logger.info(f"Processing file: {file}")
        sim_matrix = pd.read_csv(file, compression="gzip", header=0)

        # Overwrite diagonal elements with -np.inf for the reference dataset
        logger.info("Detected reference dataset. Overwriting diagonal elements.")
        
        # For each file, overwrite the diagonal based on columns_processed
        for col_index in range(sim_matrix.shape[1]):
            row_index = columns_processed + col_index  # Use columns_processed to compute the row index
            if row_index < sim_matrix.shape[0]:  # Ensure valid index
                sim_matrix.iloc[row_index, col_index] = -np.inf

        # Compute the maximum similarity for each column (max across rows)
        max_similarities_column = sim_matrix.max(axis=0)  # Max similarity per column

        # Append the max similarity values to the global list
        global_max_column_similarities.extend(max_similarities_column)

        # Update the number of columns processed
        columns_processed += sim_matrix.shape[1]
        logger.info(f"Columns processed so far: {columns_processed}")

    # Bin the column-wise maximum similarities
    logger.info("Binning the global maximum similarity values (column-wise).")
    bin_edges = np.linspace(0, 1, nbins + 1)
    binned_values_column, _ = np.histogram(global_max_column_similarities, bins=bin_edges)

    # Create output filenames
    output_dir = os.path.dirname(input_dir)
    max_similarity_output_file = os.path.join(output_dir, "max_column_similarity.csv.gz")
    binned_output_file = os.path.join(output_dir, "sum_total_row-col.csv.gz")

    # Save the global maximum similarity values for each column
    logger.info(f"Saving global max similarity per column to: {max_similarity_output_file}")
    pd.Series(global_max_column_similarities).to_csv(max_similarity_output_file, compression="gzip", index=False, header=False)

    # Save the binned column max similarity data to the output file
    logger.info(f"Saving binned column max similarities to: {binned_output_file}")
    pd.Series(binned_values_column).to_csv(binned_output_file, compression="gzip", index=False, header=False)

    logger.info("Processing complete.")


if __name__ == "__main__":
    # Argument parser setup
    parser = argparse.ArgumentParser(description="Process max similarities and bin them.")
    parser.add_argument("-i", "--input", required=True, help="Directory with similarity matrix files (gzip).")
    parser.add_argument("-nb", "--nbins", type=int, required=True, help="Number of bins for histogram.")
    args = parser.parse_args()

    # Logging setup
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    # Run the processing for reference dataset
    process_max_similarity_for_reference(args.input, args.nbins)

