import argparse
import pandas as pd
import numpy as np
import os
import re
import logging

def process_max_similarities(input_file, nbins, is_ref):
    """
    Processes a similarity matrix to compute the maximum similarity per row,
    and bins the column max similarities into histograms. It also saves the results.

    Parameters:
        input_file (str): Path to the input similarity matrix file (CSV, compressed).
        nbins (int): Number of bins for the histogram.
        is_ref (bool): Whether the data is a reference dataset (diagonal 1s).
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Reading input file: {input_file}")
    sim_matrix = pd.read_csv(input_file, compression="gzip", header=0)
    logger.info("File loaded successfully.")

    # If reference dataset, overwrite diagonal elements
    if is_ref:
        logger.info("Detected reference dataset. Overwriting diagonal elements.")

        # Extract the sub-comparison's starting index from the filename
        match = re.search(r"_(\d{4})_", os.path.basename(input_file))
        if not match:
            raise ValueError("Filename does not contain a valid sub-comparison number.")
        subcomp_index = int(match.group(1))

        # Compute the starting row index for this sub-comparison
        rows_per_subcomp = 5944
        starting_row_index = subcomp_index * rows_per_subcomp

        # Overwrite the diagonal
        for col_index in range(sim_matrix.shape[1]):
            row_index = starting_row_index + col_index
            if row_index < sim_matrix.shape[0]:  # Ensure valid index
                sim_matrix.iloc[row_index, col_index] = -np.inf

    # Compute the maximum similarity per row
    logger.info("Extracting maximum similarity values per row.")
    max_similarities_row = sim_matrix.max(axis=1)

    # Compute the maximum similarity per column 
    logger.info("Extracting maximum similarity values per column.")
    max_similarities_column = sim_matrix.max(axis=0)

    # Bin the maximum similarities (only for column-wise)
    logger.info("Binning the maximum similarity values (column-wise).")
    bin_edges = np.linspace(0, 1, nbins + 1)
    binned_values_column, _ = np.histogram(max_similarities_column, bins=bin_edges)

    # Create output filenames
    input_dir = os.path.dirname(input_file)
    base_name = os.path.basename(input_file).replace(".csv.gz", "")
    binned_output_file = os.path.join(input_dir, f"{base_name}_max_bins.csv.gz")
    row_max_output_file = os.path.join(input_dir, f"{base_name}_each_ref_max.csv.gz")

    # Save the raw maximum similarity per row (without binning)
    logger.info(f"Saving maximum similarity per row to: {row_max_output_file}")
    max_similarities_row.to_csv(row_max_output_file, compression="gzip", index=False, header=False)

    # Save the binned data (column-wise) to the output file
    logger.info(f"Saving binned data (column-wise) to output file: {binned_output_file}")
    pd.Series(binned_values_column).to_csv(binned_output_file, compression="gzip", index=False, header=False)

    logger.info("Processing complete.")

if __name__ == "__main__":
    # Argument parser setup
    parser = argparse.ArgumentParser(description="Process max similarities and bin them.")
    parser.add_argument("-i", "--input", required=True, help="Input similarity matrix file (gzip).")
    parser.add_argument("-nb", "--nbins", type=int, required=True, help="Number of bins for histogram.")
    parser.add_argument("-r", "--ref", type=bool, default=False, help="Flag indicating if dataset is a reference set.")
    args = parser.parse_args()

    # Logging setup
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    # Run the processing
    process_max_similarities(args.input, args.nbins, args.ref)

