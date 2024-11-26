import argparse
import pandas as pd
import numpy as np
import os
import logging

def process_max_similarities(input_dir, nbins):
    """
    Processes a similarity matrix to compute the global maximum similarity per row across subcomparisons
    and bins the values into a histogram.

    Parameters:
        input_dir (str): Path to the directory containing similarity matrix files (gzip).
        nbins (int): Number of bins for the histogram.
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Processing files in directory: {input_dir}")

    # Initialize an array to store the global maximum similarity values for each row
    row_max_similarities = None

    # Loop over all files in the directory and process each subcomparison file
    for filename in os.listdir(input_dir):
        if filename.endswith("_each_ref_max.csv.gz"):  # Ensure we're processing the correct files
            #logger.info(f"Processing file: {filename}")

            # Read the subcomparison file
            file_path = os.path.join(input_dir, filename)
            subcomparison_data = pd.read_csv(file_path, compression="gzip", header=None)

            # If it's the first file, initialize the row_max_similarities array
            if row_max_similarities is None:
                row_max_similarities = np.zeros(subcomparison_data.shape[0])

            # Update the global row max similarities (taking the max across all subcomparisons)
            row_max_similarities = np.maximum(row_max_similarities, subcomparison_data.max(axis=1).values)

    # After processing all files, we now have the global maximum similarity for each row
 #   logger.info(f"Global row max similarities: {row_max_similarities}")

    # Binning the global row max similarities into specified bins
    logger.info(f"Binning the global row max similarities into {nbins} bins.")
    bin_edges = np.linspace(0, 1, nbins + 1)
    
    # Bin the global row max similarities
    binned_row_max_values, _ = np.histogram(row_max_similarities, bins=bin_edges)

    # Output the binned values (for the global row max similarities)
    binned_output_file = os.path.join(input_dir, "global_row_max_bins.csv.gz")
    pd.Series(binned_row_max_values).to_csv(binned_output_file, index=False)

    logger.info(f"Saved binned global row max similarities to: {binned_output_file}")
    logger.info("Processing complete.")

if __name__ == "__main__":
    # Argument parser setup
    parser = argparse.ArgumentParser(description="Process max similarities and bin them.")
    parser.add_argument("-d", "--directory", required=True, help="Directory containing similarity matrix files.")
    parser.add_argument("-nb", "--nbins", type=int, required=True, help="Number of bins for histogram.")
    args = parser.parse_args()

    # Logging setup
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    # Run the processing
    process_max_similarities(args.directory, args.nbins)
