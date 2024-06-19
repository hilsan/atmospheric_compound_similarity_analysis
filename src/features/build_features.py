"""
-Build descriptors in the form of: 
    Molecular weight
    Atomic composition 
    Atomic ratios
    Heavy atoms
-Maccs fingerprint
-Topological fingerprint 
-Fun. group analysis
"""
#Load modules
import click
import logging
import pandas as pd
import joblib as jl 
import rdkit; 
import sys, os
import numpy as np
from load_datasets import load_datasets
from build_FunGroups import compute_functional_groups
from build_fingerprints import compute_maccs_fingerprint, compute_topological_fingerprint, split_data_set, write_gecko
from build_descriptors import build_descriptors
#---------------------------------------------------- 
@click.command()
@click.argument('input_filepath', type=click.Path())
@click.argument('output_filepath', type=click.Path())
@click.argument('reports_filepath', type=click.Path())  

def main(input_filepath, output_filepath, reports_filepath):
    """Runs data processing scripts to turn processed data from (../processed) 
    into data features ready to be analyzed (saved in ../reports)."""
    logger = logging.getLogger(__name__)
    datasets = load_datasets()
    maccs_analysis = pd.DataFrame(columns = ['name', 'feature_counts', 'nz_features', 'n_shared_features'])
    percentages = np.linspace(0.0,0.95,11)
    batch_size=6000
    for key in datasets.keys():
        
        logger.info('Starting analysis on '+key)
        dataframe = jl.load(input_filepath+datasets[key])

        logger.info('Compute MACCS fingerprint for dataset')
        maccs_analysis.loc[len(maccs_analysis)] = compute_maccs_fingerprint(dataframe, key, percentages)

        logger.info('Compute topological fingerprint for dataset')
        compute_topological_fingerprint(dataframe)

        logger.info('Splitting dataset in prep. of fingerprint comparison')
        split_data_set(dataframe, output_filepath, key, batch_size)
        if key == "Gecko" or key == "Wang":
                write_gecko(dataframe,key, output_filepath)
        

        logger.info('Building functional groups...')
        compute_functional_groups(dataframe, key,'/home/sandsth2/aprl_ssp/SMARTSpatterns/SIMPOLgroups_noring_nonnitrophenol.csv', input_filepath, output_filepath)

        logger.info('Building molecular properties...') #Needs to come last at the moment (2 level columns)
        build_descriptors(dataframe, key, output_filepath, reports_filepath)

        maccs_analysis.to_csv(reports_filepath+'maccs_compiled.csv')
    
if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    main()
