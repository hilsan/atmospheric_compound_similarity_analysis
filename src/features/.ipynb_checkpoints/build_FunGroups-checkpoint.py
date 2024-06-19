"""
-Build descriptors in the form of: 
    Molecular weight
    Atomic composition 
    Atomic ratios
    Heavy atoms
-Maccs computation 
-Topological fingerprint computation
-Fingerprint comparison 
-Fun. group analysis
"""
#Load modules
import click
import logging
import joblib as jl 
import pandas as pd
import sys, os
import numpy as np
from aprl_ssp import substructure_search
from load_datasets import load_datasets
#----------------------------------------------------    
def compute_functional_groups(dataframe, key, groupfile, input_filepath, output_path):
    """Compute functional groups with simpol substructure_search code"""
    inputfile = input_filepath+key+"_to_simpol.csv"
    outputfile = output_path+key+'_SIMPOLgroups.csv'
    dataframe['compound'] = dataframe.index
    dataframe.to_csv(inputfile, index_label='compound') ##Here i can write less
    groups = substructure_search.get_smarts_groups(groupfile)
    inp = substructure_search.get_smiles(inputfile)
    export = None 
    search = substructure_search.search_groups(groups, export) 
    output = substructure_search.count_groups(inp, search)
    substructure_search.write_output(output, outputfile)

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
    for key in datasets.keys():        
        logger.info('Starting analysis on '+key)
        dataframe = jl.load(input_filepath+datasets[key])        
        logger.info('Building functional groups...')
        compute_functional_groups(dataframe, key,'/home/sandsth2/aprl_ssp/SMARTSpatterns/SIMPOLgroups_noring_nonnitrophenol.csv', input_filepath, output_filepath)

if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    main()
