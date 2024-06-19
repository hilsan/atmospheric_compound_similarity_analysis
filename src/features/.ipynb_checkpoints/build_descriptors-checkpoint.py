"""
-Build descriptors in the form of: 
    Molecular weight
    Atomic composition 
    Atomic ratios
    Heavy atoms
"""
#Load modules
import click
import logging
import pandas as pd
import joblib as jl 
import sys, os
from MolAnalysis.descriptors import CustomDescriptors
from MolAnalysis.compileanalysis import *
from load_datasets import load_datasets
#----------------------------------------------------
def build_descriptors(dataframe, key, outpath,report_path):
    """Build properties of molecular size and composition"""
    dataframe = CustomDescriptors(dataframe)
    compile_analysis(dataframe, report_path+key)
    dataframe.to_csv(outpath+key+".descr.csv", compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1}) ###here i can change type
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
        logger.info('Building molecular properties...') #Needs to come last at the moment (2 level columns)
        build_descriptors(dataframe, key, output_filepath, reports_filepath)

if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    main()
