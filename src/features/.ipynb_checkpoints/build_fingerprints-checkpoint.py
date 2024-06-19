"""
-Build descriptors in the form of: 
-Maccs computation 
-Topological fingerprint computation
-Fingerprint comparison 
-Fun. group analysis
"""
#Load modules
import click
import logging
import pandas as pd
import joblib as jl 
import sys, os
import numpy as np
from MolAnalysis.fingerprints import MACCSAnalysis, Fingerprints
from load_datasets import load_datasets
 #-------------------------------------------------------   
def split_data_set(df, outpath, key, batch_size):
    n_parts = np.ceil(len(df)/batch_size)
    df_batch = np.array_split(df, n_parts)
    counter = 0
    outpath_key = outpath+key
    write_columns = ['SMILES']
    if 'MACCS' in df.columns: 
        write_columns.append('MACCS')
    if 'topological'  in df.columns: 
        write_columns.append('topological')
    if os.path.isdir(outpath_key):
        pass
    else:
        os.mkdir(outpath_key)   
    for i_batch in df_batch:
        
        i_batch[write_columns].to_csv(outpath_key+'/'+key+'_'+str('{num:04d}'.format(num=counter))+'.csv',compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1}) #Here i can change types
        counter +=1

def write_gecko(df,key,outpath_key):
    write_columns = ['SMILES']
    if 'MACCS' in df.columns: 
        write_columns.append('MACCS')
    if 'topological'  in df.columns: 
        write_columns.append('topological')
    df.to_csv(outpath_key+'/'+key+'_fp.csv',compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1})

#----------------------------------------------------    
def compute_maccs_fingerprint(dataframe, key, percentages = np.linspace(0.0,0.95,11)):
    """Compute maccs fingerprints for the datasets and save in reports"""
    maccs_analysis = pd.DataFrame(columns = ['name', 'feature_counts', 'nz_features', 'n_shared_features'])
    analysis = MACCSAnalysis(dataframe)
    feature_counts, nz_features = analysis.get_counts()
    n_shared_features = analysis.shared_features(nz_features, percentages)
    dataframe['MACCS'] = dataframe['MACCS'].apply(lambda x: x.ToBitString())
    return [key, feature_counts.astype('int32'), nz_features.astype('float32'), np.array(n_shared_features).astype('int32')]

def compute_topological_fingerprint(dataframe):
    """Compute topological fingerprint with Lumiari max path, bits per hash and fp_size"""
    n_bits_per_hash = 2 
    max_path = 7 
    fp_size = 2048
    dataframe['topological'] = dataframe['SMILES'].apply(lambda x:
            Fingerprints(fp_type='topological', n_bits_per_hash = 2, max_path = 7, fp_size = 2048).get_fingerprint(x).ToBitString() if(np.all(pd.notnull(x)))
            else np.NaN)

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
    batch_size=1000
    for key in datasets.keys():
        
        logger.info('Starting analysis on '+key)
        dataframe = jl.load(input_filepath+datasets[key])
        
        logger.info('Compute MACCS fingerprint for dataset')
        maccs_analysis.loc[len(maccs_analysis)] = compute_maccs_fingerprint(dataframe, key, percentages)
        
        logger.info('Compute topological fingerprint for dataset')
        compute_topological_fingerprint(dataframe)
        
        logger.info('Splitting dataset in prep. of fingerprint comparison')
        split_data_set(dataframe, output_filepath, key, batch_size)
        
        if key == 'Gecko':
            write_gecko(dataframe,key, output_filepath)
        
    maccs_analysis.to_csv(reports_filepath+'maccs_compiled.csv')

if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    main()
