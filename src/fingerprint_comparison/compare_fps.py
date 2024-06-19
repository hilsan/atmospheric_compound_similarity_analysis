#Need to read datasets in chunks to create histogram. 
import sys, os
#Load modules
import logging
import time
from joblib import Parallel, delayed
from rdkit import Chem
from rdkit.Chem import (AllChem, MACCSkeys, RDKFingerprint, DataStructs)
from rdkit.DataStructs import BulkTanimotoSimilarity
from rdkit import Chem
from argparse import ArgumentParser, RawTextHelpFormatter
import numpy as np
import pandas as pd
import joblib as jl

#------- Define command-line arguments
parser = ArgumentParser(description='''
Compare the datasets  ''',formatter_class=RawTextHelpFormatter)
###_ . Arguments
parser.add_argument('-k1','--key1',type=str,
                    help='name of dataset 1 data dump with fingerprint')
parser.add_argument('-k2','--key2',type=str,
                    help='name of dataset 2 datadump with fingerprint')
parser.add_argument('-n','--n_jobs',type=int,default=2,
                    help='Number of processes')
parser.add_argument('-fp','--fptype',type=str,default='MACCS',
                    help='Type of fingerprint')  
#----------------------------------------------------------
def similarity(molfp, df2, fptype): 
    return df2[fptype].apply(lambda x: DataStructs.FingerprintSimilarity(x, molfp, metric=DataStructs.TanimotoSimilarity))

def get_df(key): 
    name = '../'+key+'_fp.csv'
    if os.path.isfile(name):
        df = pd.read_csv(name,compression={'method':'gzip', 'compresslevel': 1,'mtime': 1})
        return df
    else:
        print('cant find dataset 1 fingerprint')
    
def get_batch_df(key): 
    name = key+'.csv'
    d = pd.read_csv(name, compression={'method':'gzip', 'compresslevel': 1,'mtime': 1})
    return d
    
def write_comparison(sim_df,fptype, key1, outfile):
    if os.path.isdir('comp_'+key1+'_'+fptype):
        pass
    else:
        os.mkdir('comp_'+key1+'_'+fptype)       
    if 'index' in sim_df.columns: 
        sim_df.drop(columns='index').astype('float32').to_csv(''.join(['comp_',key1,'_',fptype,'/',outfile,'.csv']),compression={'method': 'gzip', 'compresslevel': 1,'mtime': 1}, index=False)
    else: 
        sim_df.astype('float32').to_csv(''.join(['comp_',key1,'_',fptype,'/',outfile,'.csv']),compression={'method': 'gzip', 'compresslevel': 1,'mtime': 1}, index=False)
        
def bin_batch_data(chunk, filename, fptype, key1, key2, nbins = 200):
    if os.path.isdir('bins_'+key1+'_'+fptype):
        pass
    else:
        os.mkdir('bins_'+key1+'_'+fptype)  
    bin_edges = np.linspace(0, 1, nbins + 1)
    if 'index' in chunk.columns:
        chunk=chunk.drop(columns='index')
    # compute bin counts over the 3rd column
    subtotal, e = np.histogram(chunk.stack(), bins=bin_edges)
    pd.Series(subtotal).astype('int32').to_csv(''.join(['bins_',key1,'_',fptype,'/',filename,'_binned_total.csv']),compression={'method': 'gzip', 'compresslevel': 1,'mtime': 1}, index=False)

def tan_sim(i, d2, d1):
    return DataStructs.BulkTanimotoSimilarity(d2[i], d1)

def main():
    #------------Main-------------------------
    logger = logging.getLogger(__name__)
    logger.info('Parse data')
    args = parser.parse_args()
    logger.info('Load reference dataset fingerprint')
    df1 = get_df(args.key1)
    logger.info('Loaded:'+args.key1+'_fp.csv')
    logger.info('Load batch dataset fingerprint')
    df2 = get_batch_df(args.key2)
    logger.info('Loaded:'+args.key2+'.csv')
    df2[args.fptype]=df2[args.fptype].apply(DataStructs.CreateFromBitString)
    df1[args.fptype]=df1[args.fptype].apply(DataStructs.CreateFromBitString)
    #-------------------------------------
    logger.info('Compute comparison using '+args.fptype+' fingerprint')
    t = time.time()
    outfile = args.key2+'_'+args.key1+'_'+args.fptype
    results = Parallel(n_jobs=args.n_jobs)(delayed(BulkTanimotoSimilarity)(x, df2[args.fptype]) for x in df1[args.fptype].values)
    sim_df = pd.DataFrame(results)
    logger.info('Write similarities to file')
    write_comparison(sim_df, args.fptype, args.key1, outfile)
    logger.info('Bin similarities into a subtotal')
    bin_batch_data(sim_df, outfile, args.fptype, args.key1, args.key2, nbins = 200)
    print(time.time() - t)
        
if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    main()
        