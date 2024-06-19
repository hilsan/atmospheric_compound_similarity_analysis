import pandas as pd
import seaborn as sns
from sklearn.manifold import TSNE
from sklearn.decomposition import TruncatedSVD, PCA
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from joblib import Parallel, delayed
import mpl_toolkits.mplot3d  # noqa: F401
from sklearn import manifold, datasets
from rdkit.Chem import (DataStructs)
import os
import logging
import time
from argparse import ArgumentParser, RawTextHelpFormatter
#------- Define command-line arguments
parser = ArgumentParser(description='''
Compare the datasets  ''',formatter_class=RawTextHelpFormatter)
###_ . Arguments
parser.add_argument('-l','--locations',type=str, nargs='+', help='loc all datasets')
parser.add_argument('-t','--tags',type=str, nargs='+', help='tags all datasets')

parser.add_argument('-fp','--fptype',type=str,
                    help='fptype')
parser.add_argument('-i','--init',type=str, default='pca',
                    help='initilizer')
parser.add_argument('-o','--output',type=str,
                    help='output dir')
parser.add_argument('-n','--n_jobs',type=int, default=24,
                    help='number of processes')
parser.add_argument('-p','--perplexity',type=int, nargs='+', default=50,
                    help='tsne perplexity')


def pandas_scatter(data, title):
    fig, ax1 = plt.subplots(figsize=(7, 4), dpi=300, constrained_layout=True)
    sns.set_palette(sns.color_palette("colorblind"))
    plt.rc('font', size=10)
    sns.scatterplot(data=data, x="x", y="y", hue="tag", ax=ax1, s=1, alpha=0.85)
    ax1.xaxis.set_major_formatter(ticker.NullFormatter())
    ax1.yaxis.set_major_formatter(ticker.NullFormatter())
    ax1.set_xlabel('')
    ax1.set_ylabel('')
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 12}, fancybox=True, framealpha=0)
    plt.savefig(title+'.png', bbox_inches='tight', 
               transparent=True,
               pad_inches=0)

def fp_to_array(fp_):
    "Convert maccs string fingerprint to an array"
    try:
        bit_fp = np.frombuffer(fp_.encode(), 'u1') - ord('0')
    except:
        bit_fp = np.NaN
    return bit_fp


def get_batch_df(directory, fptype):
    fp_data = []
    for filename in os.listdir(directory):
        if '.csv' in filename:
            batch_fp = pd.read_csv(directory+'/'+filename, compression={'method': 'gzip', 'compresslevel': 1,'mtime': 1})
            fp_data.append(batch_fp[fptype].apply(fp_to_array))
    if fptype == 'MACCS':
        return np.stack(pd.concat(fp_data).values)[1:]
    else:
        return np.stack(pd.concat(fp_data).values)

def get_parallel_batch_df(directory, filename, fptype):
    if '.csv' in filename:
        batch_fp = pd.read_csv(directory+'/'+filename, compression={'method': 'gzip', 'compresslevel': 1,'mtime': 1})
        if fptype == 'MACCS':
            return batch_fp[fptype].apply(fp_to_array)[1:] 

        else:
            return batch_fp[fptype].apply(fp_to_array)

def pca_trans(data, pca_components = 50):
    pca = PCA(n_components=pca_components)
    pca_data = pca.fit_transform(np.concatenate(data))
    print('n components:',"%i" %pca_components)
    print('Variance explained:',"%.2f" %np.sum(pca.explained_variance_ratio_))
    return pca_data 

def get_color(d, tag):
    return [tag]*len(d)

def run_tsne(data, perplexity, init, n_jobs):
	n_components = 2  
	t_sne = manifold.TSNE(
		n_components=n_components,
		perplexity=perplexity,
		init=init,
		n_iter_without_progress=300,
		n_iter=5000,
		random_state=0,
		n_jobs=n_jobs)
	return t_sne.fit_transform(data)

def main():
    logger = logging.getLogger(__name__)
    args = parser.parse_args()
    locs = args.locations
    fptype = args.fptype
    init = args.init
    tags=args.tags
    d_list=[]
    c_list=[]
    t = time.time()
    logger.info('loading data..')
    for i,loc in enumerate(locs): 
        d = Parallel(n_jobs=args.n_jobs)(delayed(get_parallel_batch_df)(loc, filename, fptype) for filename in os.listdir(loc))
        d = np.stack(pd.concat(d).values)
        d_list.append(d)
        c_list.append(get_color(d,tags[i]))


    c_m = np.concatenate(c_list)
    logger.info('starting pca...')
    pca_m = pca_trans(d_list, pca_components = 50)

    for perplexity in args.perplexity:
        print('Perplexity: ', perplexity)
        logger.info('starting tsne fit..')
        f_t_t_sne = run_tsne(pca_m, perplexity, init, args.n_jobs)
        title = args.output+'/Perplexity='+str(perplexity)
        pandas_scatter(pd.concat([pd.DataFrame(f_t_t_sne), pd.DataFrame(c_m)], axis=1, ignore_index=True).rename({0:'x', 1:'y', 2:'tag'}, axis='columns'), title)
        print(time.time() - t)
        pd.concat([pd.DataFrame(f_t_t_sne), pd.DataFrame(c_m)], axis=1, ignore_index=True).rename({0:'x', 1:'y', 2:'tag'}, axis='columns').to_csv(title+'.csv', compression={'method': 'gzip', 'compresslevel': 1,'mtime': 1}, index=False)
#---------------------------------------------------- 
if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    main()
