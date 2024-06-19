#Need to read datasets in chunks to create histogram. 
import sys, os
import numpy as np
import pandas as pd
import joblib as jl
import logging
from argparse import ArgumentParser, RawTextHelpFormatter
#------- Define command-line arguments
parser = ArgumentParser(description='''
Compare the datasets  ''',formatter_class=RawTextHelpFormatter)
###_ . Arguments
parser.add_argument('-i','--input',type=str,
                    help='loc bins dataset 2')
parser.add_argument('-o','--output',type=str,
                    help='loc bins dataset 2')
parser.add_argument('-nb','--nbins',type=int, default = 200, help='loc bins dataset 2')
parser.add_argument('-r','--ref', type=bool, default=False)
def sum_bins(directory, nbins):
    bin_edges = np.linspace(0, 1, nbins + 1)
    total = np.zeros(nbins)
    count=0
    for filename in os.listdir(directory):
        if 'binned' in filename:
            subtotal = pd.read_csv(directory+'/'+filename, compression={'method': 'gzip', 'compresslevel': 1,'mtime': 1}, dtype=int)
            print('Analyzing subtotal #'+str(count))  
            count += 1  
            total += subtotal.stack().values
    return total

#---------------------------------------------------- 
def main():
    args = parser.parse_args()
    total = sum_bins(args.input, args.nbins)
    if args.ref == True: 
        total[-1] = total[-1] - np.sqrt(np.sum(total))
    pd.Series(total).to_csv(args.output+'sum_total.csv', index=False)

#---------------------------------------------------- 
if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    main()