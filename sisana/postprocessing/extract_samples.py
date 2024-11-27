import networkx as nx
import argparse
import os
import pandas as pd
import pickle
import numpy as np
from pathlib import Path
   
__author__ = 'Nolan Newman'
__contact__ = 'nolankn@uio.no'
    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Example command: python extract_samples.py -d <file.pickle> -f pickle -n samps.txt -o ./all_nws.pickle -s 1064_samps")
    ArgGroup = parser.add_argument_group('Required arguments')  
    ArgGroup.add_argument("-d", "--datafile", type=str, help="Path to pickle file created by lioness_to_pickle_df.py script", required=True) 
    ArgGroup.add_argument("-f", "--filetype", choices = ["pickle", "csv"], type=str, help="File type of the input data file", required=True)
    ArgGroup.add_argument("-n", "--sampnames", type=str, help="Path to a .txt file containing a list of sample names to extract, one name per line, no header", required=True) 
    ArgGroup.add_argument("-o", "--outdir", type=str, help="Path to directory where file will be output, which contains the input df, subset for just the specified samples", required=True) 
    ArgGroup.add_argument("-s", "--suffix", type=str, help="Suffix to add to the output file, used to accurately name the subsetting that was performed", required=True) 

    args = parser.parse_args()
    
    if args.filetype == "pickle":
        lion = pd.read_pickle(args.datafile, index_col = 0)   
    elif args.filetype == "csv":
        lion = pd.read_csv(args.datafile, index_col = 0)  
    
    sampfile = open(args.sampnames, "r")
    savedir = args.outdir

    # strip new line chars from samp names
    samp_list = [i.strip() for i in sampfile.readlines()]
    
    print(f"First few sample names: {samp_list[0:4]}") 
    print("Structure of input data file:")
    print(lion)
    print("Now extracting samples from input file, please wait...")
    
    subsetlion = lion.loc[:, samp_list]
    print(subsetlion)

    base_file_name = Path(args.datafile).stem

    if args.filetype == "pickle":
        filename = f"{base_file_name}_subset_{args.suffix}.pickle"
        save_file_path = os.path.join(savedir, filename)
        subsetlion.to_pickle(save_file_path)   
    elif args.filetype == "csv":
        filename = f"{base_file_name}_subset_{args.suffix}.csv"
        save_file_path = os.path.join(savedir, filename)
        subsetlion.to_csv(save_file_path)   

    print(f"File saved: {save_file_path}")
 
