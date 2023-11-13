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
    
    parser = argparse.ArgumentParser(description="Example command: python extract_samples.py -e <file.pickle> -o ./all_nws.pickle")
    ArgGroup = parser.add_argument_group('Required arguments')  
    ArgGroup.add_argument("-p", "--picklefile", type=str, help="Path to pickle file created by lioness_to_pickle_df.py script", required=True) 
    ArgGroup.add_argument("-s", "--sampnames", type=str, help="Path to a .txt file containing a list of sample names to extract, one name per line", required=True) 
    ArgGroup.add_argument("-o", "--outdir", type=str, help="Path to pickle file to be output, which contains the input df, subset for just the specified samples", required=True) 
    
    args = parser.parse_args()
    
    lion = pd.read_pickle(args.picklefile)   
    sampfile = open(args.sampnames, "r")
    savedir = args.outdir

    
    # strip new line chars from samp names
    samp_list = [i.strip() for i in sampfile.readlines()]
    
    print(samp_list) 
    
    print(lion)
    
    subsetlion = lion.loc[:, samp_list]
    
    base_file_name = Path(args.picklefile).stem
    print(base_file_name)

    print(samp_list[0])

    print(samp_list[len(samp_list)-1])


    filename =  f"{base_file_name}_{samp_list[0]}_to_{samp_list[len(samp_list)-1]}.pickle"
    save_file_path = os.path.join(savedir, filename)
    
    subsetlion.to_pickle(save_file_path)   
    
    print(f"File saved: {save_file_path}")
 