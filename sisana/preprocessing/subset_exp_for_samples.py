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
    
    parser = argparse.ArgumentParser(description="Example command: python subset_exp_for_samples.py -e exp_filtered.csv -o ./output/exp_subset.tsv")
    ArgGroup = parser.add_argument_group('Required arguments')  
    ArgGroup.add_argument("-e", "--expfile", type=str, help="Path to csv expression file", required=True) 
    ArgGroup.add_argument("-n", "--nameorder", type=str, help="Path to a .txt file containing the names of samples in order they appear in the expression file", required=True) 
    ArgGroup.add_argument("-s", "--sampnames", type=str, help="Path to a .txt file containing a list of sample names to extract, one name per line", required=True) 
    ArgGroup.add_argument("-o", "--outdir", type=str, help="Name of the directory to save output to", required=True) 
    
    args = parser.parse_args()
    
    # Create output directory if one does not already exist
    os.makedirs(args.outdir, exist_ok=True)
    
    exp = pd.read_csv(args.expfile, sep='\t', header=None, index_col=0)
    
    sampfile = open(args.nameorder, "r")
    samp_list = [i.strip() for i in sampfile.readlines()]
    print(samp_list)
    
    subsetfile = open(args.sampnames, "r")
    subset_list = [i.strip() for i in subsetfile.readlines()]
    print(subset_list)
    
    exp.columns = samp_list

    subset_exp = exp.loc[:, subset_list]
    
    print(subset_exp)
    
    base_file_name = Path(args.expfile).stem

    subset_exp.to_csv(f"{args.outdir}/{base_file_name}_subset.tsv", sep='\t', index=True, header=False)
    
    print(f"Expression file saved to: {args.outdir}/{base_file_name}_subset.tsv")
    
    with open(f"{args.outdir}/{base_file_name}_subset_sample_order.txt", 'w') as f:
        for samp in subset_exp.columns:
            f.write(f"{samp}\n")
    
    print(f"Order of samples in expression file were saved to: {args.outdir}/{base_file_name}_subset_sample_order.txt")
