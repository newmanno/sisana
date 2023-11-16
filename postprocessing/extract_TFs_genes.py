import argparse
import sys
import pickle
import pandas as pd
import numpy as np
from pathlib import Path
import os

if __name__ == '__main__':
    """
    Description:
        This code filters the lioness output file for only desired TFs 
    """
    
    parser = argparse.ArgumentParser(description="Example command: python extract_TFs_genes.py -q lioness_output.pickle -o ./output/lioness_df.pickle")
    requiredArgGroup = parser.add_argument_group('Required arguments')  
    requiredArgGroup.add_argument("-p", "--picklefile", type=str, help="Path to lioness output in pickle format", required=True) 
    requiredArgGroup.add_argument("-t", "--tfs", type=str, nargs = '+', help="Name of TFs to subset for", required=True)     
    requiredArgGroup.add_argument("-o", "--outdir", type=str, help="Path to output directory", required=True)    
    
    args = parser.parse_args()
    
    lion = pd.read_pickle(args.picklefile)    
    
    # Separate the row names to be two different columns, TF and Target
    lion = lion.rename_axis("TF").reset_index()
    lion[['TF','Target']] = lion['TF'].str.split('<==>',expand=True)   

    col = lion.pop("Target")
    lion.insert(1, col.name, col)
    lion_subset = lion.loc[lion['TF'].isin(args.tfs)]
    
    base_file_name = Path(args.picklefile).stem
    out_filename =  f"{base_file_name}_filtered_for_TFs.csv"
    save_file_path = os.path.join(args.outdir, out_filename)
    lion_subset.to_csv(save_file_path, index = False)

    print(f"\nFile saved: {save_file_path}\n")  
