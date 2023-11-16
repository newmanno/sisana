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
    
    base_file_name = Path(args.picklefile).stem
    out_filename =  f"{base_file_name}_filtered_for_TFs.csv"
    save_file_path = os.path.join(args.outdir, out_filename)   
    
    # Separate the row names to be two different columns, TF and Target
    lion = lion.rename_axis("TF").reset_index()
    lion[['TF','Target']] = lion['TF'].str.split('<==>',expand=True)   
    lion['TF'] = lion['TF'].str.replace('TF_', '')

    # Check if the input TF(s) is/are in the data frame
    if set(args.tfs).issubset(np.unique(lion['TF'])):
        print("Specified TFs are in data frame, now filtering...")
    else:
        # Output a list of available TFs to filter file for
        out_filename_TFs =  f"{base_file_name}_available_TFs.csv"
        save_file_path_TFs = os.path.join(args.outdir, out_filename_TFs)          
        uniq_TFs = np.unique(lion['TF'])
        with open(save_file_path_TFs, 'w') as fp:
            fp.write('\n'.join(uniq_TFs))
        raise Exception(f"\n\nError: Not all specified TFs were found in the data frame. Please check that the names of the TFs were spelled correctly. A list of available TFs to filter for has been saved in {save_file_path_TFs}\n")
    
    col = lion.pop("Target")
    lion.insert(1, col.name, col)
    lion_subset = lion.loc[lion['TF'].isin(args.tfs)]
    
    

    lion_subset.to_csv(save_file_path, index = False)

    print(f"\nFile saved: {save_file_path}\n")  
