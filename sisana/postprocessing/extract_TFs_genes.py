import argparse
import sys
import pickle
import pandas as pd
import numpy as np
from pathlib import Path
import os

#Get access to file_to_list function that is in another directory
# userdir = os.getcwd() # Use this to return to the user's directory after importing the function

# scriptPath = os.path.realpath(os.path.dirname(sys.argv[0]))
# print(scriptPath)
# os.chdir(scriptPath)
# sys.path.append("../network_analysis/")
# from .analyze import file_to_list

def file_to_list(fname):
    all_lines = open(fname, "r").read().splitlines()
    returnlist = [name for name in all_lines if name]
    return(returnlist)

def extract_tfs_genes(pickle: str, datatype: str, namefile: str, outdir: str):
    '''
    Description:
        This code filters the lioness output file for only desired TFs or genes
     
    Parameters:
    -----------
        - pickle: str, Path to lioness output in pickle format
        - datatype: str, Do you want to subset for TFs (tf) or genes (gene)?
        - namefile: str, Single column text file, no header, containing the name of genes or TFs (must be consistent with --type) to subset for
        - outdir: str, Path to output directory
    
    Returns:
    -----------
        - Nothing
    '''
    
    
    # os.chdir(userdir)
    
    # parser = argparse.ArgumentParser(description="Example command: python extract_TFs_genes.py -p lioness_output.pickle -t gene -n genes.txt -o ./output/")
    # requiredArgGroup = parser.add_argument_group('Required arguments')  
    # requiredArgGroup.add_argument("-p", "--picklefile", type=str, help="Path to lioness output in pickle format", required=True) 
    # requiredArgGroup.add_argument("-t", "--datatype", type=str, choices = ["tf", "gene"], help="Do you want to subset for TFs (tf) or genes (target)?", required=True)     
    # requiredArgGroup.add_argument("-n", "--namefile", type=str, help="Single column text file, no header, containing the name of genes or TFs (must be consistent with --type) to subset for", required=True)     
    # requiredArgGroup.add_argument("-o", "--outdir", type=str, help="Path to output directory", required=True)    
    
    # args = parser.parse_args()
    
    # Create output dir if one does not already exist
    os.makedirs(outdir, exist_ok=True)
    
    lion = pd.read_pickle(pickle)
    
    base_file_name = Path(pickle).stem
    if datatype == "tfs":
        out_filename =  f"{base_file_name}_filtered_for_TFs.csv"
    elif datatype == "genes":    
        out_filename =  f"{base_file_name}_filtered_for_genes.csv"

    save_file_path = os.path.join(outdir, out_filename)   
    
    # Separate the row names to be two different columns, TF and Target
    lion = lion.rename_axis("TF").reset_index()
    lion[['TF','Target']] = lion['TF'].str.split('<==>',expand=True)   
    lion['TF'] = lion['TF'].str.replace('TF_', '')
    
    # Import the list of TFs/genes to subset for
    names_to_find = file_to_list(namefile)

    # Check if the input TF(s)/gene(s) is/are in the data frame. If not, give the user a list of the available ones to filter for
    if datatype == "tfs":
        if set(names_to_find).issubset(np.unique(lion['TF'])):
            print("Specified TFs are all in the lioness data frame, now filtering...")
        else:
            # Output a list of available TFs to filter file for
            out_filename_TFs =  f"{base_file_name}_available_TFs.csv"
            save_file_path_TFs = os.path.join(outdir, out_filename_TFs)          
            uniq_TFs = np.unique(lion['TF'])
            with open(save_file_path_TFs, 'w') as fp:
                fp.write('\n'.join(uniq_TFs))
            raise Exception(f"\n\nError: Not all specified TFs were found in the data frame. Please check that the names of the TFs were spelled correctly. A list of available TFs to filter for has been saved in {save_file_path_TFs}\n")
    elif datatype == "genes":
        if set(names_to_find).issubset(np.unique(lion['Target'])):

            print("Specified genes are all in the lioness data frame, now filtering...")
        else:
            # Output a list of available TFs to filter file for
            out_filename_genes =  f"{base_file_name}_available_genes.csv"
            save_file_path_genes = os.path.join(outdir, out_filename_genes)          
            uniq_genes = np.unique(lion['Target'])
            with open(save_file_path_genes, 'w') as fp:
                fp.write('\n'.join(uniq_genes))
            raise Exception(f"\n\nError: Not all specified genes were found in the data frame. Please check that the names of the genes were spelled correctly. A list of available genes to filter for has been saved in {save_file_path_genes}\n")
    
    # move the "Target" column to be the second column in the df
    col = lion.pop("Target")
    lion.insert(1, col.name, col)
    
    # Subset for desired node type and save results
    if datatype == "tfs":
        lion_subset = lion.loc[lion['TF'].isin(names_to_find)]
        lion_subset.to_csv(save_file_path, index = False)
    elif datatype == "genes":
        lion_subset = lion.loc[lion['Target'].isin(names_to_find)]
        lion_subset.to_csv(save_file_path, index = False)

    print(f"\nFile saved: {save_file_path}\n")  

