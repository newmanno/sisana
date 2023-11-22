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
    requiredArgGroup.add_argument("-t", "--type", type=str, choices = ["tf", "gene"], help="Do you want to subset for TFs (tf) or genes (target)?", required=True)     
    requiredArgGroup.add_argument("-g", "--genes", type=str, nargs = '+', help="Name of genes to subset for; required if gene is used for --type (-t)", required=False)     
    requiredArgGroup.add_argument("-f", "--tfs", type=str, nargs = '+', help="Name of TFs to subset for; required if tf is used for --type (-t)", required=False)     
    requiredArgGroup.add_argument("-o", "--outdir", type=str, help="Path to output directory", required=True)    
    
    args = parser.parse_args()
    lion = pd.read_pickle(args.picklefile) 
    
    base_file_name = Path(args.picklefile).stem
    if args.type == "tf":
        out_filename =  f"{base_file_name}_filtered_for_TFs.csv"
    elif args.type == "gene":    
        out_filename =  f"{base_file_name}_filtered_for_genes.csv"
    else:
        raise Exception('Error: Please provide either "tf" or "gene" as your desired node type to subset for')

    save_file_path = os.path.join(args.outdir, out_filename)   
    
    # Separate the row names to be two different columns, TF and Target
    lion = lion.rename_axis("TF").reset_index()
    lion[['TF','Target']] = lion['TF'].str.split('<==>',expand=True)   
    lion['TF'] = lion['TF'].str.replace('TF_', '')

    # Check if the input TF(s)/gene(s) is/are in the data frame. If not, give the user a list of the available ones to filter for
    if args.type == "tf":
        if set(args.tfs).issubset(np.unique(lion['TF'])):
            print("Specified TFs are all in the lioness data frame, now filtering...")
        else:
            # Output a list of available TFs to filter file for
            out_filename_TFs =  f"{base_file_name}_available_TFs.csv"
            save_file_path_TFs = os.path.join(args.outdir, out_filename_TFs)          
            uniq_TFs = np.unique(lion['TF'])
            with open(save_file_path_TFs, 'w') as fp:
                fp.write('\n'.join(uniq_TFs))
            raise Exception(f"\n\nError: Not all specified TFs were found in the data frame. Please check that the names of the TFs were spelled correctly. A list of available TFs to filter for has been saved in {save_file_path_TFs}\n")
    elif args.type == "gene":
        if set(args.genes).issubset(np.unique(lion['Target'])):
            print("Specified genes are all in the lioness data frame, now filtering...")
        else:
            # Output a list of available TFs to filter file for
            out_filename_genes =  f"{base_file_name}_available_genes.csv"
            save_file_path_genes = os.path.join(args.outdir, out_filename_genes)          
            uniq_genes = np.unique(lion['Target'])
            with open(save_file_path_genes, 'w') as fp:
                fp.write('\n'.join(uniq_genes))
            raise Exception(f"\n\nError: Not all specified genes were found in the data frame. Please check that the names of the genes were spelled correctly. A list of available genes to filter for has been saved in {save_file_path_genes}\n")
    
    # move the "Target" column to be the second column in the df
    col = lion.pop("Target")
    lion.insert(1, col.name, col)
    
    # Subset for desired node type and save results
    if args.type == "tf":
        lion_subset = lion.loc[lion['TF'].isin(args.tfs)]
        lion_subset.to_csv(save_file_path, index = False)
    elif args.type == "gene":
        lion_subset = lion.loc[lion['Target'].isin(args.genes)]
        lion_subset.to_csv(save_file_path, index = False)

    print(f"\nFile saved: {save_file_path}\n")  

