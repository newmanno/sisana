import argparse
import os
import pandas as pd
import pickle
import numpy as np
from pathlib import Path
import time
from .analyze import indeg_calculator, outdeg_calculator

__author__ = 'Nolan Newman'
__contact__ = 'nolankn@uio.no'
    
def calculate_degree(inputfile: str, datatype: str, outdir: str):
    '''
    Description:
        This code calculates the indegree and outdegree of genes and TFs, respectively
     
    Parameters:
    -----------
        - inputfile: str, Path to lioness file, either in .csv format or the .pickle file created by lioness_to_pickle_df.py script
        - datatype: str, Type of inputfile, either "csv" or "pickle"
        - outdir: str, Path to output directory
    
    Returns:
    -----------
        - Nothing
    '''
    # start = time.time()    

    # parser = argparse.ArgumentParser(description="Example command: python lioness_df_indeg_outdeg_calculator.py -i lioness_df.pickle -t pickle -o ./output")
    # ArgGroup = parser.add_argument_group('Required arguments')  
    # ArgGroup.add_argument("-i", "--inputfile", type=str, help="Path to lioness file, either in .csv format or the .pickle file created by lioness_to_pickle_df.py script", required=True) 
    # ArgGroup.add_argument("-t", "--inputtype", choices= ['csv', 'pickle'], type=str, help="File type of the lioness file", required=True) 
    # ArgGroup.add_argument("-o", "--outdir", type=str, help="Path to directory to output file to", required=True) 
    
    # args = parser.parse_args()
    
    if datatype == 'pickle':
        nwdf = pd.read_pickle(inputfile)
    elif datatype == 'csv':
        nwdf = pd.read_csv(inputfile, index_col='TF-target')       
    
    # Separate the row names to be two different columns, TF and Target
    nwdf = nwdf.rename_axis("TF").reset_index()
    nwdf[['TF','Target']] = nwdf['TF'].str.split('<==>',expand=True)
    
    # Ensure columns of factors are strings
    nwdf['TF'] = nwdf['TF'].astype(str)
    nwdf['Target'] = nwdf['Target'].astype(str)
    
    # Remove unneeded columns or else it adds strings together and wastes memory
    nwdf_tf = nwdf.drop(columns=['Target'])
    nwdf_target = nwdf.drop(columns=['TF'])



    # Perform calculation
    outdeg = outdeg_calculator(nwdf_tf)
    indeg = indeg_calculator(nwdf_target)

    # Remove the "TF_" prefix for TFs, which was used for PANDA/LIONESS calculations 
    newind = [x[3:] for x in outdeg.index]
    outdeg.index = newind
    
    # format file names and output
    base_file_name = Path(inputfile).stem
    outdeg_filename =  f"{base_file_name}_outdegree.csv"
    indeg_filename =  f"{base_file_name}_indegree.csv"

    save_file_path_outdeg = os.path.join(outdir, outdeg_filename)
    save_file_path_indeg = os.path.join(outdir, indeg_filename)
    
    outdeg.to_csv(save_file_path_outdeg, index_label='TF')  
    indeg.to_csv(save_file_path_indeg)
    
    print("\nFinished calculating degrees!")
    print(f"In-degree output can be found here: {save_file_path_indeg}")  
    print(f"Out-degree output can be found here: {save_file_path_outdeg}")  

    # end = time.time()
    # timetaken = end - start
    # print(f"\nTime required for script: {timetaken:2.2f} seconds")
