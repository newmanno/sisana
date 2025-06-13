import os
import pandas as pd
import pickle
import numpy as np
from pathlib import Path
import sys
from .analyze import indeg_calculator, outdeg_calculator

__author__ = 'Nolan Newman'
__contact__ = 'nolankn@uio.no'
    
def calculate_panda_degree(inputfile: str):
    '''
    Description:
        This code calculates the indegree and outdegree of genes and TFs, respectively, for PANDA networks
     
    Parameters:
    -----------
        - inputfile: str, Path to panda file in .txt format 
    
    Returns:
    -----------
        - Nothing
    '''

    panda_df = pd.read_csv(inputfile, sep=" ")   
    
    # Ensure columns of factors are strings
    panda_df['tf'] = panda_df['tf'].astype(str)
    panda_df['gene'] = panda_df['gene'].astype(str)
    
    # Remove unneeded columns or else it adds strings together and wastes memory
    panda_df_tf = panda_df.drop(columns=['gene', 'motif'])
    panda_df_target = panda_df.drop(columns=['tf', 'motif'])
    
    # Perform calculation
    outdeg = outdeg_calculator(panda_df_tf)
    indeg = indeg_calculator(panda_df_target)
    
    # format file names and output
    outdeg_filename =  f"{str(inputfile)[:-4]}_outdegree.csv"
    indeg_filename =  f"{str(inputfile)[:-4]}_indegree.csv"

    outdeg.to_csv(outdeg_filename, index_label='tf')  
    indeg.to_csv(indeg_filename)

def calculate_lioness_degree(inputfile: str, datatype: str):
    '''
    Description:
        This code calculates the indegree and outdegree of genes and TFs, respectively, for LIONESS networks
     
    Parameters:
    -----------
        - inputfile: str, Path to lioness file, either in .csv format or the .pickle file created by lioness_to_pickle_df.py script
        - datatype: str, Type of inputfile, either "csv" or "pickle"
    
    Returns:
    -----------
        - Nothing
    '''

    if datatype == 'pickle':
        nwdf = pd.read_pickle(inputfile)
    elif datatype == 'csv':
        nwdf = pd.read_csv(inputfile, index_col='TF-target')       
    
    # Separate the row names to be two different columns, TF and Target
    nwdf = nwdf.rename_axis("TF").reset_index()
    nwdf[['tf','gene']] = nwdf['TF'].str.split('<==>',expand=True)
    
    # Ensure columns of factors are strings
    nwdf['tf'] = nwdf['tf'].astype(str)
    nwdf['gene'] = nwdf['gene'].astype(str)
    
    # Remove unneeded columns or else it adds strings together and wastes memory
    nwdf_tf = nwdf.drop(columns=['TF', 'gene'])
    nwdf_target = nwdf.drop(columns=['TF', 'tf'])

    # Perform calculation
    outdeg = outdeg_calculator(nwdf_tf)
    indeg = indeg_calculator(nwdf_target)

    # Remove the "TF_" prefix for TFs, which was used for PANDA/LIONESS calculations 
    newind = [x[3:] for x in outdeg.index]
    outdeg.index = newind
    
    # format file names and output
    file_path_and_stem = os.path.splitext(inputfile)[0]
    outdeg_filename =  f"{file_path_and_stem}_outdegree.csv"
    indeg_filename =  f"{file_path_and_stem}_indegree.csv"

    file_path_and_stem = os.path.splitext(inputfile)[0]

    outdeg.to_csv(outdeg_filename, index_label='tf')  
    indeg.to_csv(indeg_filename, index_label='target')
    
