
import numpy as np
import pandas as pd
import argparse
from .post import files_to_dfs

def convert_lion_to_pickle(panda: str, lion: str, type: str, names: str, outfile: str):    
    '''
    Creates data frames from the input panda and lioness files        
     
    Parameters:
    -----------
        - panda: str, Path to panda output file
        - lion: str, lioness data frame with transformed edge values
        - type: str, file type of lioness file, either npy or txt
        - names: str, File with list of sample names (one per line) in the same order that were supplied for panda/lioness
        - outfile: str, Path to output file in pickle format (e.g. lioness.pickle)
    
    Returns:
    -----------
        - Nothing
    '''
    # Create data frames from input files
    pan = pd.read_csv(panda, sep = " ", engine = "python")
    pan["TF-gene"] = "TF_" + pan["tf"] + "<==>" + pan["gene"]

    # Lioness file does not have any header or column names, needs them for t-test later    
    samps = []
    with open(names, 'r') as file:
        for line in file: 
            samps.append(line.strip())

    lion.index = pan["TF-gene"]  
    lion.columns = samps 
    
    savefile = outfile
    lion.to_pickle(savefile)   
 
    
    
