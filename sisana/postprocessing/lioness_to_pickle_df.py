
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
                         
    # parser = argparse.ArgumentParser(description="Example command: python lioness_to_pickle_df.py -p panda_output.txt -q lioness_output.npy -t npy -n sampnames.txt -o ./output/lioness_df.pickle")
    # requiredArgGroup = parser.add_argument_group('Required arguments')  
    # requiredArgGroup.add_argument("-p", "--pandaFile", type=str, help="Path to panda output produced by the run_panda.py script", required=True) 
    # requiredArgGroup.add_argument("-q", "--lionessFile", type=str, help="Path to file produced by the run_lioness.py script", required=True)
    # requiredArgGroup.add_argument("-t", "--lionessFileType", type=str, choices = ['txt', 'npy'], help="File type of lioness input (the -q file)", required=True)
    # requiredArgGroup.add_argument("-n", "--sampnames", type=str, help="File with list of sample names (one per line) in the same order that were supplied to run_lioness.py", required=True)    
    # requiredArgGroup.add_argument("-o", "--outfile", default = "./lioness_output.pickle", type=str, help="Path to output file in pickle format (e.g. lioness.pickle)", required=True)    
    
    # args = parser.parse_args()
    
    # print("Reading in data...")

    # dfs_from_files = files_to_dfs(panda, lion, type)

    # Create data frames from input files
    pan = pd.read_csv(panda, sep = " ", engine = "python")
    pan["TF-gene"] = "TF_" + pan["tf"] + "<==>" + pan["gene"]


    # Lioness file does not have any header or column names, needs them for t-test later
    sampsfile = open(names, "r")
    fileread = sampsfile.read()

    lion.index = pan["TF-gene"]  

    savefile = outfile
    lion.to_pickle(savefile)   
 
    
    
