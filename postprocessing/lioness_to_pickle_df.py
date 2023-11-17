
import numpy as np
import pandas as pd
import argparse
import sys
import pickle
from post import files_to_dfs

if __name__ == '__main__':
    """
    Description:
        This code converts the output of the lioness file to a pickle format and adds names of the samples to the df
    """
        
    parser = argparse.ArgumentParser(description="Example command: python lioness_to_pickle_df.py -p panda_output.txt -q lioness_output.npy -t npy -o ./output/lioness_df.pickle")
    requiredArgGroup = parser.add_argument_group('Required arguments')  
    requiredArgGroup.add_argument("-p", "--pandaFile", type=str, help="Path to panda output produced by the run_panda.py script", required=True) 
    requiredArgGroup.add_argument("-q", "--lionessFile", type=str, help="Path to file produced by the run_lioness.py script", required=True)
    requiredArgGroup.add_argument("-t", "--lionessFileType", type=str, choices = ['txt', 'npy'], help="File type of lioness input (the -q file)", required=True)
    requiredArgGroup.add_argument("-n", "--sampnames", type=str, help="File with list of sample names (one per line) in the same order that were supplied to run_lioness.py", required=True)    
    requiredArgGroup.add_argument("-o", "--outfile", default = "./lioness_output.pickle", type=str, help="Path to output file in pickle format (e.g. lioness.pickle)", required=True)    
    
    args = parser.parse_args()
    
    print("Reading in data...")

    dfs_from_files = files_to_dfs(args.pandaFile, args.lionessFile, args.lionessFileType)

    # Create data frames from input files
    pan = dfs_from_files[0]
    lion = dfs_from_files[1]

    pan.columns = ['TF', 'Target', 'Interaction', 'Score']
    pan["TF-target"] = f"TF_{pan["TF"]}<==>{pan["Target"]}

    # Lioness file does not have any header or column names, needs them for t-test later
    sampsfile = open(args.sampnames, "r")
    fileread = sampsfile.read()
    namelist = fileread.split("\n") 
    namelist = list(filter(None, namelist))
        
    lion.columns = namelist
    lion.index = pan["TF-target"]  
    
    savefile = args.outfile
    lion.to_pickle(savefile)   
    print(f"File saved: {savefile}")
 
    
    
