
import numpy as np
from netZooPy.lioness import Lioness
from netZooPy.lioness.analyze_lioness import AnalyzeLioness
import pandas as pd
import argparse
import sys
import pyarrow.csv as arrcsv
import pickle

if __name__ == '__main__':
    """
    Description:
        This code converts the output of the lioness file to a pickle format
    """

    def files_to_dfs(fname_panda, fname_lion, ftype):
        '''
        Creates data frames from the input panda and lioness files         
        Arguments:
            - fname_panda: str, the panda file supplied by the user
            - fname_lioness: str, the lioness file supplied by the user
            - ftype: str, file type of lioness file, either npy or txt
        '''        
        try: 
            pandaFile = pd.read_csv(fname_panda, sep = " ", engine = "pyarrow", header = None)
            
            if ftype == "txt": 
                lionFile = pd.read_csv(fname_lion, sep = "\t", engine = "pyarrow", header = None)
            elif ftype == "npy":    
                lionnpy = np.load(fname_lion)
                lionFile = pd.DataFrame(lionnpy)
                print(lionFile)
        except:
            raise Exception("There was an error reading in the data. Please make sure the file paths are correct and the lioness data is in the correct format you specified.")
            sys.exit()
            
        return([pandaFile,lionFile])
        
    parser = argparse.ArgumentParser(description="Example command: python initial_analysis_exp.py -p panda_output.txt -q lioness_output.npy -t npy -o ./output/lioness_df.pickle")
    requiredArgGroup = parser.add_argument_group('Required arguments')  
    requiredArgGroup.add_argument("-p", "--pandaFile", type=str, help="Path to panda output produced by the run_panda.py script", required=True) 
    requiredArgGroup.add_argument("-q", "--lionessFile", type=str, help="Path to file produced by the run_lioness.py script", required=True)
    requiredArgGroup.add_argument("-t", "--loinessFileType", type=str, choices = ['txt', 'npy'], help="File type of lioness input (the -q file)", required=True)
    requiredArgGroup.add_argument("-o", "--outfile", default = "./lioness_output.pickle", type=str, help="Path to output file in pickle format (e.g. lioness.pickle)", required=False)    
    
    args = parser.parse_args()
    
    print("Reading in data...")

    dfs_from_files = files_to_dfs(args.pandaFile, args.lionessFile, args.loinessFileType)

    # Create data frames from input files
    pan = dfs_from_files[0]
    lion = dfs_from_files[1]

    pan.columns = ['TF', 'Target', 'Interaction', 'Score']
    pan["TF-target"] = "TF_" + pan["TF"] + "<==>" + pan["Target"]

    # Lioness file does not have any header or column names, needs them for t-test later
    samps = ["samp" + str(i) for i in range(1, len(lion.columns)+1)]

    lion.columns = samps
    lion.index = pan["TF-target"]  
    
    savefile = args.outfile
    lion.to_pickle(savefile)   
    print(f"File saved: {savefile}")
 
    
    
