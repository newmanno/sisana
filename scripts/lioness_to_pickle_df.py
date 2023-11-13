
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
        This code converts the output of the .txt lioness file to a pickle format
    """

    def files_to_dfs(fname_panda, fname_lion, ftype):
        
        try: 
            pandaFile = pd.read_csv(fname_panda, sep = " ", engine = "pyarrow", header = None)
            
            if ftype == "txt": 
                lionFile = pd.read_csv(fname_lion, sep = "\t", engine = "pyarrow", header = None)
            elif ftype == "npy":    
                lionnpy = np.load(fname_lion)
                lionFile = pd.DataFrame(lionnpy)
                print(lionFile)
        except:
            print("There was an error reading in the data. Please make sure the file paths are correct and the data is in the correct format you specified.")
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
    
    pan = dfs_from_files[0]
    lion = dfs_from_files[1]

    pan.columns = ['TF', 'Target', 'Interaction', 'Score']
    
    pan["TF-target"] = "TF_" + pan["TF"] + "<==>" + pan["Target"]
        
    samps = ["samp" + str(i) for i in range(1, len(lion.columns)+1)]
    print(samps)
    lion.columns = samps
    lion.index = pan["TF-target"]  
    print(lion)
    
    savefile = args.outfile
    
    lion.to_pickle(savefile)   
    
    print(f"File saved: {savefile}")
 
    
    
