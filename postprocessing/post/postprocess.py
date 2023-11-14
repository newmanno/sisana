import pandas as pd
import numpy as np
import pyarrow.csv as arrcsv

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