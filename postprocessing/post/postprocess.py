import pandas as pd
import numpy as np
import os
import pickle

def files_to_dfs(fname_panda, fname_lion, ftype):
    '''
    Creates data frames from the input panda and lioness files        
     
    Parameters:
    -----------
        - fname_panda: str, the panda file supplied by the user
        - fname_lioness: str, the lioness file supplied by the user
        - ftype: str, file type of lioness file, either npy or txt
    
    Returns:
    -----------
        - data frames of the panda and lioness files
    '''        
    # pandaFile = pd.read_csv(fname_panda, sep = " ", engine = "pyarrow", header = None)
    pandaFile = pd.read_csv(fname_panda, sep = " ", engine = "python", header = None)
    
    if ftype == "txt": 
        # lionFile = pd.read_csv(fname_lion, sep = "\t", engine = "pyarrow", header = None)
        lionFile = pd.read_csv(fname_lion, sep = "\t", engine = "python", header = None)
    elif ftype == "npy":    
        lionnpy = np.load(fname_lion)
        lionFile = pd.DataFrame(lionnpy)
        
    return([pandaFile,lionFile])

def extract_edges(df):
    '''
    Retrieves edge list from the lioness data frame
    
    Parameters:
    -----------
        - df: data frame with edge pairs as index (e.g. node1<==>node2) and columns as lioness edge weights
        
    Returns:
    --------
        - edge list in tuple format (e.g. [('TF_TF1', 'g1'), ('TF_TF2', 'g1'), etc.])
    '''    
    all_edges = list(df.index)
    
    un_all_edges = len(np.unique(all_edges))
    print(f"Number of unique edges in input: {un_all_edges}")
    
    edges_for_nx = []
    for edge in all_edges:
        splitedge = edge.split("<==>")
        edges_for_nx.append(splitedge)
        
    #print(edges_for_nx)
    tupled_edges = [tuple(l) for l in edges_for_nx]
    
    return(tupled_edges)

def add_weights(nw, df, edgelist, name):
    '''
    Adds edges to patient specific network by looping through the previous edge list and 
    assigning weights based on the lioness data frame
    
    Parameters
    ----------
        - nw: Network in nx format with unweighted edges
        - df: data frame with edge pairs as index (e.g. node1<==>node2) and columns as lioness edge weights
        - edgelist: List of edges extracted from the dataframe in extract_edges()
        - name: Name of sample to create network of
        
    Returns
    -------
        - a pickled dictionary where keys are sample names and values are the network for that sample 
    '''   
                    
    # Extract weights as list from sample in df      
    col_as_list = df[name].tolist()

    # Loop through column of weights and assign to edge
    i = 0
    while i < len(col_as_list):
        for edge in edgelist:
            nw[edge[0]][edge[1]]['weight'] = col_as_list[i]
            i += 1        
    
    return(nw)

def save_results(nw, format, base_name_begin, base_name_end, outdir):
    '''
    Saves file based on the specified user format
    
    Parameters
    ----------
        - nw: df containing the networks
        - format: format to save the file in, supplied by user. Either "csv" or "pickle"
        - base_name_begin: the beginning of the base name (pre-period) of the file to be output (example: lioness_nw)
        - base_name_end: the end of the base name (pre-period) of the file to be output (example: filtered)
        - outdir: path to the output directory
        
    Returns
    -------
        - a single pickle or csv file that is written out  
    ''' 
    
    if format == "csv":
        outfile = f"{base_name_begin}_{base_name_end}.csv"
        outfile_path = os.path.join(outdir, outfile)
        nw.to_csv(outfile_path, index = False)
    elif format == "pickle":
        outfile = f"{base_name_begin}_{base_name_end}.pickle"
        outfile_path = os.path.join(outdir, outfile)
        nw.to_pickle(outfile_path)   
    
    print(f"File saved: {outfile_path}")  
