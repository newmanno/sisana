import networkx as nx
import argparse
import os
import pandas as pd
import pickle
import numpy as np
from pathlib import Path
import time

   
__author__ = 'Nolan Newman'
__contact__ = 'nolankn@uio.no'
    
if __name__ == '__main__':
    """
    Description:
        This code converts the lioness outputs to individual networkx networks
    """

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

    parser = argparse.ArgumentParser(description="Example command: python create_single_samp_nws.py -e <file.pickle> -o ./all_nws.pickle")
    ArgGroup = parser.add_argument_group('Required arguments')  
    ArgGroup.add_argument("-p", "--picklefile", type=str, help="Path to pickle file created by lioness_to_pickle_df.py script", required=True) 
    ArgGroup.add_argument("-o", "--outdir", type=str, help="Path to directory to output file to", required=True) 
    
    args = parser.parse_args()
    
    start = time.time()    

    lion = pd.read_pickle(args.picklefile)    
    end = time.time()
    timetaken = end - start
    print(f"Time required for import: {timetaken}")
    
    # Convert the edge list from the data frame to be accepted into networkx
    print("Converting to edges for nx...")
    edges_from_rnames = extract_edges(lion)
    test = edges_from_rnames[0:4]
    
    start = time.time()    
    print("Creating nx networks...")
    G = nx.DiGraph()     # One graph is originally created with all the edges that each sample/network has (since they all have same edges, but different weights)
    G.add_edges_from(list(edges_from_rnames))
    
    listsamps = lion.columns
    all_lion_nws = {}

    # Add edge weights for each sample
    for sampname in listsamps:
        H = G.copy()
        H = add_weights(H, lion, edges_from_rnames, sampname)
        all_lion_nws[sampname] = H

    end = time.time()
    timetaken = end - start
    print(f"Time required for converting nws to nx format: {timetaken}")

    # write out a pickle file so the creation step doesnt need to be repeated
    base_file_name = Path(args.picklefile).stem
    filename =  f"{base_file_name}_nx_nws.pickle"
    save_file_path = os.path.join(args.outdir, filename)
   
    start = time.time()    
    pickle.dump(all_lion_nws, open(save_file_path, "wb"))
    print(f"File saved: {save_file_path}")
    
    end = time.time()
    timetaken = end - start
    print(f"Time required for writing pickle file: {timetaken}")
    