import networkx as nx
import argparse
import os
import pandas as pd
import pickle
import numpy as np
import time
   
__author__ = 'Nolan Newman'
__contact__ = 'nolankn@uio.no'
    
if __name__ == '__main__':
    """
    Description:
        This code calculates network properties of the lioness networks
    """

    def indegree_calculation(nw_dict):
        '''
        Calculates the in-degree for each network in the input dictionary
        
        Parameters:
        -----------
            - nw_dict: dictionary keyed on sample name with sample networkx network as values
            
        Returns:
        --------
            - ??? Not sure yet, will update later
        '''    
        
        print("Starting in-degree calculation...")
        
        for k,v in nw_dict.items():
            print(f"Network {k}")
            for node in list(all_nws[k].nodes):
                if not node.startswith("TF_"): # Do not calculate in-degree for TFs
                    indeg = all_nws[k].in_degree(node, weight = 'weight')
                    #print(f"{node} in-degree: {indeg}")
            print("\n")
    
    def outdegree_calculation(nw_dict):
        '''
        Calculates the out-degree for each network in the input dictionary
        
        Parameters:
        -----------
            - nw_dict: dictionary keyed on sample name with sample networkx network as values
            
        Returns:
        --------
            - ??? Not sure yet, will update later
        '''    
        print("Starting out-degree calculation...")

        for k,v in nw_dict.items():
            print(f"Network {k}")
            for node in list(all_nws[k].nodes):
                if node.startswith("TF_"): # Only calculate out-degree for TFs
                    indeg = all_nws[k].out_degree(node, weight = 'weight')
                    #print(f"{node} out-degree: {indeg}")
            print("\n")

    parser = argparse.ArgumentParser(description="Example command: python analyze_networks.py -e <file.pickle> -o ./file.txt")
    ArgGroup = parser.add_argument_group('Required arguments')  
    ArgGroup.add_argument("-p", "--picklefile", type=str, help="Path to pickle file created by create_single_samp_nws.py script", required=True) 
    ArgGroup.add_argument("-o", "--outputfile", default = "./lioness_all_single_samp_nws.pickle", type=str, help="Path to pickle file to be output, contains all lioness networks", required=False) 
    
    args = parser.parse_args()
    
    print("Reading in pickle file... please wait...")
    start = time.time()    
    all_nws = pd.read_pickle(args.picklefile)
    end = time.time()
    timetaken = end - start
    print(f"Time required for import: {timetaken}")
    
    indegree_calculation(all_nws)
    outdegree_calculation(all_nws)

    end = time.time()
    timetaken = end - start
    print(f"Time required for whole script: {timetaken}")    