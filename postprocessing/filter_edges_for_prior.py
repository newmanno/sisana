#test
import numpy as np
import pandas as pd
import argparse
import sys
import pickle
from pathlib import Path
from post import save_results
import time
import os 

if __name__ == '__main__':
    """
    Description:
        This code takes the edges that were output by Lioness and filters them for just the edges that were in the prior motif
    """
    
    parser = argparse.ArgumentParser(description="Example command: python filter_edges_for_prior.py -p lioness_df.pickle -m motif.tsv -o ./output -f pickle")
    ArgGroup = parser.add_argument_group('Required arguments')  
    
    ArgGroup.add_argument("-p", "--picklefile", type=str, help="Path to pickle file created by lioness_to_pickle_df.py script", required=True)
    ArgGroup.add_argument("-m", "--motiffile", type=str, help="Path to the prior motif file used in the original PANDA/LIONESS script", required=True)     
    ArgGroup.add_argument("-o", "--outdir", type=str, help="Path to directory to output file to", required=True) 
    ArgGroup.add_argument("-f", "--outformat", type=str, choices = ["csv", "pickle"], help="Format of file to output the filtered network to", required=True) 
    
    start = time.time()    

    args = parser.parse_args()
    nwdf = pd.read_pickle(args.picklefile)
    
    try: 
        motfile = pd.read_csv(args.motiffile, sep = "\t", engine = "pyarrow", header = None)
    except:
        raise Exception("There was an error reading in the data. Please make sure the motif file is tab-delimited.")
        sys.exit()
        
    print("Structure of the network data frame and prior motif file")
    print(nwdf.head())
    print(motfile.head())

    print(f"There are {len(nwdf)} edges in the network prior to filtering for prior edges")
    print(f"There are {len(motfile)} edges in the prior")

    # Extract a list of prior edges from the motif file
    prior_list = motfile.iloc[:,0]

    # Format the edge names to match the string format in the network file
    prior_list_TF_prefix = ["TF_" + i for i in prior_list]
    prior_list_matched_str_to_nw = [i.replace('-', '<==>') for i in prior_list_TF_prefix]


    #df_w_tf = nwdf[nwdf.index.str.contains('ARNT')]
    #print(df_w_tf)
    #mot_w_tf = nwdf[motfile.index.str.contains('ARNT')]
    #print(df_w_tf)

    #for i in prior_list_matched_str_to_nw:
    #    for j in list(nwdf.index):
    #        print(i, j)
    #        if i == j:
    #            print("Match")

    # Perform filtering
    #filt_nw = nwdf.loc[[prior_list_matched_str_to_nw],:]
    filt_nw = nwdf.loc[nwdf.index.isin(prior_list_matched_str_to_nw), :]
    #filt_nw = nwdf.loc[prior_list_matched_str_to_nw]
    print(filt_nw)

    print(f"There are {len(filt_nw)} edges in the network after filtering for prior edges")

    # format file names and output
    base_file_name = Path(args.picklefile).stem
    # outdeg_filename =  f"{base_file_name}_outdegree.csv"

    print("Now saving, please wait...")
    save_results(filt_nw, args.outformat, base_file_name, "filtered_for_prior_edges", args.outdir)   

    #testdf = nwdf.head()
    #print(testdf)
    #filtdata = ['TF_AR<==>41157', 'TF_ARNT<==>41157']

    #filt_nw = testdf.loc[filtdata]
    
    #print(filt_nw)
