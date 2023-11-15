import argparse
import os
import pandas as pd
import pickle
import numpy as np
from pathlib import Path
import time
import csv
import re
from scipy import stats
from analyze import assign_node_type, calc_tt

__author__ = 'Nolan Newman'
__contact__ = 'nolankn@uio.no'
    
if __name__ == '__main__':
    """
    Description:
        This code calculates the in-degree and out-degree of lioness networks in a pickled df object
    """
    start = time.time()    

    parser = argparse.ArgumentParser(description="Example command: python compare_degrees.py -m map.csv -p lioness.pickle -c high low -o ./output")
    ArgGroup = parser.add_argument_group('Required arguments')  
    ArgGroup.add_argument("-m", "--mapfile", type=str, help="Path to mapping file (csv) that maps samples (first column, no header) to groups (second column, no header)", required=True) 
    ArgGroup.add_argument("-p", "--degfile", type=str, help="Path to csv file containing the degrees (in/out) of each node", required=True) 
    ArgGroup.add_argument("-c", "--compgroups", type=str, nargs=2, help="Name of groups in mapping file to compare", required=True) 
    ArgGroup.add_argument("-t", "--testtype", type=str, choices = ["tt", "mw"], help="Type of comparison to perform, either Student's t-test or Mann-Whitney U", required=True)     
    ArgGroup.add_argument("-o", "--outdir", type=str, help="Path to directory to output file to", required=True) 
    
    args = parser.parse_args()
    liondf = pd.read_csv(args.degfile, index_col = 0)
    groups = args.compgroups
    
    sampdict = assign_node_type(args.mapfile, groups[0], groups[1])
    
    for k,v in sampdict.items():
        print(f"{k} samples: {v}")                     
        
    total_samps = len(sampdict[groups[0]]) + len(sampdict[groups[1]])
    
    if len(liondf.columns) > total_samps:
        allsamps = []
        allsamps = sampdict[groups[0]] + sampdict[groups[1]]
        compdf = liondf.loc[:, allsamps]
    else:
        compdf = liondf
        
    del liondf
            
    compdf['pval'] = compdf.apply(lambda row : calc_tt(row[sampdict[groups[0]]], row[sampdict[groups[1]]], args.testtype), axis = 1)
    compdf['FDR'] = stats.false_discovery_control(compdf['pval'])
    
    save_file_path = os.path.join(args.outdir, f"comparison_{args.testtype}_between_{args.compgroups[0]}_{args.compgroups[1]}.txt")
    compdf.to_csv(save_file_path, sep = "\t")

    print(f"\nFile saved: {save_file_path}\n")
