import argparse
import os
import pandas as pd
import numpy as np
from pathlib import Path
import csv
from scipy import stats
import sys
from sklearn.metrics.cluster import adjusted_rand_score

__author__ = 'Nolan Newman'
__contact__ = 'nolankn@uio.no'
    
if __name__ == '__main__':
    """
    Description:
        This code compares the clustering found for samples using the adjusted rand score.
    """

    parser = argparse.ArgumentParser(description="Example command: python compare_clusters.py -i map.csv -o ./output/")
    ArgGroup = parser.add_argument_group('Required arguments')  
    ArgGroup.add_argument("-i", "--infile", type=str, help="Path to input file (csv). File must have a header, which is used to specify the columns used for comparison.", required=True)   
    ArgGroup.add_argument("-c", "--compcols", type=str, nargs = 2, help="Column names to use for comparison.", required=True)   
    
    ArgGroup.add_argument("-o", "--outdir", type=str, help="Path to directory to output file to", required=True) 
    
    args = parser.parse_args()
    
    indata = pd.read_csv(args.infile, engine = "python", index_col=[0])

    group1_assignment = list(indata[args.compcols[0]])
    group2_assignment = list(indata[args.compcols[1]])
    
    ar_score = adjusted_rand_score(group1_assignment, group2_assignment)
    
    print(f"Adjusted rand: {ar_score:.3f}")