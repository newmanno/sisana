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
    
    def assign_node_type(node_list_file, type1, type2):
        '''
        Function that assigns samples to groups for statistical analysis
        
            Arguments:
                - node_list_file: input file from user
                - type1: the type of nodes in group 1
                - type2: the type of nodes in group 2
        '''    

        samp_type_dict = {}

        type1_list = []
        type2_list = []            

        init_dict = {}
        samp_type_dict = {}

        # Add all node-type pairs from the input file into the node_type_dict
        with open(node_list_file) as node_file:
            node_file = csv.reader(node_file, delimiter = ',')
            
            for row in node_file:
                init_dict[row[0]] = row[1]
                        
        for key,value in init_dict.items():
            try:
                if re.search(type1, value):
                    type1_list.append(key)
                elif re.match(type2, value):
                    type2_list.append(key)
            except:
                print("Unexpected value in the 'type' column of node_type input file.")    
                   
        samp_type_dict[type1] = type1_list
        samp_type_dict[type2] = type2_list

        return (samp_type_dict)        
        
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
        
    print(compdf)
    
    def calc_tt(group1, group2, ttype):
        if ttype == 'tt':
            p = stats.ttest_ind(group1, group2)
        elif ttype == 'mw':
            p = stats.mannwhitneyu(group1, group2)

        return(p[1])
        
    print("Original DataFrame:\n", compdf)
    print(sampdict)
    compdf['pval'] = compdf.apply(lambda row : calc_tt(row[sampdict['high']], row[sampdict['low']], args.testtype), axis = 1)
    compdf['FDR'] = stats.false_discovery_control(compdf['pval'])
    print(compdf)