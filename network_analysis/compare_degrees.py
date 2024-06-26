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
from analyze import file_to_list, map_samples, calc_tt

__author__ = 'Nolan Newman'
__contact__ = 'nolankn@uio.no'
    
if __name__ == '__main__':
    """
    Description:
        This code compares the in-degrees/out-degrees from the output of lioness_df_indeg_outdeg_calculator.py. Comparisons are done between two sample groups.
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
    
    # Assign samples from mapping file to groups
    sampdict = map_samples(args.mapfile, groups[0], groups[1])
        
    # for k,v in sampdict.items():
    #     print(f"{k} samples: {v}")                     
        
    total_samps = len(sampdict[groups[0]]) + len(sampdict[groups[1]])
        
    # remove unnecessary samples to save on memory
    if len(liondf.columns) > total_samps:
        allsamps = []
        allsamps = sampdict[groups[0]] + sampdict[groups[1]]
        compdf = liondf.loc[:, allsamps]
    else:
        compdf = liondf
        
    del liondf
    
    # Add pseudocount of 1
    # print("\nhead")
    # print(compdf)
    # print("done")
    # compdf += 1
    
    # Calculate p-value/FDR
    pval = compdf.apply(lambda row : calc_tt(row[sampdict[groups[0]]], row[sampdict[groups[1]]], args.testtype), axis = 1)
    print("Comparisons finished...")
    # print(pval)
    
    pvaldf = pd.DataFrame({'Target':pval.index, 'pval':pval.values})
    newpvaldf = pd.DataFrame(pvaldf['pval'].to_list(), columns=['test_statistic','pval'])
    newpvaldf['Target'] = pvaldf['Target']
    newpvaldf = newpvaldf.set_index('Target')
    newpvaldf = newpvaldf.dropna()
    
    # Perform multiple test correction
    newpvaldf['FDR'] = stats.false_discovery_control(newpvaldf['pval'])
    
    # Create new df without pval, ranked on test statistic
    ranked = newpvaldf.sort_values('test_statistic', ascending = False)
    ranked.drop(['pval', 'FDR'], inplace=True, axis=1)
    # print(ranked)

    save_file_path = os.path.join(args.outdir, f"comparison_{args.testtype}_between_{args.compgroups[0]}_{args.compgroups[1]}.txt")
    save_file_path_ranked = os.path.join(args.outdir, f"comparison_{args.testtype}_between_{args.compgroups[0]}_{args.compgroups[1]}_ranked_test_stat.rnk")

    newpvaldf.to_csv(save_file_path, sep = "\t")
    ranked.to_csv(save_file_path_ranked, sep = "\t", header = False)

    print(f"\nFile saved: {save_file_path}\nThis file contains all the statistics results and is just for your reference.\n")
    print(f"File saved: {save_file_path_ranked}\nThis file contains only the calculated test statistic for each gene and is used for input to the GSEA analysis step.\n")
