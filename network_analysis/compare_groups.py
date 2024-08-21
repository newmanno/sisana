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
from analyze import file_to_list, map_samples, calc_tt, calc_log2_fc
import sys

__author__ = 'Nolan Newman'
__contact__ = 'nolankn@uio.no'
    
if __name__ == '__main__':
    """
    Description:
        This code compares either the in-degrees/out-degrees or expression of samples. Comparisons are done between two sample groups.
    """
    start = time.time()    

    parser = argparse.ArgumentParser(description="Example command: python compare_groups.py -m map.csv -p lioness.pickle -d degree -c high low -t mw -o ./output/")
    ArgGroup = parser.add_argument_group('Required arguments')  
    ArgGroup.add_argument("-m", "--mapfile", type=str, help="Path to mapping file (csv). If doing an unpaired test (--testtype = mw or tt) then this file maps samples (first column, no header) to groups (second column, no header). Otherwise, if doing a paired analysis (--testtype = paired_tt or wilcoxon) then the samples for one group will go in column 1 (no header) while their paired samples will go in column 2.", required=True) 
    ArgGroup.add_argument("-p", "--datafile", type=str, help="Path to csv file containing the degrees (in/out) of each node", required=True)
    ArgGroup.add_argument("-s", "--sampnames", type=str, help="Path to txt file containing the order of samples. Only required if using expression data, since degree data will already have the correct header", required=False)
    ArgGroup.add_argument("-f", "--filetype", choices = ["csv", "txt"], type=str, help="File type, either csv or txt (for tab separated files)", required=True)
    ArgGroup.add_argument("-d", "--datatype", choices = ["indegree", "outdegree", "expression"], type=str, help="Type of input data, either degree or expression", required=True)
    ArgGroup.add_argument("-c", "--compgroups", type=str, nargs=2, help="Name of groups in mapping file to compare, required if not performing a paired analysis. Please note that if comparing expression, the second group listed will be used as the numerator in calculating the log2 fold change (e.g. log2(group2/group1))", required=False) 
    ArgGroup.add_argument("-t", "--testtype", type=str, choices = ["tt", "mw", "paired_tt", "wilcoxon"], help="Type of comparison to perform, either Student's t-test, Mann-Whitney U, or a test for paired samples", required=True)     
    ArgGroup.add_argument("-o", "--outdir", type=str, help="Path to directory to output file to", required=True) 
    
    args = parser.parse_args()
    
    if args.filetype == "csv":
        datadf = pd.read_csv(args.datafile, index_col = 0)
    else:
        datadf = pd.read_csv(args.datafile, index_col = 0, sep = "\t")
        
    if args.datatype == "expression":        
        sampsfile = open(args.sampnames, "r")
        fileread = sampsfile.read()
        namelist = fileread.split("\n") 
        namelist = list(filter(None, namelist))
        
        datadf.columns = namelist 
    
    print(datadf.head())
    
    print("Performing calculations, please wait...")
    
    if args.testtype == "tt" or args.testtype == "mw":
        groups = args.compgroups
          
        # Assign samples from mapping file to groups
        sampdict = map_samples(args.mapfile, groups[0], groups[1])
        total_samps = len(sampdict[groups[0]]) + len(sampdict[groups[1]])
    
    elif args.testtype == "paired_tt" or args.testtype == "wilcoxon":
        map = pd.read_csv(args.mapfile, header = None)

        groups = {}
        groups["group1"] = map.iloc[:,0].tolist()
        groups["group2"] = map.iloc[:,1].tolist()
        total_samps = len(groups["group1"]) + len(groups["group2"])
               
    # remove unnecessary samples to save on memory
    if len(datadf.columns) > total_samps:
        allsamps = []
        allsamps = sampdict[groups[0]] + sampdict[groups[1]]
        compdf = datadf.loc[:, allsamps]
    else:
        compdf = datadf
        
    del datadf

    # Calculate p-value/FDR
    if args.testtype == "tt" or args.testtype == "mw":
        pval = compdf.apply(lambda row : calc_tt(row[sampdict[groups[0]]], row[sampdict[groups[1]]], args.testtype), axis = 1)
    elif args.testtype == "paired_tt" or args.testtype == "wilcoxon":
        # print("Comparing the following samples with a paired t-test:\n")
        # print(f"Samnples in group 1:\n{groups["group1"]})")
        # print(f"\nSamples in group 2:\n{groups["group2"]})")
        
        pval = compdf.apply(lambda row : calc_tt(row[groups["group1"]], row[groups["group2"]], args.testtype), axis = 1)
                
    # Calculate log2FC (only for expression data, since you can have negative degree values)
    if args.datatype == "expression":
        fc = compdf.apply(lambda row : calc_log2_fc(row[sampdict[groups[0]]], row[sampdict[groups[1]]]), axis = 1)                
    print("Comparisons finished...")
    
    # Format the output data frame
    pval_colname = args.testtype + "_pvalue"
    pvaldf = pd.DataFrame({'Target':pval.index, pval_colname:pval.values})
    newpvaldf = pd.DataFrame(pvaldf[pval_colname].to_list(), columns=['test_statistic',pval_colname])
    newpvaldf['Target'] = pval.index
    newpvaldf = newpvaldf.set_index('Target')  
    print(newpvaldf)    
    
    if args.datatype == "expression":
        fc_colname = f"mean_log2FC_{args.compgroups[1]}/{args.compgroups[0]}"
        newpvaldf[fc_colname] = fc
        
    # Perform multiple test correction
    newpvaldf['FDR'] = stats.false_discovery_control(newpvaldf[pval_colname])
    newpvaldf = newpvaldf.sort_values(pval_colname, ascending = True)
   
    # Create new df without pval, ranked on test statistic
    ranked = newpvaldf.sort_values('test_statistic', ascending = False)
    ranked.drop([pval_colname, 'FDR'], inplace=True, axis=1)
    
    # Write to disk
    if args.testtype == "tt" or args.testtype == "mw":
        save_file_path = os.path.join(args.outdir, f"comparison_{args.testtype}_between_{args.compgroups[0]}_{args.compgroups[1]}_{args.datatype}.txt")
        save_file_path_ranked = os.path.join(args.outdir, f"comparison_{args.testtype}_between_{args.compgroups[0]}_{args.compgroups[1]}_{args.datatype}_ranked_test_stat.rnk")
    if args.testtype == "paired_tt" or args.testtype == "wilcoxon":
        save_file_path = os.path.join(args.outdir, f"comparison_{args.testtype}_{args.datatype}.txt")
        save_file_path_ranked = os.path.join(args.outdir, f"comparison_{args.testtype}_{args.datatype}_ranked_test_stat.rnk")
    
    # Make output directory if it does not already exist
    Path(args.outdir).mkdir(parents=True, exist_ok=True)    
    
    newpvaldf.to_csv(save_file_path, sep = "\t")
    ranked.to_csv(save_file_path_ranked, sep = "\t", header = False)

    print(f"\nFile saved: {save_file_path}\nThis file contains all the statistics results and is just for your reference.\n")
    print(f"File saved: {save_file_path_ranked}\nThis file contains only the calculated test statistic for each gene and is used for input to the GSEA analysis step.\n")
