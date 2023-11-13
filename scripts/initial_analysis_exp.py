# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 16:23:26 2023

@author: Nolan K Newman

This script takes as input an expression matrix and finds the number of samples
each gene is expressed in, as well as their stdev and variance
"""
import numpy as np
import pandas as pd
import argparse

if __name__ == '__main__':
    """
    Description:
        Filters expression data for parameters (e.g. genes) that are only present in at least m samples.
    """

    parser = argparse.ArgumentParser(description="Example command: python initial_analysis_exp.py -e <expression_file.tsv>")
    requiredArgGroup = parser.add_argument_group('Required arguments')        
    requiredArgGroup.add_argument("-e", "--exp", type=str, help="Path to file containing the gene expression data. By default, the expression file does not have a header, and the cells ares separated by a tab", required=True)
    requiredArgGroup.add_argument("-m", "--minsamps", type=int, help="Minimum number of samples a gene must be expressed in; expression data will be filtered for only genes that pass", required=True)
    requiredArgGroup.add_argument("-o", "--out-file", type=str, dest = "outfile", help="prefix to output file (i.e. '/path/to/prefix', which gives /path/to/prefix.txt)", required=True)    
    
    args = parser.parse_args()
    
    exp_file = args.exp
    outdir = args.outfile
    
    expdf = pd.read_csv(exp_file, sep='\t', header=None, index_col=0)
    num_originalgenes = len(expdf)
    nsamps = len(expdf.columns)
    header_samps = ["sample_"+str(i) for i in list(range(1,nsamps + 1))] # Create the headers for samples
    
    expdf["stdev"] = expdf.std(axis=1) # Calculate the standard deviation for each gene

    num_nonzeros = nsamps - np.count_nonzero(expdf.iloc[:,0:nsamps], axis=1) # Count the number of non zeros for each gene in the df
    expdf["num_zeros"] = num_nonzeros # add the count to a new col in the df
    
    # Subset df for only genes that appear in at least k samples (user-specified)
    #print(expdf.iloc[:,[len(expdf.columns)-1]] < args.minsamps)
    samples_removed = len(expdf[expdf["num_zeros"] >= args.minsamps])
    
    cutdf = expdf[expdf["num_zeros"] < args.minsamps]
    cutdf = cutdf.drop(columns=['stdev', 'num_zeros'])
    cutdf.columns = header_samps
    
    # Print summary statistics 
    dist = expdf["num_zeros"].value_counts()
    dist.index.name = None
    print("\n\nNumber of samples with zero expression | Number of instances")
    print(dist)    
    
    print("\n" + str(samples_removed) + " out of " + str(num_originalgenes) + " parameters were removed that did not pass the minsamps threshold. " + str(len(cutdf)) + " parameters remain. Results have been stored in " + outdir + "_filtered.txt\n")
    
    # Export
    expdf.to_csv(outdir + "_w_summary_stats.txt", sep = "\t")
    cutdf.to_csv(outdir + "_filtered_low_abundant.txt", sep = "\t", columns = header_samps, header = None)
    
    # Now filter the motif prior for only the TFs and target genes that are in the resulting filtered expression file
    
    # Then filter the ppi data for the same 
    print(min(cutdf["sample_1"]))