import scipy 
from scipy import stats
import scipy.stats
import csv
import re
import pandas 
from statistics import mean
import math
import numpy as np
import sys

def file_to_list(fname):
    """
    This function takes as input a text file with one object per line and returns the contents of that file as a list

    Args:
        fname: the name of the file to convert to a list

    Returns:
        returnlist: list of objects from the text file
    """
    all_lines = open(fname, "r").read().splitlines()
    returnlist = [name for name in all_lines if name]

    return (returnlist)   

def map_samples(mapfile, type1, type2):
    '''
    Function that assigns samples to groups for statistical analysis

        Arguments:
            - mapfile: input file from user
            - type1: the name of the group in the first set of samples, must be present in the mapfile
            - type2: the name of the group in the second set of samples, must be present in the mapfile
    '''

    samp_type_dict = {}

    type1_list = []
    type2_list = []

    init_dict = {}
    samp_type_dict = {}

    # Add all node-type pairs from the input file into the node_type_dict
    with open(mapfile) as samp_file:
        samp_file = csv.reader(samp_file, delimiter = ',')

        for row in samp_file:
            init_dict[row[0]] = row[1]

    for key,value in init_dict.items():
        try:
            if re.search(type1, value):
                type1_list.append(key)
            elif re.match(type2, value):
                type2_list.append(key)
        except:
            print("Unexpected value found for the sample group name.")

    samp_type_dict[type1] = type1_list
    samp_type_dict[type2] = type2_list

    return (samp_type_dict)

def calc_tt(group1, group2, ttype):
    '''
    Performs either a students t-test or mann-whitney test between two groups
    
        Arguments:
            - group1: list of samples in first group
            - group2: list of samples in second group
            - ttype: str, the type of test to perform, either tt or mw
    '''
    
    if ttype == 'tt':
        p = stats.ttest_ind(group1, group2)       
        statistic = p[0]
        pval = p[1]
    elif ttype == 'mw':
        p = stats.mannwhitneyu(group1, group2)
        
        # Note that since mann whitney technically has two test statistics (one for U1 and the other for U2), we have to compute both and take the minimum.
        # Wilcoxon also has two test statistics, but by default it returns the smaller of the two, so no need for SiSaNA to calculateNo
        u1 = p[0]
        u2 = group1.shape[0] * group2.shape[0] - p[0]
        
        statistic = min(u1,u2)
        if mean(group1) < mean(group2):
            statistic = statistic * -1

        pval = p[1]
    elif ttype == 'paired_tt':
        p = stats.ttest_rel(group1, group2)
        statistic = p[0]
        pval = p[1]
    elif ttype == 'wilcoxon':
        p = stats.wilcoxon(group1, group2)
        statistic = p[0]
        pval = p[1]

    return(statistic, pval)
    
def transform_edge_to_positive_val(edgeval):
    '''
    Transform degree values to be positive so fold change can be calculated. This transformation 
    is described in the paper "Regulatory Network of PD1 Signaling Is Associated with Prognosis 
    in Glioblastoma Multiforme"
        
        Arguments:
            - edgeval: float, value of edge
            
        Returns:
            - float, transformed value of edge
    '''
    # We get an inf if we have too large of values, but since transforming a large value
    # does not change the resulting value (i.e. ln(e^1000) + 1) = 1000, we can avoid the
    # inf values by just not transforming those edges
    if edgeval > 700:
        newval = edgeval
    else:
        newval = math.log(np.exp(edgeval) + 1)
    if newval == float("inf"):
        raise Exception(f"val {edgeval} gives a result of inf")
    return(newval)

def calc_log2_fc(group1, group2):
    '''
    Calculates the log2 fold change of means across two groups
        
        Arguments:
            - group1: list of samples in first group
            - group2: list of samples in second group
    '''
    
    # print("\n")
    # print(f"group1 values:")
    # print(group1)
    # print(f"Mean of group 1: {scipy.mean(group1)}")
    # print(f"group2 values:")
    # print(group2)
    # print(f"Mean of group 2: {scipy.mean(group2)}")

    ### Note: For the following calculations, assuming the user has followed the previous SiSaNA steps,
    ### neither group will ever have a mean of 0 since the preprocess.py step filters out any genes of 
    ### low abundance. So to simplify the logic for this step, I have removed any checks for group means
    ### of 0. May need to implement this change later though if people are not running preprocess.py.

    # if scipy.mean(group1) == 0 and scipy.mean(group2) == 0:
    #     log2FC = 0
    # elif scipy.mean(group1) == 0:
        
    # elif scipy.mean(group2) != 0:
    
    avg_g1 = mean(group1)
    avg_g2 = mean(group2)

    if avg_g1 == 0 or avg_g2 == 0:
        log2FC = "NA"
        # print("log2FC is NA")

    elif all(i >= 0 for i in group1) and all(i >= 0 for i in group2):
        log2FC = math.log2(avg_g2/avg_g1)
        # print(f"log2FC is {log2FC}\n")

        # print(log2FC)
    else:
        print(group1, group2)
        raise Exception("\n\nError: Negative values found in data. The log2 fold change can only be calculated on expression data. Negative values indicate that degree was likely used as an input type instead.\n")
    # elif scipy.mean(group2) == 0:
    #     # Take the smallest non-zero value, divide by 10, and use that as the average value for the denominator
    #     no_zeros = [i for i in group1 if i != 0]
    #     denom = min(no_zeros)/10
    #     log2FC = math.log2(scipy.mean(group1)/denom)        
    # else:
    #     log2FC = "NA"
    
    return(log2FC)
