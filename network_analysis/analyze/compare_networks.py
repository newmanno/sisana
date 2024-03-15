import scipy 
from scipy import stats
import csv
import re
import pandas 
import math

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
    elif ttype == 'mw':
        p = stats.mannwhitneyu(group1, group2)
    print("\n")
    print(scipy.mean(group1))
    print(scipy.mean(group2))

    if scipy.mean(group1) == 0 and scipy.mean(group2) == 0:
        log2FC = 0
    elif scipy.mean(group1) != 0:

    elif scipy.mean(group2) != 0:
        log2FC = math.log2(scipy.mean(group1)/scipy.mean(group2))
    elif scipy.mean(group2) == 0:
        print(group2)
        # Take the smallest non-zero value, divide by 10, and use that as the average value for the denominator
        no_zeros = [i for i in group1 if i != 0]
        denom = min(no_zeros)/10
        log2FC = math.log2(scipy.mean(group1)/denom)        
    else:
        log2FC = "NA"

    return([p[1], log2FC])
