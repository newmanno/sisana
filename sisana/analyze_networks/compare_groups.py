import argparse
import os
import pandas as pd
import numpy as np
from pathlib import Path
import time
import csv
from scipy import stats
from .analyze import file_to_list, map_samples, calc_tt, calc_group_difference
import sys
# from statistics import mean

__author__ = 'Nolan Newman'
__contact__ = 'nolankn@uio.no'
    
def compare_bw_groups(datafile: str, mapfile: str, datatype: str, groups: list, testtype: str, filetype: str, outdir: str):
    '''
    Description:
        This code compares the degree or expression of two sample groups
     
    Parameters:
    -----------
        - datafile: str, Path to the data file
        - mapfile: str, Path to the mapping file, which maps sample name to sample group
        - datatype: str, The type of data being used ("expression" or "degree")
        - groups: str, Names of the two groups (from the second column of mapfile) to be compared. The second group listed will be used as the numerator in the fold change calulation.
        - testtype: str, Type of comparison to perform, either "tt" for Student's t-test, "mw" for Mann-Whitney U, "paired_tt", or "wilcoxon"
        - filetype: str, The type of data file ("csv" or "txt" or "tsv") being used
        - outdir: str, The directory to save the output to
        
    Returns:
    -----------
        - Nothing
    '''

    # parser = argparse.ArgumentParser(description="Example command: python compare_groups.py -m map.csv -p lioness.pickle -d degree -c high low -t mw -o ./output/")
    # ArgGroup = parser.add_argument_group('Required arguments')  
    # ArgGroup.add_argument("-p", "--datafile", type=str, help="Path to csv file containing the expression or degree (in/out) of each node", required=True)
    # ArgGroup.add_argument("-m", "--mapfile", type=str, help="Path to mapping file (csv). If doing an unpaired test (--testtype = mw or tt) then this file maps samples (first column, no header) to groups (second column, no header). Otherwise, if doing a paired analysis (--testtype = paired_tt or wilcoxon) then the samples for one group will go in column 1 (no header) while their paired samples will go in column 2.", required=True) 
    # ArgGroup.add_argument("-s", "--sampnames", type=str, help="Path to txt file containing the order of samples. Required if input file does not have a header", required=False)
    # ArgGroup.add_argument("-f", "--filetype", choices = ["csv", "txt"], type=str, help="Data file type of --datafile, either csv or txt (for tab separated files)", required=True)
    # ArgGroup.add_argument("-d", "--datatype", type=str, help="Type of input data, used as a suffix for the file names that are created", required=True)
    # ArgGroup.add_argument("-c", "--compgroups", type=str, nargs=2, help="Name of groups in mapping file to compare, required if not performing a paired analysis. Please note that if comparing expression, the second group listed will be used as the numerator in calculating the log2 fold change (e.g. log2(group2/group1))", required=False) 
    # ArgGroup.add_argument("-t", "--testtype", type=str, choices = ["tt", "mw", "paired_tt", "wilcoxon"], help="Type of comparison to perform, either Student's t-test, Mann-Whitney U, or a test for paired samples", required=True)     
    # ArgGroup.add_argument("-x", "--foldchange", action="store_true", help="Flag for whether to calculate the fold change between samples. Note: Do not use for degrees due to degrees having a negative metric", required=False)    
    # ArgGroup.add_argument("-o", "--outdir", type=str, help="Path to directory to output file to", required=True) 
    
    # args = parser.parse_args()
    
    # Create output directory if one does not already exist    
    os.makedirs(outdir, exist_ok=True)
    
    if filetype == "csv" and (datafile[-3:] == "txt" or datafile[-3:] == "tsv"):
        raise Exception("Error: The supplied data file has the extension 'csv', not the expected 'txt' or 'tsv' based on the supplied value to --filetype. If this is desired, please change the extension of the data file to match the extension of --filetype")
    if filetype == "txt" and datafile[-3:] == "csv":
        raise Exception(f"Error: The supplied data file has the extension 'txt', not the expected 'csv' based on the supplied value to --filetype. If this is desired, please change the extension of the data file to match the extension of --filetype")
    
    if filetype == "csv":
        datadf = pd.read_csv(datafile, index_col = 0)
    else:
        datadf = pd.read_csv(datafile, index_col = 0, sep = "\t")
        
    print("Performing calculations, please wait...")

    if testtype == "tt" or testtype == "mw":
          
        # Assign samples from mapping file to groups
        sampdict = map_samples(mapfile, groups[0], groups[1])
        total_samps = len(sampdict[groups[0]]) + len(sampdict[groups[1]])
    
    elif testtype == "paired_tt" or testtype == "wilcoxon":
        map = pd.read_csv(mapfile, header = None)

        groups = {}
        groups["group1"] = mapfile.iloc[:,0].tolist()
        groups["group2"] = mapfile.iloc[:,1].tolist()
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
    if testtype != "mw":
        if testtype == "tt": 
            pval = compdf.apply(lambda row : calc_tt(row[sampdict[groups[1]]], row[sampdict[groups[0]]], testtype), axis = 1)
        elif testtype == "paired_tt" or testtype == "wilcoxon":
            pval = compdf.apply(lambda row : calc_tt(row[groups["group1"]], row[groups["group2"]], testtype), axis = 1) 
            
        # Format the output data frame
        pval_column = testtype + "_pvalue"
        pvaldf = pd.DataFrame({'Target':pval.index, pval_column:pval.values})
        newpvaldf = pd.DataFrame(pvaldf[pval_column].to_list(), columns=['test_statistic',pval_column])
        newpvaldf['Target'] = pval.index
        newpvaldf = newpvaldf.set_index('Target')
        test_stat_column = "test_statistic"
         
    else:
        compdf.apply(lambda row : calc_tt(row[sampdict[groups[1]]], row[sampdict[groups[0]]], testtype), axis = 1)
        # mwu_uval, mwu_pval, mwu_cles = compdf.apply(lambda row : calc_tt(row[sampdict[groups[1]]], row[sampdict[groups[0]]], testtype), axis = 1)
        mwu_calculations = compdf.apply(lambda row : calc_tt(row[sampdict[groups[1]]], row[sampdict[groups[0]]], testtype), axis = 1)
    
        # Format the output data frame
        pval_column = testtype + "_pvalue"
        test_stat_column = "mw_uvalue"
        neglogp_column = "mw_-log(pvalue)"
        cles_column = "mw_CLES"
        
        pvaldf = pd.DataFrame({'Target':mwu_calculations.index, pval_column:mwu_calculations})
        newpvaldf = pd.DataFrame(pvaldf[pval_column].to_list(), columns=[test_stat_column, pval_column, neglogp_column, cles_column])
        newpvaldf['Target'] = mwu_calculations.index
        newpvaldf = newpvaldf.set_index('Target')

    # For expression, do regular log2FC
    # For degrees, first need to transform the edge value by doing ln(e^w + 1),
    # then calculate degrees. Then you can do the log2FC of degrees
    # Calculate log2FC (only for expression data, since you can have negative 
    # degree values). This transformation is described in the paper "Regulatory Network 
    # of PD1 Signaling Is Associated with Prognosis in Glioblastoma Multiforme"
    # if datatype == "expression":
    #     fc = compdf.apply(lambda row : calc_log2_fc(row[sampdict[groups[1]]], row[sampdict[groups[0]]]), axis = 1)

    mean_diff = compdf.apply(lambda row : calc_group_difference(row[sampdict[groups[0]]], row[sampdict[groups[1]]], difftype="mean"), axis = 1)
    median_diff = compdf.apply(lambda row : calc_group_difference(row[sampdict[groups[0]]], row[sampdict[groups[1]]], difftype="median"), axis = 1)

    print("Comparisons finished...") 
    
    # Calcuate means per group
    mean_g2_colname = f"mean_{groups[1]}"    
    mean_g1_colname = f"mean_{groups[0]}"      
    newpvaldf[mean_g2_colname] = compdf[sampdict[groups[1]]].mean(axis=1)
    newpvaldf[mean_g1_colname] = compdf[sampdict[groups[0]]].mean(axis=1)
    meandiff_colname = f"difference_of_means_({groups[1]}-{groups[0]})"      
    newpvaldf[meandiff_colname] = mean_diff
    newpvaldf["abs(difference_of_means)"] = abs(mean_diff)

    # Calcuate medians per group    
    median_g2_colname = f"median_{groups[1]}"    
    median_g1_colname = f"median_{groups[0]}"      
    newpvaldf[median_g2_colname] = compdf[sampdict[groups[1]]].median(axis=1)
    newpvaldf[median_g1_colname] = compdf[sampdict[groups[0]]].median(axis=1)    
    mediandiff_colname = f"difference_of_medians_({groups[1]}-{groups[0]})"      
    newpvaldf[mediandiff_colname] = median_diff
    newpvaldf["abs(difference_of_medians)"] = abs(median_diff)

    # Perform multiple test correction
    FDR_colname = "FDR"
    newpvaldf[FDR_colname] = stats.false_discovery_control(newpvaldf[pval_column])
    newpvaldf = newpvaldf.sort_values(pval_column, ascending = True)
    
    # Create new df without pval, ranked on test statistic
    if testtype != "mw":
        ranked = newpvaldf.sort_values(test_stat_column, ascending = False)
        ranked.drop([pval_column, FDR_colname], inplace=True, axis=1)
        ranked = ranked[test_stat_column]
    else:
        ranked = newpvaldf.sort_values(neglogp_column, ascending = False)
        ranked.drop([pval_column, FDR_colname], inplace=True, axis=1)
        ranked = ranked[neglogp_column]        
    
    # Rearrange column order so that FDR calculations comes after p-value
    if testtype != "mw":
        colorder = [test_stat_column, pval_column, FDR_colname,
                    mean_g2_colname, mean_g1_colname,
                    meandiff_colname, median_g2_colname, median_g1_colname,
                    mediandiff_colname]
        newpvaldf = newpvaldf.loc[:, colorder] 
    else:
        colorder = [test_stat_column, pval_column, neglogp_column, FDR_colname,
                    cles_column, mean_g2_colname, mean_g1_colname,
                    meandiff_colname, median_g2_colname, median_g1_colname,
                    mediandiff_colname]
        newpvaldf = newpvaldf.loc[:, colorder]
    
    # Write to disk
    if testtype == "tt" or testtype == "mw":
        save_file_path = os.path.join(outdir, f"comparison_{testtype}_between_{groups[0]}_{groups[1]}_{datatype}.txt")
        save_file_path_ranked = os.path.join(outdir, f"comparison_{testtype}_between_{groups[0]}_{groups[1]}_{datatype}_ranked_test_stat.rnk")
    if testtype == "paired_tt" or testtype == "wilcoxon":
        save_file_path = os.path.join(outdir, f"comparison_{testtype}_{datatype}.txt")
        save_file_path_ranked = os.path.join(outdir, f"comparison_{testtype}_{datatype}_ranked_test_stat.rnk")
    
    # Make output directory if it does not already exist
    Path(outdir).mkdir(parents=True, exist_ok=True)    
    
    # Remove the abs() columns (previously for sorting) for exporting the file
    columns_to_keep = [x for x in newpvaldf.columns if x[0:3] != "abs"]
    
    newpvaldf.to_csv(save_file_path, columns = columns_to_keep, sep = "\t")
    ranked.to_csv(save_file_path_ranked, sep = "\t", header = False)

    print(f"\nFile saved: {save_file_path}\nThis file contains all the statistics results and is just for your reference.\n")
    print(f"File saved: {save_file_path_ranked}\nThis file contains only the calculated test statistic for each gene and is used for input to the GSEA analysis step.\n")
