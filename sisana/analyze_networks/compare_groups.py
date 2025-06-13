import os
import pandas as pd
import numpy as np
from pathlib import Path
import csv
from scipy import stats
from .analyze import file_to_list, map_samples, calc_tt, calc_group_difference
import sys
from numpy import log

__author__ = 'Nolan Newman'
__contact__ = 'nolankn@uio.no'
    
def compare_bw_groups(datafile: str, mapfile: str, datatype: str, groups: list, testtype: str, filetype: str, rankby_col: str, outdir: str):
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
        - rankby_col: str, Choices: ["mediandiff", "mwu", "neglogp", "meandiff"]. The statistic to rank the .rnk output file by for GSEA. 
        - outdir: str, The directory to save the output to
        
    Returns:
    -----------
        - list of the output file paths
    '''
    
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
        
    if testtype == "tt" or testtype == "mw":
          
        # Assign samples from mapping file to groups
        mapfile = pd.read_csv(mapfile, index_col=0)
        sampdict = map_samples(mapfile, groups[0], groups[1])
        total_samps = len(sampdict[groups[0]]) + len(sampdict[groups[1]])
    
    elif testtype == "paired_tt" or testtype == "wilcoxon":
        mapfile = pd.read_csv(mapfile)

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
        
    print("Performing comparisons, please wait...")
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
        neglogp_column = "mw_signed_-log(pvalue)"
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
    
    if testtype == "mw": 
        if rankby_col == "mwu":
            sortcol = "mw_uvalue"
        elif rankby_col == "mediandiff":
            sortcol = f"difference_of_medians_({groups[1]}-{groups[0]})"
        elif rankby_col == "meandiff":
            sortcol = f"difference_of_means_({groups[1]}-{groups[0]})"
        elif rankby_col == "neglogp":
            sortcol = "mw_signed_-log(pvalue)"
    else:
        sortcol = test_stat_column
    
    # Create new df without pval, ranked on test statistic (as chosen by user)
    ranked = newpvaldf.sort_values(sortcol, ascending = False)
    ranked.drop([pval_column, FDR_colname], inplace=True, axis=1)
    ranked = ranked[sortcol]
    
    
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
        
        if testtype == "tt":
            save_file_path_ranked = os.path.join(outdir, f"comparison_{testtype}_between_{groups[0]}_{groups[1]}_{datatype}_ranked_test_stat.rnk")
        else:
            save_file_path_ranked = os.path.join(outdir, f"comparison_{testtype}_between_{groups[0]}_{groups[1]}_{datatype}_ranked_{rankby_col}.rnk")
            
    if testtype == "paired_tt" or testtype == "wilcoxon":
        save_file_path = os.path.join(outdir, f"comparison_{testtype}_{datatype}.txt")
        save_file_path_ranked = os.path.join(outdir, f"comparison_{testtype}_{datatype}_ranked_{rankby_col}.rnk")
    
    # Make output directory if it does not already exist
    Path(outdir).mkdir(parents=True, exist_ok=True)    
    
    # Remove the abs() columns (previously for sorting) for exporting the file
    columns_to_keep = [x for x in newpvaldf.columns if x[0:3] != "abs"]
    
    newpvaldf.to_csv(save_file_path, columns = columns_to_keep, sep = "\t")
    ranked.to_csv(save_file_path_ranked, sep = "\t", header = False)

    print(f"\nFile saved: {save_file_path}\nThis file contains all the statistics results and is just for your reference.\n")
    print(f"File saved: {save_file_path_ranked}\nThis file contains only the calculated test statistic for each gene and is used for input to the GSEA analysis step.\n")
    
    return([save_file_path, save_file_path_ranked])
