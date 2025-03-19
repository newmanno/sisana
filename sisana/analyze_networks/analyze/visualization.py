import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import warnings

def filter_for_user_defined_genes(datafile: pd.DataFrame, genes: list, verbose: bool=True):
    '''
    Description:
        This code filters a data file for the genes defined by a user
     
    Parameters:
    -----------
        - datafile: str, Data frame containing the expression or indegrees of each gene (rows) per sample (column)
        - genes: list, List of genes that the user wishes to fileter for
    
    Returns:
    -----------
        - The input data frame, filtered for the genes supplied by the user
    '''
    datafile = datafile.T
    filtered_datafile = datafile.loc[:,genes]
    filtered_datafile = filtered_datafile.T
    
    if verbose:
        print(f"\nData file BEFORE filtering for user-defined genes")
        print(datafile.T)

        print(f"\nData file AFTER filtering for user-defined genes")
        print(filtered_datafile)
    
    return(filtered_datafile)

def filter_for_top_genes(datafile: pd.DataFrame, statsfile: pd.DataFrame, number: int, pval_threshold: float=0.05, returnvalue: str="df", verbose: bool=True):
    '''
    Description:
        This code filters a data file for significant genes, then further sorts based on effect size
        and filters the data frame for the top X largest effect sizes that are also significant
     
    Parameters:
    -----------
        - datafile: str, Data frame containing the expression or indegrees of each gene (rows) per sample (column)
        - statsfile: str, Data frame containing the comparison statistics between two groups
        - number: int, Number of top genes to find.
        - pval_threshold: float, The value to use as a p-value threshold
        - returnvalue: str, Which type of data to return, either "df" for the filtered data frame, or "genelist" for a 
          list of the top genes
    
    Returns:
    -----------
        - The input data frame, filtered for the top genes
    '''
    datafile = datafile.T

    if number==None:  
        number=10

    # Filter out any non-significant genes
    signif_indata = statsfile[statsfile['FDR'] < pval_threshold]
    
    sorted_indata = signif_indata.sort_values('FDR', ascending=True)

    topgenes = sorted_indata.index[0:number]

    try:
        filtered_datafile = datafile.loc[:,topgenes]
        filtered_datafile = filtered_datafile.T

    except KeyError:
        topgenes = list(set(datafile.columns) & set(topgenes))
        filtered_datafile = datafile.loc[:,topgenes]
        filtered_datafile = filtered_datafile.T

    if verbose:
        print(f"\nData file BEFORE filtering for top genes")
        print(datafile.T)

        print(f"\nData file AFTER filtering for top genes")
        print(filtered_datafile)

    if returnvalue == "df":
        return(filtered_datafile)
    else:
        return(filtered_datafile.index)
    
def create_label_list(statsfile: pd.DataFrame, difference_column: str, adjp_column: str):
    '''
    Description:
        This code creates a list of labels in the format [Text(value, value, label), etc.], which is necessary for 
        drawing the labels on the volcano plots
     
    Parameters:
    -----------
        - statsfile: str, Data frame containing the gene name as index and statistics for plotting as columns
        - difference_column: str, Name of the column containing the difference (to be plotted on the x-axis)
        - adjp_column: str, Name of the column containing the adjust p-value/FDR (to be plotted on the y-axis)
    
    Returns:
    -----------
        - A list of labels in the format [Text(value, value, label), etc.]
    '''
    label_list = []
    for i,r in statsfile.iterrows():
        label_list.append(plt.text(x=r[difference_column],y=-np.log10(r[adjp_column]),s=i))
    return(label_list)

def find_sample_overlap(samp_data_list, samp_meta_list):
    '''
    Description:
        This code Finds the overlap between the samples in the data df and the metadata df. Note that ideally there 
        should be a perfect overlap, but this ensures the code will not crash due to issues later 
        
    Parameters:
    -----------
        - samp_data_list: str, Data frame containing the gene name as index and statistics for plotting as columns
        - samp_meta_list: str, Name of the column containing the difference (to be plotted on the x-axis)
    
    Returns:
    -----------
        - A list of samples that are common between the two data frames
    '''
    dif1 = np.setdiff1d(np.array(samp_data_list), np.array(samp_meta_list))
    dif2 = np.setdiff1d(np.array(samp_meta_list), np.array(samp_data_list))
    
    diffsamps = np.concatenate((dif1, dif2))    
        
    if len(diffsamps) > 0:  
        warnings.warn("Warning: Your samples in your metadata file are not a perfect match to the samples in the data file. The heatmaps will only be made on the overlapping samples of the two.") 
        print(f"\nThe following samples are DIFFERENT between the metadata and data table, and thus will not be plotted: {diffsamps}\n")
    
    samp_overlaps = list(set(samp_data_list) & set(samp_meta_list)) # Get the overlapping samples between metadata and data df        
    return(samp_overlaps)