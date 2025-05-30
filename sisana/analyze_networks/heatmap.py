import seaborn as sns 
import pandas as pd
import sys
import pickle
from pathlib import Path
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.stats import zscore
import matplotlib.pyplot as plt
import warnings
from .analyze import file_to_list, filter_for_top_genes, filter_for_user_defined_genes, find_sample_overlap
import csv
import os
import warnings
import numpy as np
import scipy

def plot_heatmap(datafile: str, filetype: str, statsfile: str, metadata: str, genelist: str, 
                 groups: list, prefix: str, plotnames: bool,
                 outdir: str, top: bool=True, additional_metadata_categories: list=[]):
    '''
    Description:
        This code creates a heatmap of either the expression or degrees from LIONESS networks
     
    Parameters:
    -----------
        - datafile: str, Path to file containing the expression or indegrees of each gene per sample
        - filetype: str, Type of inputfile, either "csv" for comma separated files or "txt" or "tsv" for tab-delimited
        - statsfile: str, Path to the file that contains the comparison output of the gene expression or degree
        - metadata: str, Path to the csv metadata file mapping samples to groups (groups must match names of the groups arg), must have a header of the format 'name,group'
        - genelist: str, .txt file containing a list of genes to plot
        - groups: list, The names of the groups to plot (from the --metadata file). Groups will be plotted as column groups in the order they are written in this argument
        - prefix: str, Prefix to use for the output file; note that the output file will automatically be generated with the suffix '_heatmap.png'
        - plotnmaes: str, Flag for whether to plot the names of the genes on the heatmap        
        - outdir: str, Path to output directory
        - top: Flag for whether to automatically plot the top 10 genes. Does not use the genelist in this case, but rather finds the top genes
               based on FDR and fold change.
    
    Returns:
    -----------
        - Nothing
    '''
    # parser = argparse.ArgumentParser(description="Example command: python heatmap.py -d <indegree.csv> -t csv -f -m <metadata.csv> -g <genelist.txt> -n group1 group2 group3 -o ./output/")
    # ArgGroup = parser.add_argument_group('Required arguments')  
    
    # ArgGroup.add_argument("-d", "--datafile", type=str, help="Path to file containing the expression or indegrees of each gene per sample", required=True)
    # ArgGroup.add_argument("-t", "--filetype", choices = ["csv", "txt"], help="Type of delimiter used for --datafile", required=True)    
    # # ArgGroup.add_argument("-f", "--fileheader", action='store_true', help="Flag for if --datafile has a header already. If not, requires --sampleorder", required=False)    
    # # ArgGroup.add_argument("-s", "--sampleorder", type=str, help=".txt file containing the order of the samples in the data file", required=False)   
    # ArgGroup.add_argument("-m", "--metadata", type=str, help="Path to the csv metadata file mapping samples to groups (groups must match names of the --groupnames arg), must have a header of the format 'name,group'", required=True) 
    # ArgGroup.add_argument("-g", "--genelist", type=str, help=".txt file containing a list of genes to plot", required=True)   
    # ArgGroup.add_argument("-c", "--hierarchicalcluster", action='store_true', help="Flag for if you wish to perform hierarchical clustering on the genes", required = False)    
    # ArgGroup.add_argument("-n", "--groupnames", type=str, nargs = "+", help="The names of the groups to plot (from the --metadata file). Groups will be plotted as column groups in the order they are written in this argument", required=True)   
    # ArgGroup.add_argument("-p", "--prefix", type=str, help="Prefix to use for the output file; note that the output file will automatically be generated with the suffix '_heatmap.png'", required=True)   
    # ArgGroup.add_argument("-o", "--outdir", type=str, help="Path to directory to output file to", required=True) 
    
    # args = parser.parse_args()
    
    if filetype == "csv":
        datadf = pd.read_csv(datafile, index_col = 0)
    elif filetype == "txt" or filetype == "tsv":
        datadf = pd.read_csv(datafile, index_col = 0, sep = "\t")
        
    os.makedirs(outdir, exist_ok=True)
     
    def _assign_samples(mapfile, selected_groups):
        '''
        Function that assigns samples to groups for statistical analysis

            Arguments:
                - mapfile: pandas dataframe of the mapping file from the user
                - selected_groups: the list of groups the user wants to plot in the heatmap
        '''
    
        samp_type_dict = {}
        
        # Create empty lists inside of a dictionary keyed on group names to store sample names as values
        for i in selected_groups:
            samp_type_dict[i] = []
        
        # Add all node-type pairs from the input file into the node_type_dict
        for index, row in mapfile.iterrows():
            samp_type_dict[row.iloc[1]].append(row.iloc[0])

        return(samp_type_dict)
  
    # Add column names if they are not already there
    # Note: Removed after changing upstream steps to ensure that only headered files
    # can be used in the pipeline
    # if not fileheader:
    #     samps = file_to_list(sampleorder)
    #     indata.columns = samps     


    # Find the overlap between the samples in the data df and the metadata df
    samp_data_list = datadf.columns
    samp_meta_file = pd.read_csv(metadata)
    samp_meta_list = samp_meta_file.iloc[:,0].tolist()    
    sample_overlap = find_sample_overlap(samp_data_list, samp_meta_list)
  
    dat = datadf[datadf.columns.intersection(sample_overlap)] # remove non-overlap samples from data df
    samp_meta_file = samp_meta_file[samp_meta_file.iloc[:,0].isin(sample_overlap)] # remove non-overlap samples from metadata df
    
    # Assign samples from mapping file to groups
    sampdict = _assign_samples(samp_meta_file, groups)
    print("Number of samples in each group:")
    for k,v in sampdict.items():
        print(f"{k}: {len(v)}")
    
    compare_df = pd.read_csv(statsfile, sep = "\t", index_col=0)
    
    # If the user has supplied the gene list, plot just the genes they supplied. Otherwise, plot the top genes based on FDR and fold change
    if not top:
        # Filter the data frame containing the degree/expression values for just the genes in the supplied gene list
        genes_to_plot = file_to_list(genelist)
        filtered = dat.filter(items = genes_to_plot, axis=0)

    else:
        filtered = filter_for_top_genes(datafile=dat, 
                       statsfile=compare_df,                     
                       number=50)
   
    out_filtered_dat_path = os.path.join(outdir, f"{prefix}_filtered_data_file_for_heatmap_genes.csv")
    filtered.to_csv(out_filtered_dat_path)

    # Calculate z-score. Normally this could just be done in the call to the heatmap
    # method, but in this case we need grouped heatmaps with hierarchical clustering,
    # so I am calculating z-scores first, then performing the clustering, and 
    # finally appending heatmaps to one another on a per-group basis
    filtered_z = filtered.apply(zscore, axis=1)
    Z = linkage(filtered_z, method='ward')
    leaf_order = dendrogram(Z, no_plot=True)['ivl']
    ordered_df = filtered_z.iloc[map(int, leaf_order),:]
    
    out_filtered_z_path = os.path.join(outdir, f"{prefix}_filtered_data_file_for_heatmap_genes_zscore.csv")
    filtered_z.to_csv(out_filtered_z_path)    
    ordered_df.to_csv(os.path.join(outdir, f"{prefix}_filtered_data_file_for_heatmap_genes_ordered_df.csv"))   


    # Split the z-score data frame per sample group
    zdfs = {}
    all_input_samps = ordered_df.columns
    
    for i in groups:
        zdfs[i] = ordered_df[sampdict[i]]
        
    print(zdfs)
    sys.exit(0)
    
    # fig, axs =plt.subplots(nrows = 1, ncols = len(groups), constrained_layout = True)
    fig, axs = plt.subplots(nrows = 1, ncols = len(groups), figsize=(len(groups) * 2, 2), width_ratios=[len(df.columns) for df in zdfs.values()], constrained_layout = True)
        
    print("\nCreating heatmap(s), please wait...")
    
    # Manually set the min and max values of the heatmaps so the same scale applies to all heatmaps,
    # otherwise the color bar only applies to the last heatmap that gets plotted
    arr = ordered_df.to_numpy().flatten() # an array of all values in the dat_z df

    min_val = np.percentile(arr, 5) # vals will be the 1st and 99th percentile, so extreme (outlier) values will not influence the scale
    max_val = np.percentile(arr, 95)
    
    # min_val = ordered_df.to_numpy().min()
    # max_val = ordered_df.to_numpy().max()
    
    cmap_col = "RdBu_r"
    square_choice = False
    xticklabels_choice = True
    yticklabels_choice = plotnames
    center_choice = 0

    i = 0 # Counter to keep track of position in axs array
    # sys.exit(0)
    if len(groups) < 2:
        ax = axs

        sns_plot = sns.heatmap(zdfs[groups[0]], ax = axs, xticklabels = xticklabels_choice, yticklabels = yticklabels_choice, square = square_choice, cmap = cmap_col, cbar = False, vmin = -3, vmax = 3, center = center_choice)
        sns_plot.set_yticklabels(sns_plot.get_yticklabels(), rotation = 0, fontsize = 2)
        sns_plot.set_xticklabels(sns_plot.get_xticklabels(), fontsize = 2)
        sns_plot.set_ylabel('',size=6)
        ax.yaxis.set_tick_params(labelsize = 2)

        outname = os.path.join(outdir, f"{prefix}_heatmap.png")
        plt.savefig(outname, dpi = 600)

    else:
        
        # Parse through the zdfs dict (keyed on group name, values are dfs) and make a single heatmap plot of multiple heatmap subplots
        for k,v in zdfs.items():
            ax = axs[i] 

            col_options_counter = 0

            # Need to create each heatmap separately depending on its position in the 1D array of heatmaps (axs)
            if i == 0:
                
                sns_plot = sns.heatmap(v, ax=axs[i], xticklabels = xticklabels_choice, yticklabels = yticklabels_choice, square = square_choice, cmap = cmap_col, cbar = False, vmin = -3, vmax = 3, center = center_choice)
                sns_plot.set_yticklabels(sns_plot.get_yticklabels(), rotation = 0, fontsize = 2)
                sns_plot.set_xticklabels(sns_plot.get_xticklabels(), fontsize = 2)
                sns_plot.set_ylabel('',size=6)
                ax.yaxis.set_tick_params(labelsize = 2)

            # for any middle position (not first or last) in the axs array
            elif (i > 0) and (i < len(groups) - 1):
                sns_plot = sns.heatmap(v, ax=axs[i], xticklabels= xticklabels_choice, square = square_choice, cmap = cmap_col, cbar = False, yticklabels = False, vmin = -3, vmax = 3, center = center_choice) 
                sns_plot.set_ylabel('')    
                sns_plot.set_xticklabels(sns_plot.get_xticklabels(), fontsize = 2)
            
            # for the last position in the axs array
            else:            
                sns_plot = sns.heatmap(v, ax=axs[i], xticklabels= xticklabels_choice, square = square_choice, cmap = cmap_col, cbar = True, yticklabels = False, vmin = -3, vmax = 3, center = center_choice, cbar_kws={'label': 'z-score', "shrink": 0.5}) 
                sns_plot.set_xticklabels(sns_plot.get_xticklabels(), fontsize = 2)
                sns_plot.set_ylabel('') 
                ax.figure.axes[-1].yaxis.label.set_size(5)
                ax.figure.axes[-1].tick_params(labelsize=3, width=0.3, length=1)

            # import matplotlib
            # cmap = matplotlib.colormaps['RdBu_r']
            # norm = matplotlib.colors.Normalize(min_val, max_val)
            # cbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm, cmap), ax)
            
            ax.set_title(k, fontsize = 8)
            ax.tick_params(axis='both', which='both', length=0)

            # ax.figure.tight_layout()
            
            i += 1
        
    plt.yticks(rotation=0)  
   
    outname = os.path.join(outdir, f"{prefix}_heatmap.png")
    plt.savefig(outname, dpi = 600)
    
    print(f"\nFile created: {outname}")
    print(f"File created: {out_filtered_dat_path}")
    print(f"File created: {out_filtered_z_path}")
