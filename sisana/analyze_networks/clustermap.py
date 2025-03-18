from .analyze import file_to_list, filter_for_top_genes, filter_for_user_defined_genes, find_sample_overlap
import seaborn as sns 
import pandas as pd
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib 
import warnings
from .analyze import WrongAmountOfColorsError
def pnq(obj):
        print(obj)
        sys.exit(0)
      
def plot_clustermap(datafile: str, filetype: str, metadata: str, genelist: str, column_cluster: bool,
                 row_cluster: bool, prefix: str, outdir: str, plot_gene_names: bool, plot_sample_names: bool, top: bool=True, 
                 category_label_columns: list=[], category_column_colors: list=[], statsfile: str=""):
    '''
    Description:
        This code creates a heatmap of either the expression or degrees from LIONESS networks
     
    Parameters:
    -----------
        - datafile: str, Path to file containing the expression or indegrees of each gene per sample
        - filetype: str, Type of inputfile, either "csv" for comma separated files or "txt" or "tsv" for tab-delimited
        - statsfile: str, Path to the file that contains the comparison output of the gene expression or degree
        - metadata: str, Path to the csv metadata file mapping samples to groups (groups must match names of the groups arg)
            Example:
                name group metastasis  stage
                TCGA_E2_A10E_01A_21R_A10J_07  LumA        yes      4
                TCGA_D8_A1JB_01A_11R_A13Q_07  LumA        yes      4
                TCGA_AR_A5QN_01A_12R_A28M_07  LumB        yes      3
                TCGA_LL_A50Y_01A_11R_A266_07  LumB        no      1
                TCGA_B6_A0I5_01A_11R_A034_07  LumB        no      2
            
        - genelist: str, .txt file containing a list of genes to plot. Required if "top" is not set
        - rowcluster: bool, Flag for if you wish to cluster the rows
        - prefix: str, Prefix to use for the output file; note that the output file will automatically be generated with the suffix '_clustermap.png'
        - outdir: str, Path to output directory
        - top: Flag for whether to automatically plot the top 50 genes. Does not use the genelist in this case, but rather finds the top genes
               based on FDR and fold change.
        - category_label_columns: Name of columns that include the groups used for labeling individuals (e.g. group, metastasis, stage, etc.)
    
    Returns:
    -----------
        - Nothing
    '''
    
    if filetype == "csv":
        datadf = pd.read_csv(datafile, index_col = 0)
    elif filetype == "txt" or filetype == "tsv":
        datadf = pd.read_csv(datafile, index_col = 0, sep = "\t")
        
    os.makedirs(outdir, exist_ok=True)
    
    # Find the overlap between the samples in the data df and the metadata df
    samp_data_list = datadf.columns
    samp_meta_file = pd.read_csv(metadata)
    samp_meta_list = samp_meta_file.iloc[:,0].tolist()  

    sample_overlap = find_sample_overlap(samp_data_list, samp_meta_list)
    dat = datadf[datadf.columns.intersection(sample_overlap)] # remove non-overlap samples from data df
    samp_meta_file = samp_meta_file[samp_meta_file.iloc[:,0].isin(sample_overlap)] # remove non-overlap samples from metadata df
    
    # If the user has supplied the gene list, plot just the genes they supplied. Otherwise, plot the top genes based on FDR and fold change
    if not top:
        # Filter the data frame containing the degree/expression values for just the genes in the supplied gene list
        genes_to_plot = file_to_list(genelist)
        filtered = dat.filter(items = genes_to_plot, axis=0)

    else:
        filtered = filter_for_top_genes(datafile=dat, 
                       statsfile=compare_df,                     
                       number=50)
   
    #### Note: This code will manually calculate the z-scores, as it is copied from the heatmap code. However, the built-in z-score ####
    ####       option in the clustermap function does this as well. I have kept this manual calculation, however, as it allows for  ####
    ####       the export of the z-score file that a user could use for plotting in other software.                                 ####                 
    
    from scipy.cluster.hierarchy import linkage, dendrogram
    from scipy.stats import zscore
    filtered_z = filtered.apply(zscore, axis=1)
    Z = linkage(filtered_z, method='ward')
    leaf_order = dendrogram(Z, no_plot=True)['ivl']
    ordered_df = filtered_z.iloc[map(int, leaf_order),:]
    
    out_filtered_z_path = os.path.join(outdir, f"{prefix}_filtered_data_file_for_heatmap_genes_zscore.csv")
    ordered_df.to_csv(out_filtered_z_path)    
        
    min_val = np.percentile(filtered_z, 5) # vals will be the 1st and 99th percentile, so extreme (outlier) values will not influence the scale
    max_val = np.percentile(filtered_z, 95)     
        
    def _create_column_colors(data_df, metadata_df, cat_label_columns, column_colors):
        '''
        Description
            This function takes the metadata and zscore data frames and assigns color values to the
            samples based on the sample's sub-category. 

        Parameters:
        -----------
            - data_df: pd.DataFrame, pandas df with sample names as columns, degrees of genes or TFs as rows
            - metadata_df: pd.DataFrame, metadata pandas df, containing the following:
                1. A column called "name", which contains the name of the sample
                2. A column containing the group the sample belongs to 
                3. Additional columns that are used for classifying the samples into sub-categories (the 'metadata_column' parameter)

                Example:
                
                name group metastasis  stage
                TCGA_E2_A10E_01A_21R_A10J_07  LumA        yes      4
                TCGA_D8_A1JB_01A_11R_A13Q_07  LumA        yes      4
                TCGA_AR_A5QN_01A_12R_A28M_07  LumB        yes      3
                TCGA_LL_A50Y_01A_11R_A266_07  LumB        no      1
                TCGA_B6_A0I5_01A_11R_A034_07  LumB        no      2
            
            - cat_label_columns: list, Name of the columns in the metadata_df containing the addational sub-categories for classification 
              of samples
            - column_colors: list, List of dictionaries of the colors wished to use for each sub-category, as defined in the params.yml file
         
        Returns:
        -----------   
            A list, where the first item is the data frame that was originally input, but sorted by the sample category, so that
            samples of the same primary category (the first category listed by the user) will cluster together. The second objet of the 
            list is a dictionary of panda series, keyed on category label and values are panda series that are the order of samples, 
            mapped to their groups per category
        '''
        column_cat_level = 0 # Counter that keeps track of if the primary (first) category is being plotted, which is used to group same-type samples in the resulting plot

        assigned_colors = {} # keys are category titles, values are pandas series that are the order of samples, mapped to categories
        luts_dict = {} # For making the legend of the plot later

        # loop through column names ("Cancer_subtype", "Tumor_grade", etc.) and create a pandas series for each category to map the 
        # sample name to the column color in the resulting clustermap
        for category in cat_label_columns: 
            num_unique_subcategories = len(metadata_df[category].unique())
            num_unique_color_codes = len(column_colors[column_cat_level])
            if num_unique_subcategories != num_unique_color_codes:
                raise WrongAmountOfColorsError(category, num_unique_subcategories, num_unique_color_codes)
            # elif num_unique_subcategories < num_unique_color_codes:
            #     warnings.warn(f"\nWarning: For your '{category}' category, you have entered more colors ({num_unique_color_codes}) than sub-categories ({num_unique_subcategories}). Only the first {num_unique_subcategories} colors given will be used.\n")

            lut = column_colors[column_cat_level]
            luts_dict[category] = column_colors[column_cat_level]
            metadata_df['col_colors'] = metadata_df[category].map(lut)
            
            # column_cat_level is the counter. This code checks if we are looking at the first category the user gives in the params.yml 
            # file, and if so then it sorts the data frame so it "groups" the samples by that category. We can only group based off one
            # category, so we need this counter
            if column_cat_level == 0:
                data_df = data_df.T
                data_df['color'] = data_df.index.map(metadata_df.set_index('name').col_colors.squeeze())
                sorted_data_df = data_df.sort_values(by=['color'])
                assigned_colors[category] = sorted_data_df.pop("color")
            else:
                sorted_data_df['color'] = sorted_data_df.index.map(metadata_df.set_index('name').col_colors.squeeze())
                assigned_colors[category] = sorted_data_df.pop("color")
            
            column_cat_level += 1         
        
        return([sorted_data_df.T, assigned_colors, luts_dict])

    df_for_plotting, col_color_dict, luts = _create_column_colors(ordered_df, samp_meta_file, category_label_columns, category_column_colors)

    col_colors_df = pd.DataFrame.from_dict(col_color_dict)
    cbar_kws={"orientation": "vertical", "pad":0.02, "use_gridspec":True, "label":'z-score', "ticks":[-2,-1,0,1,2]}
    # kws = dict(cbar_kws=dict(label='z-score',, orientation='horizontal', center=0))
    # sns_plot = sns.clustermap(df_for_plotting, col_colors=col_colors_df, z_score=0, method='ward', row_cluster=True, col_cluster=False, 
    #                          dendrogram_ratio=0.05, vmin = -2, vmax = 2, cmap = matplotlib.colormaps['RdBu_r'], cbar_kws=cbar_kws)
    cbar_pos = (1, 0.5, 0.02, 0.25)

    print(row_cluster)
    print(df_for_plotting.head(10))
    df_for_plotting = df_for_plotting.reindex(genes_to_plot).reset_index()
    df_for_plotting = df_for_plotting.set_index('index')
    print(df_for_plotting.head(10))
    
    sns_plot = sns.clustermap(df_for_plotting, col_colors=col_colors_df, z_score=None, row_cluster=row_cluster, col_cluster=column_cluster, 
                              dendrogram_ratio=0.05, vmin = -2, vmax = 2, cmap = matplotlib.colormaps['RdBu_r'], cbar_kws=cbar_kws, cbar_pos=cbar_pos,
                              xticklabels=plot_sample_names, yticklabels=plot_gene_names)

    # Add legend to plot
    from matplotlib.patches import Patch
    
    #### I have absolutely no idea why the following code doesn't visualize properly. 
    # legend_counter = 0
    # for category in category_label_columns: # group, Grade, etc.
    #     xx = []
    #     for label in samp_meta_file[category].unique(): #group1, group2, etc
    #         x = sns_plot.ax_col_dendrogram.bar(0, 0, color=luts[category][label], label=label, linewidth=0)
    #         xx.append(x)
    #     if legend_counter == 0:
    #         legend = plt.legend(xx, samp_meta_file[category].unique(), loc="upper center", ncol=2, title=category, bbox_to_anchor=(.25, 1.05), bbox_transform=plt.gcf().transFigure)
    #     elif legend_counter == 1:
    #         legend = plt.legend(xx, samp_meta_file[category].unique(), loc="upper center", ncol=2, title=category, bbox_to_anchor=(.5, 1.05), bbox_transform=plt.gcf().transFigure)
    #     elif legend_counter == 2:
    #         legend = plt.legend(xx, samp_meta_file[category].unique(), loc="upper center", ncol=2, title=category, bbox_to_anchor=(.75, 1.05), bbox_transform=plt.gcf().transFigure)
    #     plt.gca().add_artist(legend)
    #     legend_counter += 1
    
    # Calulate x axis position of legends based on how many categories are given by the user
    x_positions = [x/(len(category_label_columns)+1) for x in range(len(category_label_columns)+2) if x != 0 and x != len(category_label_columns)+1]

    x_positions = [0.2, 0.4, 0.6, 0.85]

    for label in samp_meta_file[category_label_columns[0]].unique(): #group1, group2, etc
        sns_plot.ax_col_dendrogram.bar(0, 0, color=luts[category_label_columns[0]][label], label=label, linewidth=0)
        legend = sns_plot.ax_col_dendrogram.legend(title=category_label_columns[0], loc="upper center", ncol=2, bbox_to_anchor=(x_positions[0], 1.05), bbox_transform=plt.gcf().transFigure)
    
    xpos_counter = 1
    # legend_counter = 1
    for category in category_label_columns[1:len(category_label_columns)]: # group, Grade, etc.
        artists = []
        for label in samp_meta_file[category].unique(): #group1, group2, etc
            x = sns_plot.ax_col_dendrogram.bar(0, 0, color=luts[category][label], label=label, linewidth=0)
            artists.append(x)
        legend = plt.legend(artists, samp_meta_file[category].unique(), loc="upper center", ncol=2, title=category, bbox_to_anchor=(x_positions[xpos_counter], 1.05), bbox_transform=plt.gcf().transFigure)

        plt.gca().add_artist(legend)
        
        xpos_counter += 1
            
    outname = os.path.join(outdir, f"{prefix}_clustermap.png")
    sns_plot.savefig(outname, dpi = 600, bbox_inches='tight')
    
    print(f"\nFile created: {outname}")
    print(f"\nFile created: {out_filtered_z_path}")

    # Save the order of rows if desired, so that the same order can be applied on another dataset
    if row_cluster:
        row_cluster_order = [df_for_plotting.index[x] for x in sns_plot.dendrogram_row.reordered_ind]
        row_cluster_order_outname = os.path.join(outdir, f"{prefix}_clustered_row_order.txt")

        with open(row_cluster_order_outname, "w") as outfile:
            outfile.write("\n".join(row_cluster_order))

        print(f"\nFile created: {row_cluster_order_outname}")

