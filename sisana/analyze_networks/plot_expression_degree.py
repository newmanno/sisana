import seaborn as sns 
import pandas as pd
import argparse
import sys
import pickle
from pathlib import Path
import matplotlib.pyplot as plt
import warnings
from .analyze import file_to_list, NotASubsetError, filter_for_top_genes, filter_for_user_defined_genes, IncorrectHeaderError
import os

def plot_expression_degree(datafile: str, filetype: str, statsfile: str, metadata: str, 
                  plottype: str, groups: list, colors: list, prefix: str, yaxisname: str, outdir: str,
                  top: bool=True, numgenes: int=10, genelist: str=""):
    """
    Description:
        This code creates either a violin plot or a box plot of gene expression or degree data. The plots are annotated using the previously
        calculated statistics from the "sisana compare means" step
        
    Parameters:
    -----------
        - datafile: str, Path to file containing the expression or indegrees of each gene per sample.
        - filetype: str, Type of input file, must be either "csv", "txt", or "tsv", where csv implies comma separated values and txt/tsv 
                        implies tab-separated
        - statsfile: str, Path to the file that contains the comparison output of the gene expression or degree
        - metadata: str, Path to the csv metadata file mapping samples to groups (groups must match names of the groups arg), 
                        must have a header of the format 'name,group'
        - genelist: str, .txt file containing a list of genes to plot, must match the name of genes in the datafile. Recommended 
                        to not use more than 10 genes, otherwise use a heatmap.
        - plottype: str, The type of plot to create. Choices are "boxplot" or "violin"
        - groups: list(str), The names of the groups to plot, should ideally be a list of 5 names or less. Groups will be plotted 
                        in the order they are written in this argument
        - colors: str, The colors for each group, in the same order the groups appear in the groups arg
        - prefix: str, Prefix to use for the output figures
        - yaxisname: str, Name to use for the y-axis
        - outdir: str, Path to directory to output file to
        - top: Flag for whether to automatically plot the top 10 values. Does not use the genelist in this case, but rather finds the top genes
                        based on FDR and fold change.
        
    Returns:
    -----------
        - Nothing
    """
    
    # parser = argparse.ArgumentParser(description="Example command: python plot_expression_degree.py -d <indegree.csv> -f csv -m <metadata.csv> -g <genelist.txt> -s <samporder.txt> -p violin -n group1 group2 group3 -o ./output/")
    # ArgGroup = parser.add_argument_group('Required arguments')  
    
    # ArgGroup.add_argument("-d", "--datafile", type=str, help="Path to file containing the expression or indegrees of each gene per sample. If file has a header, denote this with the -f argument. Otherwise, supply the header using the -s argument", required=True)
    # ArgGroup.add_argument("-t", "--filetype", choices = ["csv", "txt"], help="Type of delimiter used for --datafile", required=True)    
    # # ArgGroup.add_argument("-f", "--fileheader", action='store_true', help="Flag for if --datafile has a header already. If not, requires --sampleorder", required=False)    
    # ArgGroup.add_argument("-m", "--metadata", type=str, help="Path to the csv metadata file mapping samples to groups (groups must match names of the --groupnames arg), must have a header of the format 'name,group'", required=True) 
    # ArgGroup.add_argument("-g", "--genelist", type=str, help=".txt file containing a list of genes to plot", required=True)   
    # # ArgGroup.add_argument("-s", "--sampleorder", type=str, help=".txt file containing the order of the samples in the data file", required=False)   
    # ArgGroup.add_argument("-p", "--plottype", type=str, choices = ["boxplot","violin"], help="The type of plot to create", required=True)   
    # ArgGroup.add_argument("-n", "--groupnames", type=str, nargs = "+", help="The names of the groups to plot, should ideally be a list of 5 names or less. Groups will be plotted in the order they are written in this argument", required=True)   
    # ArgGroup.add_argument("-c", "--colors", type=str, nargs = "+", help="The colors for each group, in the same order the groups appear in --groupnames", required=False) 
    # ArgGroup.add_argument("-x", "--prefix", type=str, help="Prefix to use for the output figures", required=True)         
    # ArgGroup.add_argument("-y", "--yaxis_name", type=str, help="Name to use for the y-axis", required=True)         
    # ArgGroup.add_argument("-o", "--outdir", type=str, help="Path to directory to output file to", required=True) 
    
    # args = parser.parse_args()
    
    # Check that there are no more than 5 group names, otherwise warn user
    if len(groups) > 5:
        warnings.warn("Warning: Supplying more than 10 groups at once may cause the created graphs to be unreadable. Consider reducing the number of groups.")
    
    # Get data and metadata
    if filetype == "csv":
        # indata = pd.read_csv(datafile, engine = "pyarrow", index_col=[0])
        indata = pd.read_csv(datafile, engine = "python", index_col=[0])
    elif filetype == "txt" or filetype == "tsv":
        # indata = pd.read_csv(datafile, sep='\t', engine = "pyarrow", header=None, index_col=[0])
        indata = pd.read_csv(datafile, sep='\t', engine = "python", index_col=[0])
        
    # meta = pd.read_csv(metadata, engine = "pyarrow", index_col=[0])
    meta = pd.read_csv(metadata, engine = "python", header=None, index_col=0)
    # if meta.columns[0] != "name":
    #     raise IncorrectHeaderError(meta)
    # else:
    #     meta = meta.set_index('name')
    meta.index.name = "name"
    meta.columns = ["group"]

    # Create list of genes if the user has supplied a gene list 
    if genelist != None:
        user_gene_list = file_to_list(genelist)

        # Check if the user-supplied gene list is a subset of the genes in the data
        data_genes = list(indata.index)
        
        if not set(user_gene_list).issubset(set(data_genes)):
            raise NotASubsetError(user_gene_list, data_genes, "genes")

    indata = indata.T
    indata.index.name = "name"
    indata = indata.T

    compare_df = pd.read_csv(statsfile, sep = "\t", index_col=0)
    
    # If the user has supplied the gene list, plot just the genes they supplied. Otherwise, plot the top genes based on FDR and fold change
    if genelist != None:
        # Filter the data frame containing the degree/expression values for just the genes in the supplied gene list
        filtered_indata_genelist = filter_for_user_defined_genes(datafile=indata, genes=user_gene_list)
        # filtered_indata_genelist = indata.loc[:,user_gene_list]
        filtered_indata_genelist = filtered_indata_genelist # Need to transform here, do not remove without testing

    else:
        filtered_indata_genelist = filter_for_top_genes(datafile=indata, 
                statsfile=compare_df,                     
                number=numgenes)
        topgenes = filtered_indata_genelist.index

        # # Find top genes
        # sorted_indata = compare_df.sort_values(['abs(difference_of_means)', 'FDR'], ascending=[False, True])
        # topgenes = sorted_indata.index[0:10]
        
        # try:
        #     # filter the compfile for just the topgenes
        #     filtered_indata_genelist = indata.loc[:,topgenes]
        # except KeyError:
        #     topgenes = list(set(indata.columns) & set(topgenes))
        #     filtered_indata_genelist = indata.loc[:,topgenes]

        
        # # Filter the data frame containing the degree/expression values for just the top genes 

    filtered_indata_genelist = filtered_indata_genelist.T

    # Remove samples that are not part of the user supplied groups to be plotted
    filtered_indata_genelist['dupename'] = list(filtered_indata_genelist.index)

    meta['dupename'] = list(meta.index)  

    filtered_indata_genelist['group'] = filtered_indata_genelist['dupename'].map(meta.set_index('dupename')['group'])

    subdata = filtered_indata_genelist[filtered_indata_genelist['group'].notnull()] 
    subdata = subdata.drop(['dupename','group'], axis=1) # Remove columns for the melt    

    # Melt to get into long format, which is required for plotting
    subdata_melt = subdata.melt(ignore_index=False).reset_index()
    subdata_melt.columns.values[1] = "gene"

    # Add group names to the melted data frame for plotting
    subdata_melt['dupename'] = list(subdata_melt['name'])
    subdata_melt['group'] = subdata_melt['dupename'].map(meta.set_index('dupename')['group'])
    subdata_melt = subdata_melt.drop('dupename', axis=1)
    
    # Sort group names based on the order that they appear in groups arg so the user can control the order the groups are plotted in
    subdata_melt['group'] = pd.Categorical(subdata_melt.group, ordered=True, categories=groups)
    subdata_melt = subdata_melt.sort_values('group')   
   
    # Group samples by gene in the melted data frame, required for plotting
    if genelist != None:
        subdata_melt['gene'] = pd.Categorical(subdata_melt.gene, ordered=True, categories=user_gene_list)
    else:
        subdata_melt['gene'] = pd.Categorical(subdata_melt.gene, ordered=True, categories=topgenes)
        
    subdata_melt = subdata_melt.sort_values(['group','gene'])
    subdata_melt.set_index('name')
    
    # Set colors for plotting if the user has specified colors      
    if colors is not None:
        assert len(colors) == len(groups)
        custom_colors  = colors
        sns.set_palette(custom_colors)
    
    # Make output directory if it does not already exist
    Path(outdir).mkdir(parents=True, exist_ok=True)   
    
    # Create plots
    plt.figure(figsize=(6, 6), dpi = 600) 
    plt.xticks(rotation=45)

    # Add statistical annotations to plot    
    hue_plot_params = {
    'data': subdata_melt,
    'x': 'gene',
    'y': 'value',
    "hue": "group",
    }
            
    # Create the nested list (i.e. "[[('gene1', 'group1'), ('gene1', 'group2')], [('gene2', 'group1'), ('gene2', 'group2')], etc.] 
    # for telling the Annotator how to structure the plots. Note that the Annotator does not do the comparison itself, however, since 
    # the FDRs have already been calculated previously, using all genes for the FDR calculation
    def _create_comparison_list(genelist: list):
        pairs = []
        
        for i in genelist:
            num = 0
            temp_list = []
            for j in subdata_melt['group'].unique():
                comparison = (i,j)
                temp_list.append(comparison)
                
                if num == 0:
                    num +=1
                else:
                    pairs.append(temp_list)
                    temp_list = []
        return(pairs)

    # Subset the data frame for just the samples belonging to the user-defined groups
    subdata_melt = subdata_melt[subdata_melt['group'].isin(groups)]

    if genelist != None:
        stat_annotation_pairs = _create_comparison_list(user_gene_list)
    else:
        stat_annotation_pairs = _create_comparison_list(topgenes)

    # print(stat_annotation_pairs)
    
    from statannotations.Annotator import Annotator
    
    if plottype == "violin":
        ax = sns.violinplot(**hue_plot_params, inner = None)
        plt.ylabel(yaxisname)
        # plt.tight_layout()
        outname = os.path.join(outdir, f"{prefix}_violin_plot.png")
        annotator = Annotator(ax, stat_annotation_pairs, **hue_plot_params, plot="violinplot", hide_non_significant = True)

    elif plottype == "boxplot":
        ax = sns.boxplot(**hue_plot_params, fliersize = 2)
        plt.ylabel(yaxisname)
        # plt.tight_layout()
        outname = os.path.join(outdir, f"{prefix}_box_plot.png")
        annotator = Annotator(ax, stat_annotation_pairs, **hue_plot_params, hide_non_significant = True)
    
    ax.set(xlabel=None)
        
    # Annotate with the FDRs that were calculated previously
    if genelist != None:
        pval_list = [compare_df.loc[gene, "FDR"] for gene in user_gene_list]
    else:
        pval_list = [compare_df.loc[gene, "FDR"] for gene in topgenes]
        
    annotator.configure(hide_non_significant=True)
    annotator.set_pvalues_and_annotate(pvalues = pval_list)    
    
    plt.figtext(0.5, 0.01, '* p < 0.05, ** p < 0.01, *** p < 0.001, **** p < 0.0001', horizontalalignment='center')
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.18) 
    plt.savefig(outname)
    print(f"\nFile saved: {outname}\n")