import seaborn as sns 
import pandas as pd
import argparse
import sys
import pickle
from pathlib import Path
import matplotlib.pyplot as plt
import warnings
from analyze import file_to_list
import os

if __name__ == '__main__':
    """
    Description:
        This code plots the expression or degrees of genes used in the Lioness pipeline
    """
    
    parser = argparse.ArgumentParser(description="Example command: python plot_expression_degree.py -d <indegree.csv> -f csv -m <metadata.csv> -g <genelist.txt> -s <samporder.txt> -p violin -n group1 group2 group3 -o ./output/")
    ArgGroup = parser.add_argument_group('Required arguments')  
    
    ArgGroup.add_argument("-d", "--datafile", type=str, help="Path to file containing the expression or indegrees of each gene per sample, has no header", required=True)
    ArgGroup.add_argument("-t", "--filetype", choices = ["csv", "txt"], help="Type of delimiter used for --datafile", required=True)    
    ArgGroup.add_argument("-f", "--fileheader", action='store_true', help="Flag for if --datafile has a header already. If not, requires --sampleorder", required=False)    
    ArgGroup.add_argument("-m", "--metadata", type=str, help="Path to the csv metadata file mapping samples to groups (groups must match names of the --groupnames arg), must have a header of the format 'name,group'", required=True) 
    ArgGroup.add_argument("-g", "--genelist", type=str, help=".txt file containing a list of genes to plot", required=True)   
    ArgGroup.add_argument("-s", "--sampleorder", type=str, help=".txt file containing the order of the samples in the data file", required=False)   
    ArgGroup.add_argument("-p", "--plottype", type=str, choices = ["boxplot","violin"], help="The type of plot to create", required=True)   
    ArgGroup.add_argument("-n", "--groupnames", type=str, nargs = "+", help="The names of the groups to plot, should ideally be a list of 5 names or less. Groups will be plotted in the order they are written in this argument", required=True)   
    ArgGroup.add_argument("-c", "--colors", type=str, nargs = "+", help="The colors for each group, in the same order the groups appear in --groupnames", required=False) 
    ArgGroup.add_argument("-x", "--prefix", type=str, help="Prefix to use for the output figures", required=True)         
    ArgGroup.add_argument("-o", "--outdir", type=str, help="Path to directory to output file to", required=True) 
    
    args = parser.parse_args()
    
    # Check that there are no more than 5 group names, otherwise warn user
    if len(args.groupnames) > 5:
        warnings.warn("Warning: Supplying more than 5 groups at once may cause the created graphs to be unreadable. Consider reducing the number of groups.")
    
    # Get data and metadata
    if args.filetype == "csv":
        # indata = pd.read_csv(args.datafile, engine = "pyarrow", index_col=[0])
        indata = pd.read_csv(args.datafile, engine = "python", index_col=[0])
    elif args.filetype == "txt":
        # indata = pd.read_csv(args.datafile, sep='\t', engine = "pyarrow", header=None, index_col=[0])
        indata = pd.read_csv(args.datafile, sep='\t', engine = "python", header=None, index_col=[0])
        
    # meta = pd.read_csv(args.metadata, engine = "pyarrow", index_col=[0])
    meta = pd.read_csv(args.metadata, engine = "python", index_col=[0])

    # Create list of genes
    user_gene_list = file_to_list(args.genelist)
  
    # If the input file does not have a header, then assign the data table a header based on the list of sample names supplied by the user
    if not args.fileheader:
        samps = file_to_list(args.sampleorder)
        indata.columns = samps  

    print(f"There are {len(user_gene_list)} genes to be plotted.")
    print(f"There are {len(indata.columns)} samples in the input file. Only those belonging to the supplied group names in --groupnames will be plotted.")
    
    # Check if the user-supplied gene list is a subset of the genes in the data
    data_genes = list(indata.index)
    
    if not set(user_gene_list).issubset(set(data_genes)):
        genes_missing = list(set(user_gene_list).difference(data_genes)) # Find genes not in the data provided
        print(genes_missing)
        print("\nError: The following genes in the supplied gene list are not present in the data provided:")
        [print(i) for i in genes_missing]
        print("\n")

        raise Exception("Please ensure the genes in the supplied list are a subset of the genes in the dataset.")

    indata = indata.T
    indata.index.name = "name"

    # Filter the data frame containing the degree/expression values for just the genes in the supplied gene list
    filtered_indata_genelist = indata.loc[:,user_gene_list]

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
    subdata_melt['dupename'] = list(subdata_melt.name)
    subdata_melt['group'] = subdata_melt['dupename'].map(meta.set_index('dupename')['group'])
    subdata_melt = subdata_melt.drop('dupename', axis=1)

    # Sort group names based on the order that they appear in --groupnames so the user can control the order the groups are plotted in
    subdata_melt['group'] = pd.Categorical(subdata_melt.group, ordered=True, categories=args.groupnames)
    subdata_melt = subdata_melt.sort_values('group')
    
    subdata_melt['gene'] = pd.Categorical(subdata_melt.gene, ordered=True, categories=user_gene_list)
    subdata_melt = subdata_melt.sort_values(['group','gene'])
    
    print(subdata_melt)

    # Set colors for plotting if the user has specified colors      
    if args.colors is not None:
        assert len(args.colors) == len(args.groupnames)
        custom_colors  = args.colors
        sns.set_palette(custom_colors)
    
    # Create plots
    plt.figure(figsize=(6, 6), dpi = 600) 
    plt.xticks(rotation=45)
    
    if args.plottype == "violin":
        sns.violinplot(data = subdata_melt, x = "gene", y = "value", hue = "group")
        outname = os.path.join(args.outdir, f"{args.prefix}_violin_plot.png")
        plt.savefig(outname)
        print(f"File saved: {outname}")
    elif args.plottype == "boxplot":
        sns.boxplot(data = subdata_melt, x = "gene", y = "value", hue = "group", fliersize = 2)
        outname = os.path.join(args.outdir, f"{args.prefix}_box_plot.png")
        plt.savefig(outname)
        print(f"File saved: {outname}")
