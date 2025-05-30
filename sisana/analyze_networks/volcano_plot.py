import pandas as pd
import argparse
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text
import os
from .analyze import file_to_list, filter_for_top_genes, filter_for_user_defined_genes, create_label_list
import warnings

def plot_volcano(statsfile: str, diffcol: str, adjpcol: str, adjpvalthreshold: str, xaxisthreshold: float,  difftype: str,
                 outdir: str, top: bool=True, numlabels: int=15, genelist: str=""):
    """
    Description:
        This code performs a survival analysis between two user-defined groups and outputs
        both the survival plot and the statistics for the comparison(s)
        
    Parameters:
    -----------
        - statsfile: str, Path to tab delimited file containing the fold change, p-value, FDR, and mean 
          degree/statsression for each gene. This is reported with the compare_groups.py script
        - fccol: str, The name of the column containing the difference in medians or means
        - adjpcol: str, The name of the column containing the adj. p-value
        - adjpvalthreshold: str, Threshold to use for the adjusted p-value
        - genelist: str, Path to a .txt file containing a list of genes to plot. Alternatively, the top {numlabels} genes can be plotted instead if top=True.
        - outdir: str, Path to directory to output file to
        - difftype: str, The type of difference to use for the x-axis. "mean" will be difference in means and "median"
          refers to difference in medians
        - top: Flag for whether to automatically label the top 10 values. Does not use the genelist in this case, but rather finds the top genes
          based on FDR and fold change.
        - numlabels: int, Number of top values to label. Can only be used if top=True.
        
    Returns:
    -----------
        - Nothing 
    """
    
    # parser = argparse.ArgumentParser(description="Example command: python volcano_plot.py -d indegree.csv -f csv -m metadata.csv -g genelist.txt -s samporder.txt -p violin -n group1 group2 -o ./output/")
    # ArgGroup = parser.add_argument_group('Required arguments')  
    
    # ArgGroup.add_argument("-d", "--datafile", type=str, help="Path to file containing the fold change, p-value, FDR, and mean degree/statsression for each gene. This is reported with the compare_groups.py script.", required=True)
    # ArgGroup.add_argument("-f", "--fcthreshold", type=float, help="Fold change threshold to use (value will also be applied to the negative end of the fold change axis, meaning a value of 1.5 really means +/- 1.5)", required=True)
    # ArgGroup.add_argument("-p", "--adjpvalthreshold", type=float, help="Threshold to use for the adjusted p-value", required=True)    
    # ArgGroup.add_argument("-t", "--topvalstoplot", type=int, help="Number of top values to label", required=True)   
    # ArgGroup.add_argument("-o", "--outdir", type=str, help="Path to directory to output file to", required=True) 
    
    
    # Create output directory if one does not already exist
    os.makedirs(outdir, exist_ok=True)
    
    # args = parser.parse_args()
    stats = pd.read_csv(statsfile, index_col = 0, sep = "\t")    
    
    # Create initial plot 
    plt.figure(figsize=(6, 6), dpi = 1200) 
    plt.scatter(x=stats[diffcol],y=stats[adjpcol].apply(lambda x:-np.log10(x)), s=1)

    down = stats[(stats[diffcol] <= xaxisthreshold) & (stats[adjpcol] <= adjpvalthreshold)]
    up = stats[(stats[diffcol] > xaxisthreshold) & (stats[adjpcol] <= adjpvalthreshold)]
    not_down_or_up = stats[(stats[diffcol] < xaxisthreshold) & (stats[diffcol] > -1 * xaxisthreshold) & (stats[adjpcol] < adjpvalthreshold)]
    notsig = stats[stats[adjpcol] > adjpvalthreshold]
    
    if len(down) + len(up) == 0:
        raise Exception("Error: No significant values found to plot.")

    plt.scatter(x=down[diffcol], y=down[adjpcol].apply(lambda x:-np.log10(x)), s=3, label="Down-regulated",color="blue")
    plt.scatter(x=up[diffcol], y=up[adjpcol].apply(lambda x:-np.log10(x)), s=3, label="Up-regulated",color="orange")
    plt.scatter(x=notsig[diffcol], y=notsig[adjpcol].apply(lambda x:-np.log10(x)), s=3, label="Not significant",color="gainsboro")
    plt.scatter(x=not_down_or_up[diffcol], y=not_down_or_up[adjpcol].apply(lambda x:-np.log10(x)), s=3, color="darkgrey")

    # Label the user-defined genes if given, otherwise plot just the top genes
    if genelist != None:
        # Find the overlap in the user gene list and the up/down genes 
        genelist = file_to_list(genelist)
        
        user_defined_upgenes = list(set(genelist).intersection(up.index))
        user_defined_dngenes = list(set(genelist).intersection(down.index))
        
        if len(user_defined_upgenes) + len(user_defined_dngenes) != len(genelist):
            warnings.warn("\nWarning: Some genes in the input list are not significant and will not be labeled.")

        up_subset = filter_for_user_defined_genes(datafile=up, genes=user_defined_upgenes, verbose=False)
        dn_subset = filter_for_user_defined_genes(datafile=down, genes=user_defined_dngenes, verbose=False)
        
        uptexts = create_label_list(statsfile=up_subset, difference_column=diffcol, adjp_column=adjpcol)
        dntexts = create_label_list(statsfile=dn_subset, difference_column=diffcol, adjp_column=adjpcol)

    else: # Plot the top genes
        print("No list of genes to label was given, so labeling the top genes instead. Adjust the 'numlabels' parameter to change the number of genes labeled.")
        if len(up) >= numlabels:
            up_subset = up.sort_values(adjpcol, ascending=True).head(numlabels)
        else:
            warnings.warn("Warning: There are not enough significant 'up' values that pass the threshold. Only those that pass will be labeled.")
            up_subset = up.sort_values(adjpcol, ascending=True).head(len(up))

        uptexts = create_label_list(statsfile=up_subset, difference_column=diffcol, adjp_column=adjpcol)

        if len(down) >= numlabels:
            dn_subset = down.sort_values(adjpcol, ascending=True).head(numlabels)   
            print(dn_subset)
        else:
            warnings.warn("Warning: There are not enough significant 'down' values that pass the threshold. Only those that pass will be labeled.")
            dn_subset = down.sort_values(adjpcol, ascending=True).head(len(up))
            
        dntexts = create_label_list(statsfile=dn_subset, difference_column=diffcol, adjp_column=adjpcol)
        # print(dntexts)
        # sys.exit(0)


    # Add labels to plot
    adjust_text(uptexts,arrowprops=dict(arrowstyle="-", color='black', lw=0.5))        
    adjust_text(dntexts,arrowprops=dict(arrowstyle="-", color='black', lw=0.5))
    
    difftype= "median"
    if difftype == "mean":
        plt.xlabel("Difference in mean degree")
    elif difftype == "median":
        plt.xlabel("Difference in median degree")

    plt.ylabel("-log10(FDR)")

    fcthreshold = 50
    plt.axvline(-1 * fcthreshold, color="grey", linestyle="--")
    plt.axvline(fcthreshold, color="grey", linestyle="--")

    plt.axhline(-np.log10(adjpvalthreshold), color="grey", linestyle="--")
    plt.legend()
    
    if not top:
        outname = os.path.join(outdir, f"volcano_plot_adjp_{adjpvalthreshold}.png")
    else:
        outname = os.path.join(outdir, f"volcano_plot_adjp_{adjpvalthreshold}_top_{numlabels}.png")
    plt.savefig(outname)
        
    print(f"File saved: {outname}")     
    
