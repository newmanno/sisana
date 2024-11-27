import pandas as pd
import argparse
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text
import os
from .analyze import file_to_list


def plot_volcano(datafile: str, fccol: str, adjpcol: str, fcthreshold: str, adjpvalthreshold: str, numlabels: str, outdir: str):
    """
    Description:
        This code performs a survival analysis between two user-defined groups and outputs
        both the survival plot and the statistics for the comparison(s)
        
    Parameters:
    -----------
        - datafile: str, Path to tab delimited file containing the fold change, p-value, FDR, and mean 
                         degree/expression for each gene. This is reported with the compare_groups.py script
        - fccol: str, The name of the column containing the log2 fold change values
        - adjpcol: str, The name of the column containing the adj. p-value
        - fcthreshold: str, Fold change threshold to use (value will also be applied to the negative end of the 
                            fold change axis, meaning a value of 1.5 really means +/- 1.5)
        - adjpvalthreshold: str, Threshold to use for the adjusted p-value
        - numlabels: str, Number of top values to label
        - outdir: str, Path to directory to output file to
        
    Returns:
    -----------
        - Nothing
    """
    
    # parser = argparse.ArgumentParser(description="Example command: python volcano_plot.py -d indegree.csv -f csv -m metadata.csv -g genelist.txt -s samporder.txt -p violin -n group1 group2 -o ./output/")
    # ArgGroup = parser.add_argument_group('Required arguments')  
    
    # ArgGroup.add_argument("-d", "--datafile", type=str, help="Path to file containing the fold change, p-value, FDR, and mean degree/expression for each gene. This is reported with the compare_groups.py script.", required=True)
    # ArgGroup.add_argument("-f", "--fcthreshold", type=float, help="Fold change threshold to use (value will also be applied to the negative end of the fold change axis, meaning a value of 1.5 really means +/- 1.5)", required=True)
    # ArgGroup.add_argument("-p", "--adjpvalthreshold", type=float, help="Threshold to use for the adjusted p-value", required=True)    
    # ArgGroup.add_argument("-t", "--topvalstoplot", type=int, help="Number of top values to label", required=True)   
    # ArgGroup.add_argument("-o", "--outdir", type=str, help="Path to directory to output file to", required=True) 
    
    # Create output directory if one does not already exist
    os.makedirs(outdir, exist_ok=True)
    
    # args = parser.parse_args()
    exp = pd.read_csv(datafile, index_col = 0, sep = "\t")
    print(exp)

    # plot 
    plt.figure(figsize=(6, 6), dpi = 1200) 
    plt.scatter(x=exp[fccol],y=exp[adjpcol].apply(lambda x:-np.log10(x)), s=1)
    
    down = exp[(exp[fccol] <= -1 * fcthreshold) & (exp[adjpcol] <= adjpvalthreshold)]
    up = exp[(exp[fccol] >= fcthreshold) & (exp[adjpcol] <= adjpvalthreshold)]
    not_down_or_up = exp[(exp[fccol] < fcthreshold) & (exp[fccol] > -1 * fcthreshold) & (exp[adjpcol] < adjpvalthreshold)]
    notsig = exp[exp[adjpcol] > adjpvalthreshold]

    plt.scatter(x=down[fccol], y=down[adjpcol].apply(lambda x:-np.log10(x)), s=3, label="Down-regulated",color="blue")
    plt.scatter(x=up[fccol], y=up[adjpcol].apply(lambda x:-np.log10(x)), s=3, label="Up-regulated",color="orange")
    plt.scatter(x=notsig[fccol], y=notsig[adjpcol].apply(lambda x:-np.log10(x)), s=3, label="Not significant",color="gainsboro")
    plt.scatter(x=not_down_or_up[fccol], y=not_down_or_up[adjpcol].apply(lambda x:-np.log10(x)), s=3, color="darkgrey")

    uptexts=[]
    if len(up) >= numlabels:
        up_subset = up.sort_values(adjpcol, ascending=True).head(numlabels)
        
        for i,r in up_subset.iterrows():
            uptexts.append(plt.text(x=r[fccol],y=-np.log10(r[adjpcol]),s=i))
            
    else:
        up_subset = up.sort_values(adjpcol, ascending=True).head(len(up))

        for i,r in up.iterrows():
            uptexts.append(plt.text(x=r[fccol],y=-np.log10(r[adjpcol]),s=i))
        
    adjust_text(uptexts,arrowprops=dict(arrowstyle="-", color='black', lw=0.5))

    dntexts=[]
    if len(down) >= numlabels:
        dn_subset = down.sort_values(adjpcol, ascending=True).head(numlabels)
        
        for i,r in dn_subset.iterrows():
            dntexts.append(plt.text(x=r[fccol],y=-np.log10(r[adjpcol]),s=i))
            
    else:
        dn_subset = down.sort_values(adjpcol, ascending=True).head(len(up))

        for i,r in dn_subset.iterrows():
            dntexts.append(plt.text(x=r[fccol],y=-np.log10(r[adjpcol]),s=i))
        
    adjust_text(dntexts,arrowprops=dict(arrowstyle="-", color='black', lw=0.5))

    plt.xlabel("log(FC)")
    plt.ylabel("-log10(FDR)")
    plt.axvline(-1 * fcthreshold, color="grey", linestyle="--")
    plt.axvline(fcthreshold, color="grey", linestyle="--")
    plt.axhline(-np.log10(adjpvalthreshold), color="grey", linestyle="--")
    #plt.legend()
    
    outname = os.path.join(outdir, f"volcano_plot_FC_{fcthreshold}_adjp_{adjpvalthreshold}_top_{numlabels}.png")
    plt.savefig(outname)
        
    print(f"File saved: {outname}")     
    
