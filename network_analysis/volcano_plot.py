import pandas as pd
import argparse
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text
import os
from analyze import file_to_list


if __name__ == '__main__':
    """
    Description:
        This code creates a volcano plot from the p-value/FDR and log2 fold change
    """
    
    parser = argparse.ArgumentParser(description="Example command: python plot_expression_degree.py -d <indegree.csv> -f csv -m <metadata.csv> -g <genelist.txt> -s <samporder.txt> -p violin -n group1 group2 group3 -o ./output/")
    ArgGroup = parser.add_argument_group('Required arguments')  
    
    ArgGroup.add_argument("-d", "--datafile", type=str, help="Path to file containing the expression or indegrees of each gene per sample, has no header", required=True)
    ArgGroup.add_argument("-f", "--fcthreshold", type=float, help="Fold change threshold to use (value will also be applied to the negative end of the fold change axis)", required=True)
    ArgGroup.add_argument("-p", "--adjpvalthreshold", type=float, help="Threshold to use for the adjusted p-value", required=True)    
    ArgGroup.add_argument("-t", "--topvalstoplot", type=int, help="Number of top values to plot", required=True)   
    ArgGroup.add_argument("-o", "--outdir", type=str, help="Path to directory to output file to", required=True) 
    
    args = parser.parse_args()
    exp = pd.read_csv(args.datafile, index_col = 0)

    # plot 
    plt.figure(figsize=(6, 6), dpi = 1200) 

    plt.scatter(x=exp['logFC'],y=exp['adjpval'].apply(lambda x:-np.log10(x)), s=1)
    
    down = exp[(exp['logFC'] <= -1 * args.fcthreshold) & (exp['adjpval'] <= args.adjpvalthreshold)]
    up = exp[(exp['logFC'] >= args.fcthreshold) & (exp['adjpval'] <= args.adjpvalthreshold)]
    not_down_or_up = exp[(exp['logFC'] < args.fcthreshold) & (exp['logFC'] > -1 * args.fcthreshold) & (exp['adjpval'] < args.adjpvalthreshold)]
    notsig = exp[exp['adjpval'] > args.adjpvalthreshold]

    plt.scatter(x=down['logFC'], y=down['adjpval'].apply(lambda x:-np.log10(x)), s=3, label="Down-regulated",color="blue")
    plt.scatter(x=up['logFC'], y=up['adjpval'].apply(lambda x:-np.log10(x)), s=3, label="Up-regulated",color="orange")
    plt.scatter(x=notsig['logFC'], y=notsig['adjpval'].apply(lambda x:-np.log10(x)), s=3, label="Not significant",color="gainsboro")
    plt.scatter(x=not_down_or_up['logFC'], y=not_down_or_up['adjpval'].apply(lambda x:-np.log10(x)), s=3, color="darkgrey")

    uptexts=[]
    if len(up) >= args.topvalstoplot:
        up_subset = up.sort_values('adjpval', ascending=True).head(args.topvalstoplot)
        
        for i,r in up_subset.iterrows():
            uptexts.append(plt.text(x=r['logFC'],y=-np.log10(r['adjpval']),s=i))
            
    else:
        up_subset = up.sort_values('adjpval', ascending=True).head(len(up))

        for i,r in up.iterrows():
            uptexts.append(plt.text(x=r['logFC'],y=-np.log10(r['adjpval']),s=i))
        
    adjust_text(uptexts,arrowprops=dict(arrowstyle="-", color='black', lw=0.5))

    dntexts=[]
    if len(down) >= args.topvalstoplot:
        dn_subset = down.sort_values('adjpval', ascending=True).head(args.topvalstoplot)
        
        for i,r in dn_subset.iterrows():
            dntexts.append(plt.text(x=r['logFC'],y=-np.log10(r['adjpval']),s=i))
            
    else:
        dn_subset = down.sort_values('adjpval', ascending=True).head(len(up))

        for i,r in dn_subset.iterrows():
            dntexts.append(plt.text(x=r['logFC'],y=-np.log10(r['adjpval']),s=i))
        
    adjust_text(dntexts,arrowprops=dict(arrowstyle="-", color='black', lw=0.5))

    plt.xlabel("log(FC)")
    plt.ylabel("-log10(FDR)")
    plt.axvline(-1 * args.fcthreshold, color="grey", linestyle="--")
    plt.axvline(args.fcthreshold, color="grey", linestyle="--")
    plt.axhline(-np.log10(args.adjpvalthreshold), color="grey", linestyle="--")
    #plt.legend()
    
    outname = os.path.join(args.outdir, f"volcano_plot_FC_{args.fcthreshold}_adjp_{args.adjpvalthreshold}_top_{args.topvalstoplot}.pdf")
    plt.savefig(outname)
        
    print(f"File saved: {outname}")     
    
