import argparse
import os
import pandas as pd
import pickle
import numpy as np
from .analyze import file_to_list, map_samples
import gseapy as gp
from matplotlib import pyplot as plt

__author__ = 'Nolan Newman'
__contact__ = 'nolankn@uio.no'
    
def perform_gsea(genefile: str, gmtfile: str, geneset: str, outdir: str):
    """
    Description:
        This code performs a survival analysis between two user-defined groups and outputs
        both the survival plot and the statistics for the comparison(s)
        
    Parameters:
    -----------
        - genefile: str, Path to file (.rnk format, which is two column, tab delimited, no header) 
                         containing the genes and test statistics to do enrichment on
        - gmtfile: str, Path to the gene set file in gmt format
        - geneset: str, The gene set type used for gmtfile
        - outdir: str, Path to directory to output file to
        
    Returns:
    -----------
        - Nothing
    """

    # parser = argparse.ArgumentParser(description="Example command: python gsea.py -f file.rnk -g genesets.gmt -s Hallmarks -o ./output")
    # ArgGroup = parser.add_argument_group('Required arguments') 
    # ArgGroup.add_argument("-f", "--genefile", type=str, help="Path to file (.rnk format, which is two column, tab delimited, no header) containing the genes and test statistics to do enrichment on", required=True) 
    # ArgGroup.add_argument("-g", "--gmtfile", type=str, help="Path to the gene set file in gmt format", required=True)     
    # ArgGroup.add_argument("-s", "--geneset", type=str, help="The gene set type you input into --gmtfile", required=True)            
    # ArgGroup.add_argument("-o", "--outdir", type=str, help="Path to directory to output file to", required=True) 
    
    # args = parser.parse_args()
    
    # Create output directory if one does not already exist
    os.makedirs(outdir, exist_ok=True)
    
    # Make user specified directory if it does not already exist
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    rnk = pd.read_csv(genefile, header=None, index_col=0, sep="\t")
    print("rnk format")
    print(rnk)
    
    # Run GSEA on a pre-ranked list of genes
    pre_res = gp.prerank(rnk=genefile,
                     gene_sets=gmtfile,
                     threads=4,
                     min_size=5,
                     max_size=1000,
                     permutation_num=1000, # reduce number to speed up testing
                     outdir=outdir,
                     seed=6,
                     verbose=True, # see what's going on behind the scenes
                    )
  
    # Perform enrichment       
    # print("Performing enrichment...")
    # enr = gp.enrichr(gene_list = genes, 
    #                 gene_sets = gset,
    #                 organism = 'human',
    #                 outdir = outdir,
    #                 )
    
    # print(pre_res.res2d)
    
    pre_res_df = pd.DataFrame(pre_res.res2d)
    pre_res_df = pre_res_df.sort_values('NOM p-val', ascending = True)
    print(pre_res_df)

    res_file_name = os.path.join(outdir, f"Prerank_GSEA_{geneset}_results.txt")
    pre_res_df.to_csv(res_file_name, sep = "\t", index = False)
    
    terms = pre_res.res2d.Term
    
    # Plot top 5 terms in basic GSEA plot
    plt.subplots_adjust(left=0.3, right=0.9, bottom=0.3, top=0.9)
    ax = pre_res.plot(terms=terms[1:6],
                   show_ranking=True,
                   figsize=(15,20)
                  )
    
    gsea_plot_name = os.path.join(outdir, f"GSEA_{geneset}_basic_enrichment_plot.png")
    ax.figure.savefig(gsea_plot_name, bbox_inches = "tight")
    
    # Plot significant GSEA terms
    from gseapy import dotplot
    ax = dotplot(pre_res.res2d,
                column="FDR q-val",
                title=geneset,
                cmap=plt.cm.viridis,
                size=10,
                figsize=(10,20), 
                cutoff=0.25, 
                show_ring=False)

    dotplot_name = os.path.join(outdir, f"GSEA_{geneset}_basic_enrichment_dotplot.png")
    ax.figure.savefig(dotplot_name, bbox_inches = "tight")

    print("\nDone!")
    
    print(f"Files created:\n{res_file_name}\n{gsea_plot_name}\n{dotplot_name}\n")
