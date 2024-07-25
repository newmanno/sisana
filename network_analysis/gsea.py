import argparse
import os
import pandas as pd
import pickle
import numpy as np
from analyze import file_to_list, map_samples
import gseapy as gp
from matplotlib import pyplot as plt

__author__ = 'Nolan Newman'
__contact__ = 'nolankn@uio.no'
    
if __name__ == '__main__':
    """
    Description:
        This code performs gene set enrichment analysis on a set of pre-ranked genes 
    """

    parser = argparse.ArgumentParser(description="Example command: python gsea.py -f file.rnk -g genesets.gmt -s Hallmarks -o ./output")
    ArgGroup = parser.add_argument_group('Required arguments') 
    ArgGroup.add_argument("-f", "--genefile", type=str, help="Path to file (.rnk format, which is two column, tab delimited, no header) containing the genes and test statistics to do enrichment on", required=True) 
    ArgGroup.add_argument("-g", "--gmtdir", type=str, help="Path to the directory containing the gene set files in gmt format", required=True)     
    ArgGroup.add_argument("-s", "--geneset", type=str, choices = ["KEGG", "Hallmarks", "Reactome"], help="The gene set type to use", required=True)         
    ArgGroup.add_argument("-o", "--outdir", type=str, help="Path to directory to output file to", required=True) 
    
    args = parser.parse_args()
    
    # Make user specified directory if it does not already exist
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    
    # Get the correct gene set based on user input
    if args.geneset == "KEGG":
        gset = os.path.join(args.gmtdir, "c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt")
    elif args.geneset == "Hallmarks":
        gset = os.path.join(args.gmtdir, "Hallmark.v2023.2.Hs.symbols.gmt")
    elif args.geneset == "Reactome":
        gset = os.path.join(args.gmtdir, "c2.cp.reactome.v2023.2.Hs.symbols.gmt")
    
    # Run GSEA on a pre-ranked list of genes
    pre_res = gp.prerank(rnk=args.genefile,
                     gene_sets=gset,
                     threads=4,
                     min_size=5,
                     max_size=1000,
                     permutation_num=1000, # reduce number to speed up testing
                     outdir=args.outdir,
                     seed=6,
                     verbose=True, # see what's going on behind the scenes
                    )
  
    # Perform enrichment       
    # print("Performing enrichment...")
    # enr = gp.enrichr(gene_list = genes, 
    #                 gene_sets = gset,
    #                 organism = 'human',
    #                 outdir = args.outdir,
    #                 )
    
    # print(pre_res.res2d)
    
    pre_res_df = pd.DataFrame(pre_res.res2d)
    pre_res_df = pre_res_df.sort_values('NOM p-val', ascending = True)
    print(pre_res_df)

    res_file_name = os.path.join(args.outdir, f"Prerank_GSEA_{args.geneset}_results.txt")
    pre_res_df.to_csv(res_file_name, sep = "\t", index = False)
    
    terms = pre_res.res2d.Term
    
    # Plot top 5 terms in basic GSEA plot
    plt.subplots_adjust(left=0.3, right=0.9, bottom=0.3, top=0.9)
    ax = pre_res.plot(terms=terms[1:6],
                   show_ranking=True,
                   figsize=(15,20)
                  )
    
    gsea_plot_name = os.path.join(args.outdir, f"GSEA_{args.geneset}_plot.pdf")
    ax.figure.savefig(gsea_plot_name, bbox_inches = "tight")
    
    # Plot significant GSEA terms
    from gseapy import dotplot
    ax = dotplot(pre_res.res2d,
                column="FDR q-val",
                title=args.geneset,
                cmap=plt.cm.viridis,
                size=10,
                figsize=(10,20), 
                cutoff=0.25, 
                show_ring=False)

    dotplot_name = os.path.join(args.outdir, f"GSEA_dotplot_{args.geneset}.pdf")
    ax.figure.savefig(dotplot_name, bbox_inches = "tight")

    print("Done!")
    
    print("Files created:")
    print(f"Files created:\n{res_file_name}\n{gsea_plot_name}\n{dotplot_name}\n")
