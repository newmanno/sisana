import seaborn as sns 
import pandas as pd
import argparse
import sys
import pickle
from pathlib import Path
import matplotlib.pyplot as plt
import warnings

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
    ArgGroup.add_argument("-n", "--groupnames", type=str, nargs = "+", help="The names of the groups to plot, should ideally be a list of 5 names or less. Groups will be plotted in the order they are written in this argument", required=True)   
    ArgGroup.add_argument("-o", "--outdir", type=str, help="Path to directory to output file to", required=True) 
    
    args = parser.parse_args()
    
    lut = dict(zip(species.unique(), "rbg"))
    row_colors = species.map(lut)
    sns.clustermap(iris, row_colors=row_colors)