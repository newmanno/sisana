# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import argparse
import sys
import os
from pathlib import Path
import yaml

def preprocess_data(exp: str, motif: str, ppi: str, number: int, outdir: str):

    # parser = argparse.ArgumentParser(description="Example command: python preprocess.py -e <expression_file.tsv> -m <motif_file.txt> -p <ppi_file.txt> -n 10 -o ./output/")
    # requiredArgGroup = parser.add_argument_group('Required arguments')        
    # requiredArgGroup.add_argument("-e", "--exp", type=str, help="Path to file containing the gene expression data. Row names must be genes, the expression file does not have a header, and the cells are separated by a tab", required=True)
    # requiredArgGroup.add_argument("-m", "--motif", type=str, help="Path to motif file, which gets filtered to only contain genes that pass the minimum number of samples threshold", required=True)
    # requiredArgGroup.add_argument("-p", "--ppi", type=str, help="Path to ppi file, which gets filtered to only contain genes that pass the minimum number of samples threshold", required=True)   
    # requiredArgGroup.add_argument("-n", "--number", type=int, help="Minimum number of samples a gene must be expressed in; expression data will be filtered for only genes that pass", required=True)
    # requiredArgGroup.add_argument("-o", "--outdir", type=str, help="Path to output directory", required=True)    
    
    # args = parser.parse_args()
    
    # exp_file = args.exp
    # motif_file = args.motif
    # ppi_file = args.ppi
    
    # Create output file prefixes by removing the .txt suffixes
    expoutfile = Path(exp).stem
    motoutfile = Path(motif).stem
    ppioutfile = Path(ppi).stem
    
    # Create output directory if one does not already exist
    os.makedirs(outdir, exist_ok=True)
    
    print("Now filtering, please wait as this can take some time depending on the size of the input dataset...\n")
    
    expdf = pd.read_csv(exp, sep='\t', header=None, index_col=0)
    num_originalgenes = len(expdf)
    nsamps = len(expdf.columns)
    header_samps = ["sample_"+str(i) for i in list(range(1,nsamps + 1))] # Create the headers for samples
    
    expdf["stdev"] = expdf.std(axis=1) # Calculate the standard deviation for each gene

    num_nonzeros = np.count_nonzero(expdf.iloc[:,0:nsamps], axis=1) # Count the number of non zeros for each gene in the df
    expdf["num_samps_expressed"] = num_nonzeros # add the count to a new col in the df
    #print(num_nonzeros)
    
    # Subset df for only genes that appear in at least k samples (user-specified)
    #print(expdf.iloc[:,[len(expdf.columns)-1]] < args.number)
    samples_removed = len(expdf[expdf["num_samps_expressed"] < number])
    
    cutdf = expdf[expdf["num_samps_expressed"] >= number]
    cutdf = cutdf.drop(columns=['stdev', 'num_samps_expressed'])
    cutdf.columns = header_samps
    
    # Print summary statistics for the filtering of exp files
    dist = expdf["num_samps_expressed"].value_counts()
    dist = dist.rename_axis('Number of samples with expression').reset_index(name='Number of instances')
    dist = dist.sort_values(by='Number of samples with expression', ascending=False)
    #dist.index.name = None

    
    
    
    
    
    
    
    
    
    # Now filter the motif prior for only the TFs and target genes that are in the resulting filtered expression file
    mot = pd.read_csv(motif, sep='\t', header=None)
    mot.columns = ["TF", "target", "val"]
    
    remaining_genes = cutdf.index
    #print(remaining_genes)
    
    uniq_TFs = np.unique(mot["TF"]) # list of unique TFs in the motif prior
    uniq_targets = np.unique(mot["target"]) # list of unique TF targets in the motif prior
    
    mot_filtered_TFs = mot[mot["TF"].isin(remaining_genes)]
    
    mot_filtered_TFs_targets = mot_filtered_TFs[mot_filtered_TFs["target"].isin(remaining_genes)]

    # Print summary statistics for the filtering of motif files
    n_motif_TF = len(np.unique(mot["TF"])) 
    n_motif_target = len(np.unique(mot["target"])) 
    n_motif_TF_filtered = n_motif_TF - len(np.unique(mot_filtered_TFs_targets["TF"]))
    n_motif_TF_kept = len(np.unique(mot_filtered_TFs_targets["TF"]))

    n_motif_target_filtered = n_motif_target - len(np.unique(mot_filtered_TFs_targets["target"]))
    n_motif_target_kept = len(np.unique(mot_filtered_TFs_targets["target"]))
            
    uniq_remain_genes_motif = np.unique(list(mot_filtered_TFs_targets["target"]) + list(mot_filtered_TFs_targets["TF"]))

    save_file_path_motif = os.path.join(outdir, motoutfile + "_filtered.txt")

    mot_filtered_TFs_targets.to_csv(save_file_path_motif, sep = "\t", header = None, index = False)    
    
    
    
    
    
    
    
    # Then filter the ppi data for the same 
    ppi = pd.read_csv(ppi, sep='\t', header=None, index_col=None)
    ppi.columns = ["source", "targetTF", "edge_weight"]
    
    ppi_filtered_source = ppi[ppi["source"].isin(remaining_genes)]
    ppi_filtered_source_TFs = ppi_filtered_source[ppi_filtered_source["targetTF"].isin(remaining_genes)]
    
    # Print summary statistics for the filtering of ppi files
    n_ppi_source = len(np.unique(ppi["source"])) 
    n_ppi_target = len(np.unique(ppi["targetTF"])) 
    n_ppi_source_filtered = n_ppi_source - len(np.unique(ppi_filtered_source_TFs["source"]))
    n_ppi_source_kept = len(np.unique(ppi_filtered_source_TFs["source"]))

    n_ppi_target_filtered = n_ppi_target - len(np.unique(ppi_filtered_source_TFs["targetTF"]))
    n_ppi_target_kept = len(np.unique(ppi_filtered_source_TFs["targetTF"]))

    uniq_remain_genes_ppi = np.unique(list(ppi_filtered_source_TFs["targetTF"]) + list(ppi_filtered_source_TFs["source"]))

    save_file_path_ppi = os.path.join(outdir, ppioutfile + "_filtered.txt")

    ppi_filtered_source_TFs.to_csv(save_file_path_ppi, sep = "\t", header = None, index = False)    









    # Finally, filter the exp data for just the genes that are in the motif
    all_motif_TF_targets = list(uniq_TFs) + list(uniq_targets)
    
    cutdf_filtered_on_filtered_motif = cutdf[cutdf.index.isin(all_motif_TF_targets)]
    
    n_genes_removed_from_exp_based_on_motif = len(cutdf) - len(cutdf_filtered_on_filtered_motif)


    #print(cutdf_filtered_on_filtered_motif)
    cutdf_filtered_on_filtered_motif += 1 # add pseudocount of 1
    #print(cutdf_filtered_on_filtered_motif)    
    save_file_path_exp = os.path.join(outdir, expoutfile + "_filtered.txt")

    cutdf_filtered_on_filtered_motif.to_csv(save_file_path_exp, sep = "\t", columns = header_samps, header = None)

    
    
    
    print(f"Files saved:\n{save_file_path_motif}\n{save_file_path_ppi}\n{save_file_path_exp}")
    
    
    
    
    # Print summary of run
    print(dist)    
    
    print("\nExp file filtering info based on low abundance genes:")
    print(str(samples_removed) + " out of " + str(num_originalgenes) + " genes were removed that did not pass the -n threshold. " + str(len(cutdf)) + " genes remain.")
    
    print("\nMotif file filtering info:")
    print(str(n_motif_TF_filtered) + " out of " + str(n_motif_TF) + " TFs were filtered out")
    print(str(n_motif_target_filtered) + " out of " + str(n_motif_target) + " target genes were filtered out.")  
    print(str(len(uniq_remain_genes_motif)) + " unique genes remain in the motif.")
    
    print("\nPPI file filtering info:")
    print(str(n_ppi_source_filtered) + " out of " + str(n_ppi_source) + " source genes were filtered out and " + str(n_ppi_source_kept) + " genes remain.")
    print(str(n_ppi_target_filtered) + " out of " + str(n_ppi_target) + " target genes were filtered out and " + str(n_ppi_target_kept) + " genes remain.")    
    print(str(len(uniq_remain_genes_ppi)) + " unique genes remain in the PPI file.\n")
   
    print("Exp file filtering info based on only genes remaining in the motif file after filtering out TFs and genes for low abundance genes:")
    print(str(n_genes_removed_from_exp_based_on_motif) + " genes were filtered in this step. " + str(len(cutdf_filtered_on_filtered_motif)) + " genes remain in the filtered exp file.\n")

    save_file_path_filter_stats = os.path.join(outdir, expoutfile + "_filtering_statistics.txt")

    data_paths = {'exp': save_file_path_exp, 'motif': save_file_path_motif, 'ppi': save_file_path_ppi}    
    with open('./tmp/processed_data_paths.yml', 'w') as yaml_file:
        yaml.dump(data_paths, yaml_file, default_flow_style=False)
    
    with open(save_file_path_filter_stats, 'w') as sys.stdout:
        print("\n\nNumber of samples with expression | Number of instances")
        print(dist)    
        
        print("\nExp file filtering info based on low abundance genes:")
        print(str(samples_removed) + " out of " + str(num_originalgenes) + " genes were removed that did not pass the -n threshold. " + str(len(cutdf)) + " genes remain.")
        
        print("\nMotif file filtering info:")
        print(str(n_motif_TF_filtered) + " out of " + str(n_motif_TF) + " TFs were filtered out")
        print(str(n_motif_target_filtered) + " out of " + str(n_motif_target) + " target genes were filtered out.")  
        print(str(len(uniq_remain_genes_motif)) + " unique genes remain in the motif.")
        
        print("\nPPI file filtering info:")
        print(str(n_ppi_source_filtered) + " out of " + str(n_ppi_source) + " source genes were filtered out and " + str(n_ppi_source_kept) + " genes remain.")
        print(str(n_ppi_target_filtered) + " out of " + str(n_ppi_target) + " target genes were filtered out and " + str(n_ppi_target_kept) + " genes remain.")    
        print(str(len(uniq_remain_genes_ppi)) + " unique genes remain in the PPI file.\n")
    
        print("Exp file filtering info based on only genes remaining in the motif file after filtering out TFs and genes for low abundance genes:")
        print(str(n_genes_removed_from_exp_based_on_motif) + " genes were filtered in this step. " + str(len(cutdf_filtered_on_filtered_motif)) + " genes remain in the filtered exp file.\n")
    

    
