# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import argparse
import sys
import os
from pathlib import Path
import yaml

def preprocess_data(exp: str, number: int, outdir: str):
    """
    Description:
        This code performs a survival analysis between two user-defined groups and outputs
        both the survival plot and the statistics for the comparison(s)
        
    Parameters:
    -----------\
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
        - str of the output file location 
    """
    
    # Create output file prefix by removing the .txt suffix
    expoutfile = Path(exp).stem
    
    # Create output directory if one does not already exist
    os.makedirs(outdir, exist_ok=True)
        
    expdf = pd.read_csv(exp, sep='\t', index_col=0)  
    
    with open("./tmp/samples.txt", "w") as f:
        for col in expdf.columns:
            f.write(col + "\n")
      
    # num_originalgenes = len(expdf)
    nsamps = len(expdf.columns)

    num_nonzeros = np.count_nonzero(expdf.iloc[:,0:nsamps], axis=1) # Count the number of non zeros for each gene in the df
    expdf["num_samps_expressed"] = num_nonzeros # add the count to a new col in the df
    
    # Subset df for only genes that appear in at least k samples (user-specified)
    samples_removed = len(expdf[expdf["num_samps_expressed"] < number])
    
    cutdf = expdf[expdf["num_samps_expressed"] >= number]
    cutdf = cutdf.drop(columns=['num_samps_expressed'])

    cutdf.columns = expdf.columns[:-1]
    
    # Print summary statistics for the filtering of exp files
    dist = expdf["num_samps_expressed"].value_counts()
    dist = dist.rename_axis('Number of samples with expression').reset_index(name='Number of instances')
    dist = dist.sort_values(by='Number of samples with expression', ascending=False)
    dist['Removed?'] = dist['Number of samples with expression'].apply(lambda x: 'yes' if x < number else 'no')
    print(dist.to_string(index=False))
    sum_removed = dist[dist['Removed?'] == 'yes']['Number of instances'].sum()
    print(f"Number of genes removed in total: {sum_removed}")

    basename = f"{expoutfile}_preprocessed.txt"
    file_outloc = os.path.join(outdir, basename)

    cutdf.to_csv(file_outloc, sep = "\t")    
    print(f"\nFile saved: {file_outloc}")
    
    return(file_outloc)
    
    
    
    
    
    
    # # Now filter the motif prior for only the TFs and target genes that are in the resulting filtered expression file
    # mot = pd.read_csv(motif, sep='\t', header=None)
    # mot.columns = ["TF", "target", "val"]
    
    # remaining_genes = cutdf.index
    # #print(remaining_genes)
    
    # uniq_TFs = np.unique(mot["TF"]) # list of unique TFs in the motif prior
    # uniq_targets = np.unique(mot["target"]) # list of unique TF targets in the motif prior
    
    # mot_filtered_TFs = mot[mot["TF"].isin(remaining_genes)]
    
    # mot_filtered_TFs_targets = mot_filtered_TFs[mot_filtered_TFs["target"].isin(remaining_genes)]

    # # Print summary statistics for the filtering of motif files
    # n_motif_TF = len(np.unique(mot["TF"])) 
    # n_motif_target = len(np.unique(mot["target"])) 
    # n_motif_TF_filtered = n_motif_TF - len(np.unique(mot_filtered_TFs_targets["TF"]))
    # n_motif_TF_kept = len(np.unique(mot_filtered_TFs_targets["TF"]))

    # n_motif_target_filtered = n_motif_target - len(np.unique(mot_filtered_TFs_targets["target"]))
    # n_motif_target_kept = len(np.unique(mot_filtered_TFs_targets["target"]))
            
    # uniq_remain_genes_motif = np.unique(list(mot_filtered_TFs_targets["target"]) + list(mot_filtered_TFs_targets["TF"]))

    # save_file_path_motif = os.path.join(outdir, motoutfile + "_filtered.txt")

    # mot_filtered_TFs_targets.to_csv(save_file_path_motif, sep = "\t", header = None, index = False)    
    
    
    
    
    
    
    
    # # Then filter the ppi data for the same 
    # ppi = pd.read_csv(ppi, sep='\t', header=None, index_col=None)
    # ppi.columns = ["source", "targetTF", "edge_weight"]
    
    # ppi_filtered_source = ppi[ppi["source"].isin(remaining_genes)]
    # ppi_filtered_source_TFs = ppi_filtered_source[ppi_filtered_source["targetTF"].isin(remaining_genes)]
    
    # # Print summary statistics for the filtering of ppi files
    # n_ppi_source = len(np.unique(ppi["source"])) 
    # n_ppi_target = len(np.unique(ppi["targetTF"])) 
    # n_ppi_source_filtered = n_ppi_source - len(np.unique(ppi_filtered_source_TFs["source"]))
    # n_ppi_source_kept = len(np.unique(ppi_filtered_source_TFs["source"]))

    # n_ppi_target_filtered = n_ppi_target - len(np.unique(ppi_filtered_source_TFs["targetTF"]))
    # n_ppi_target_kept = len(np.unique(ppi_filtered_source_TFs["targetTF"]))

    # uniq_remain_genes_ppi = np.unique(list(ppi_filtered_source_TFs["targetTF"]) + list(ppi_filtered_source_TFs["source"]))

    # save_file_path_ppi = os.path.join(outdir, ppioutfile + "_filtered.txt")

    # ppi_filtered_source_TFs.to_csv(save_file_path_ppi, sep = "\t", header = None, index = False)    









    # # Finally, filter the exp data for just the genes that are in the motif
    # all_motif_TF_targets = list(uniq_TFs) + list(uniq_targets)
    
    # cutdf_filtered_on_filtered_motif = cutdf[cutdf.index.isin(all_motif_TF_targets)]
    
    # n_genes_removed_from_exp_based_on_motif = len(cutdf) - len(cutdf_filtered_on_filtered_motif)


    # #print(cutdf_filtered_on_filtered_motif)
    # cutdf_filtered_on_filtered_motif += 1 # add pseudocount of 1
    # #print(cutdf_filtered_on_filtered_motif)    
    # save_file_path_exp = os.path.join(outdir, expoutfile + "_filtered.txt")

    # cutdf_filtered_on_filtered_motif.to_csv(save_file_path_exp, sep = "\t", columns = header_samps, header = None)

    
    
    
    # print(f"Files saved:\n{save_file_path_motif}\n{save_file_path_ppi}\n{save_file_path_exp}")
    
    
    
    
    # # Print summary of run
    # print(dist)    
    
    # print("\nExp file filtering info based on low abundance genes:")
    # print(str(samples_removed) + " out of " + str(num_originalgenes) + " genes were removed that did not pass the -n threshold. " + str(len(cutdf)) + " genes remain.")
    
    # print("\nMotif file filtering info:")
    # print(str(n_motif_TF_filtered) + " out of " + str(n_motif_TF) + " TFs were filtered out")
    # print(str(n_motif_target_filtered) + " out of " + str(n_motif_target) + " target genes were filtered out.")  
    # print(str(len(uniq_remain_genes_motif)) + " unique genes remain in the motif.")
    
    # print("\nPPI file filtering info:")
    # print(str(n_ppi_source_filtered) + " out of " + str(n_ppi_source) + " source genes were filtered out and " + str(n_ppi_source_kept) + " genes remain.")
    # print(str(n_ppi_target_filtered) + " out of " + str(n_ppi_target) + " target genes were filtered out and " + str(n_ppi_target_kept) + " genes remain.")    
    # print(str(len(uniq_remain_genes_ppi)) + " unique genes remain in the PPI file.\n")
   
    # print("Exp file filtering info based on only genes remaining in the motif file after filtering out TFs and genes for low abundance genes:")
    # print(str(n_genes_removed_from_exp_based_on_motif) + " genes were filtered in this step. " + str(len(cutdf_filtered_on_filtered_motif)) + " genes remain in the filtered exp file.\n")

    # save_file_path_filter_stats = os.path.join(outdir, expoutfile + "_filtering_statistics.txt")

    # data_paths = {'exp': save_file_path_exp, 'motif': save_file_path_motif, 'ppi': save_file_path_ppi}    
    # with open('./tmp/processed_data_paths.yml', 'w') as yaml_file:
    #     yaml.dump(data_paths, yaml_file, default_flow_style=False)
    
    # with open(save_file_path_filter_stats, 'w') as sys.stdout:
    #     print("\n\nNumber of samples with expression | Number of instances")
    #     print(dist)    
        
    #     print("\nExp file filtering info based on low abundance genes:")
    #     print(str(samples_removed) + " out of " + str(num_originalgenes) + " genes were removed that did not pass the -n threshold. " + str(len(cutdf)) + " genes remain.")
        
    #     print("\nMotif file filtering info:")
    #     print(str(n_motif_TF_filtered) + " out of " + str(n_motif_TF) + " TFs were filtered out")
    #     print(str(n_motif_target_filtered) + " out of " + str(n_motif_target) + " target genes were filtered out.")  
    #     print(str(len(uniq_remain_genes_motif)) + " unique genes remain in the motif.")
        
    #     print("\nPPI file filtering info:")
    #     print(str(n_ppi_source_filtered) + " out of " + str(n_ppi_source) + " source genes were filtered out and " + str(n_ppi_source_kept) + " genes remain.")
    #     print(str(n_ppi_target_filtered) + " out of " + str(n_ppi_target) + " target genes were filtered out and " + str(n_ppi_target_kept) + " genes remain.")    
    #     print(str(len(uniq_remain_genes_ppi)) + " unique genes remain in the PPI file.\n")
    
    #     print("Exp file filtering info based on only genes remaining in the motif file after filtering out TFs and genes for low abundance genes:")
    #     print(str(n_genes_removed_from_exp_based_on_motif) + " genes were filtered in this step. " + str(len(cutdf_filtered_on_filtered_motif)) + " genes remain in the filtered exp file.\n")
    

    
