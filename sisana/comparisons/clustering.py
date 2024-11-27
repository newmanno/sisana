from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import argparse
from analyze import file_to_list, scale_data, pca_fit_transform, km, assign_to_cluster, plot_pca, NotASubsetError
import pandas as pd
import os 
import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':
    """
    Description:
        This code creates a volcano plot from the p-value/FDR and log2 fold change
    """

    parser = argparse.ArgumentParser(description="Example command: python clustering.py -d indegree.csv -f csv -o ./output/")
    ArgGroup = parser.add_argument_group('Required arguments')  
    
    ArgGroup.add_argument("-d", "--datafile", type=str, help="Path to file containing the indegrees or expression data. If using expression data, file should not have a header and instead sample order should be supplied with --sampnames", required=True)
    ArgGroup.add_argument("-f", "--filetype", choices = ["csv", "txt"], type=str, help="File type of the input data file", required=True)
    ArgGroup.add_argument("-c", "--compnum", type=int, help="Maximum value of components to use in PCA. Must be less than or equal to number of samples", required=True)
    ArgGroup.add_argument("-k", "--clusterchoice", choices = ["kmeans", "hkmeans", "hierarchical", "consensus"], type=str, help="Type of clustering to perform, either basic kmeans or kmeans using Hartigan's method (hkmeans)", required=True)
    ArgGroup.add_argument("-n", "--numclus", type=int, help="Number of clusters to attempt to separate samples into", required=False)
    ArgGroup.add_argument("-s", "--sampnames", type=str, help="Path to file containing the order of samples in the expression data", required=False)
    ArgGroup.add_argument("-p", "--plotpca", action='store_true', help="Flag; Do you want to create a PCA plot of the first two primary components?", required=False)
    ArgGroup.add_argument("-m", "--metadata", type=str, help="Metadata file with a header in the format 'name,group', where name is the sample name and group is the category of the sample (e.g. Control, Disease, etc.). Required if --plotpca is given as an argument.", required=False)
    ArgGroup.add_argument("-o", "--outdir", type=str, help="Path to directory to output results to", required=True) 
    
    args = parser.parse_args()

    # print input arguments to stdout for reference later
    print(f"Supplied data file: {args.datafile}")
    print(f"Maximum number of components for PCA: {args.compnum}")
    print(f"Clustering method choice: {args.clusterchoice}")
    if args.numclus is not None: print(f"Number of clusters to separate data into: {args.numclus}") 
    if args.sampnames is not None: print(f"Supplied sample name file: {args.numclus}") 
    if args.metadata is not None: print(f"Supplied metadata file: {args.numclus}") 

    if args.filetype == "csv":
        data = pd.read_csv(args.datafile, index_col = 0)
    elif args.filetype == "txt":
        data = pd.read_csv(args.datafile, index_col = 0, sep = "\t")

    # If user is using expression data, add the sample names
    if args.sampnames is not None:
        sampsfile = open(args.sampnames, "r")
        fileread = sampsfile.read()
        namelist = fileread.split("\n") 
        namelist = list(filter(None, namelist))
            
        data.columns = namelist

    # kmeans code cannot be run if user inputs too small of value for the maximum components
    if args.compnum > len(data.columns):
        raise Exception("Error: Please ensure the --compnum value supplied is larger than the number of samples in your dataset.")

    # Compute PCA
    data_t = data.T
    data_scaled = scale_data(data_t)
    PCA_out = pca_fit_transform(data_scaled, args.compnum, data_t.index) 
    PCA_df = PCA_out["df"] 
    PCA_pcs = PCA_out["pcs"] 
    
    if args.plotpca:
        # Plot first two components
        pca_df = pd.DataFrame()
        pca_df["PC1"] = [x[0] for x in PCA_pcs]
        pca_df["PC2"] = [x[1] for x in PCA_pcs]
        pca_df.index = data_t.index
        # print(pca_df)

        # Add metadata for plotting
        meta = pd.read_csv(args.metadata)
        meta_cols = meta.columns
        
        try:
            meta_w_index = meta.set_index("name")
        except KeyError:
            print("\nError: Metadata file could not be read properly. Please make sure the metadata file has the header 'name,group'")
            exit(0)
            
        try:
            meta_pc = pd.concat([PCA_df, meta_w_index], axis = 1) # merge all three tables based on their indeces
        except pd.errors.InvalidIndexError:
            raise NotASubsetError(meta_w_index.index, PCA_df.index, "samples")
            exit(0)
        # print(meta_pc)

        # plt.scatter(pca_df["PC1"], pca_df["PC2"])
        # plt.xlabel("PC1")
        # plt.ylabel("PC2")
        
        if not os.path.exists(args.outdir):
            os.makedirs(args.outdir)    
                
        outpath = os.path.join(args.outdir, f"pca.png")

        if meta_cols[0] == "name" and meta_cols[1] == "group": # This check helps ensure no uninterpretable error is thrown
            plot_pca(meta_pc, "group", outpath)
        else:
            raise Exception("Error: Please ensure your metadata file has the header 'name,group'")

    # Assign samples to cluster, calculating the optimal kmax if needed
    print("Performing clustering. Please wait...")
    
    if args.clusterchoice == "consensus": # do not find optimal k, since the ConsensusCluster class already does this
        assigned_samps = assign_to_cluster(PCA_pcs, data_t, 11, args.clusterchoice)
        uniq_clusters = len(np.unique(assigned_samps["cluster"]))
        assigned_samps_fname = os.path.join(args.outdir, f"{args.clusterchoice}_sample_cluster_assignment_{uniq_clusters}_clusters.txt")
                
    # if user did not specify the number of clusters, find the optimal number based on silhouette score and perform clustering using the resulting k value
    elif args.numclus is None: 
        # Choose a value for kmax. Normally we will use 11, but if user has 11 samples or less, this will not work so we need to adjust kmax to be one less than the number of samples
        if len(data.columns) <= 5:
            raise Exception("Error: Must have at least 6 samples for finding clusters")
        elif len(data.columns) <= 10:
            kmax = len(data.columns) - 1
        else:
            kmax = 10
        
        # Calculate silhouette scores to find optimal number of clusters
        print(f"Calculated kmax value to use in clustering: {kmax}")
        sil_scores = km(PCA_pcs, kmax, len(data.columns), args.clusterchoice)

        # Report results to user
        print(f"\nNumber of clusters : Silhouette score")
        for i in range(2, len(sil_scores)+2):
            print(f"{i}: {sil_scores[i-2]:.2f}")
        
        largest_score = max(sil_scores)
        num_clusters = sil_scores.index(largest_score) + 2

        if largest_score > 0.4:
            print(f"\nThe largest silhouette score of {largest_score:.2f} was found for {num_clusters} clusters. Thus, cluster assignment will be performed for {num_clusters} clusters.")
        else:
            import warnings
            warnings.warn("\nWARNING: Silhouette score is small!", RuntimeWarning)
            print(f"The largest silhouette score of {largest_score:.2f} was found for {num_clusters} clusters... Considering this value is smaller than 0.4, you may want to try alternative clustering methods. Regardless, cluster assignment will be performed for {num_clusters} clusters.")
        
        # Assign samples to clusters based on the k value with highest silhouette score
        assigned_samps = assign_to_cluster(PCA_pcs, data_t, num_clusters, args.clusterchoice)
        
        
        assigned_samps_fname = os.path.join(args.outdir, f"{args.clusterchoice}_sample_cluster_assignment_{num_clusters}_clusters.txt")
                
        # Plot resulting silhouette scores
        plt.style.use("fivethirtyeight")
        plt.plot(range(2, kmax + 1), sil_scores)
        plt.xticks(range(2, kmax + 1))
        plt.xlabel("Number of Clusters")
        plt.ylabel("Silhouette Coefficient")
        plt.show() 
        
        outname = os.path.join(args.outdir, f"kmeans_silhouette_scores.png")
        plt.savefig(outname, bbox_inches='tight')
        
        print(f"Calculated silhouette scores for each number of clusters saved here: {outname}")     

    else:
        assigned_samps = assign_to_cluster(PCA_pcs, data_t, args.numclus, args.clusterchoice)
        assigned_samps_fname = os.path.join(args.outdir, f"{args.clusterchoice}_sample_cluster_assignment_{args.numclus}_clusters.txt")
    
    if len(assigned_samps) <= 20: # too messy to print to stdout if there are too many samples        
        print("\nSample cluster assignment can be seen below.")
        print(assigned_samps.to_csv(sep='\t', index=False))
    
    assigned_samps.to_csv(assigned_samps_fname, sep = "\t", index = False)   
    
    print(f"\nSample cluster assignment saved here: {assigned_samps_fname}")     
    
    if args.plotpca:
        print(f"PCA plot for scaled input data saved here: {outpath}\n")     
