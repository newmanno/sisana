from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import argparse
from analyze import file_to_list, scale_data, pca_fit_transform, km, assign_to_cluster
import pandas as pd
import os 
import matplotlib.pyplot as plt

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
    ArgGroup.add_argument("-s", "--sampnames", type=str, help="Path to file containing the order of samples in the expression data", required=False)
    ArgGroup.add_argument("-o", "--outdir", type=str, help="Path to directory to output results to", required=True) 
    
    args = parser.parse_args()

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

    # Choose a value for kmax. Normally we will use 11, but if user has 11 samples or less, this will not work so we need to adjust kmax 
    if len(data.columns) <= 2:
        raise Exception("Error: Must have at least 3 samples for finding clusters")
    elif len(data.columns) <= 10:
        kmax = len(data.columns) - 1
    else:
        kmax = 10
        
    # Perform kmeans clustering
    print(f"Calculated kmax value to use in kmeans calculation: {kmax}")
    sil_scores = km(PCA_pcs, kmax, len(data.columns))
            
    # Report results to user
    print(f"\nNumber of clusters : Silhouette score")
    
    for i in range(2, len(sil_scores)+2):
        print(f"{i}: {sil_scores[i-2]:.2f}")
    
    largest_score = max(sil_scores)
    num_clusters = sil_scores.index(largest_score) + 2
    
    if largest_score > 0.4:
        print(f"\nThe largest silhouette score of {largest_score:.2f} was found for {num_clusters} clusters.")
    else:
        print(f"\nThe largest silhouette score of {largest_score:.2f} was found for {num_clusters} clusters... Considering this value is smaller than 0.4, you may want to try alternative clustering methods.")
    
    # Assign samples to clusters based on the k value with highest silhouette score
    assigned_samps = assign_to_cluster(PCA_pcs, data_t, num_clusters)

    assigned_samps_fname = os.path.join(args.outdir, f"kmeans_sample_cluster_assignment_{num_clusters}_clusters.txt")
    assigned_samps.to_csv(assigned_samps_fname, sep = "\t", index = False)   
    print(f"\nFile saved: {assigned_samps_fname}")     
        
    # Plot resulting silhouette scores
    plt.style.use("fivethirtyeight")
    plt.plot(range(2, kmax), sil_scores)
    plt.xticks(range(2, kmax))
    plt.xlabel("Number of Clusters")
    plt.ylabel("Silhouette Coefficient")
    plt.show() 
    
    outname = os.path.join(args.outdir, f"kmeans_silhouette_scores.png")
    plt.savefig(outname, bbox_inches='tight')

    print(f"\nFile saved: {outname}")     
