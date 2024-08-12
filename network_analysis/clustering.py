from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import argparse
from analyze import file_to_list, scale_data, pca_fit_transform, km
import pandas as pd
import os 

# Load in data
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

    if args.sampnames is not None:
        # Lioness file does not have any header or column names, needs them for t-test later
        sampsfile = open(args.sampnames, "r")
        fileread = sampsfile.read()
        namelist = fileread.split("\n") 
        namelist = list(filter(None, namelist))
            
        data.columns = namelist

    if args.compnum > len(data.columns):
        raise Exception("Error: Please ensure the kmax value supplied is larger than the number of samples in your dataset.")

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
        
    print(f"calculated kmax value to use in kmeans calculation: {kmax}")
    print(PCA_pcs)
    sil_scores = km(PCA_pcs, kmax, len(data.columns))
        
    import matplotlib.pyplot as plt
    plt.style.use("fivethirtyeight")
    plt.plot(range(2, kmax), sil_scores)
    plt.xticks(range(2, kmax))
    plt.xlabel("Number of Clusters")
    plt.ylabel("Silhouette Coefficient")
    plt.show() 
    
    outname = os.path.join(args.outdir, f"kmeans_silhouette_scores.png")
    plt.savefig(outname, bbox_inches='tight')

    print(f"File saved: {outname}")     
