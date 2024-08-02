from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import umap
from sklearn.manifold import TSNE
import pandas as pd

def run_pca(df):
    """
    This function performs PCA dimensionality reduction on the expression or degrees, which then gets passed into the UMAP or t-SNE function
    
    Args:
        df: the data frame to perform PCA on

    Returns:
        df_pca: a data frame with the resulting primary components
    """
    
    print("Performing PCA on input data...")
 
    # First need to scale the data
    scaler = StandardScaler()
    scaled_df = scaler.fit_transform(df)
    
    # Then perform PCA
    model = PCA(n_components=2)
    principal_components = model.fit_transform(scaled_df)
    
    plt.figure(figsize=(8,6))
    plt.scatter(principal_components[:,0],principal_components[:,1])
    plt.xlabel('PC 1')
    plt.ylabel('PC 2')

    # choose number of components
    explained_variance = []

    for n in range(1,len(df.columns)):
        pca = PCA(n_components=n)
        pca.fit(scaled_X)
        
        explained_variance.append(np.sum(pca.explained_variance_ratio_))
        
    plt.plot(range(1,len(df.columns)),explained_variance)
    plt.xlabel("Number of Components")
    plt.ylabel("Variance Explained")

def run_umap(df):
    """
    This function performs UMAP dimensionality reduction on the expression or degrees
    
    Args:
        df: the data frame to perform clustering on

    Returns:
        df_umap: a data frame with the resulting UMAP calculations
    """
    
    
def run_tsne(df):


fname = "/storage/kuijjerarea/nolan/sisana/tests/pca_umap_tsne/input/lioness_df_indegree.csv"
meta = "/storage/kuijjerarea/nolan/sisana/tests/pca_umap_tsne/input/brca_samples_subtype_info.csv"
pca = True
choice = "tsne"

datafile = pd.read_csv(fname, engine = "pyarrow", header = None)

if pca:
    datafile = run_pca(datafile)

if choice == "tsne":
    run_tsne(datafile)
elif choice == "umap":
    run_umap(datafile)
