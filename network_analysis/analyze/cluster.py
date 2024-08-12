from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

def scale_data(df):
    '''
    Function that scales data to be used for PCA decomposition later

    Parameters 
    ----------
        df : panda data frame
            A df of just numeric values
    '''   
    scaler = StandardScaler()
    scaled = scaler.fit_transform(df)
    return(scaled)


def pca_fit_transform(scaled, nc, df_index):
    '''
    Function that performs a PCA decomposition on scaled data

    Parameters 
    ----------
        scaled : panda data frame
            A df of scaled values (output of scale_data())
        nc : int
            The number of components to keep from the PCA decomposition. Must match the same number used in pca_cum_var()
        df_index : list
            A list of sample names to give the output data frame as index
    '''
    
    pca = PCA(n_components = nc)
    pca_table = pca.fit_transform(scaled)
    # print("pca table:")
    # print(pca_table)
    colnames = [f"PC{x}" for x in range(1, nc+1)]
    pca_df = pd.DataFrame(data = pca_table, columns = colnames)
    pca_df = pca_df.set_index(df_index)
    pca_dict = {}
    pca_dict["df"] = pca_df
    pca_dict["pcs"] = pca_table
    return(pca_dict)


def km(df, kmax, numsamps):
    '''
    Function that performs kmeans clustering on data for multiple values of k
    
    Adapted from https://realpython.com/k-means-clustering-python/

    Parameters 
    ----------
        df : panda data frame
            A df of PCA values output by pca_fit_transform
        kmax : int
            The maximum value of k
        numsamps : int
            The number of samples in the data frame
    ''' 
    # First need to find the best number of clusters
    silhouette_scores = []

    kmeans_kwargs = {
        "init": "random",
        "n_init": 10,
        "max_iter": 300,
        "random_state": 42,
    }

    for k in range(2, kmax):
        kmeans = KMeans(n_clusters=k, **kmeans_kwargs)
        kmeans.fit(df)
        if len(np.unique(kmeans.labels_)) <= numsamps - 1: 
            score = silhouette_score(df, kmeans.labels_)
            silhouette_scores.append(score)

    print(silhouette_scores)
    
    return(silhouette_scores)