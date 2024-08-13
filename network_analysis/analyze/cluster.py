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
    Function that performs a PCA decomposition on scaled data and returns a dictionary of the calculated components

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


def km(pca_df, kmax, numsamps, clusttype):
    '''
    Function that performs kmeans clustering on data for multiple values of k and returns the silhouette scores for each
    
    Adapted from https://realpython.com/k-means-clustering-python/

    Parameters 
    ----------
        pca_df : panda data frame
            A df of PCA values output by pca_fit_transform
        kmax : int
            The maximum value of k
        numsamps : int
            The number of samples in the data frame
        clusttype : str [choices are kmeans or hkmeans]
            The type of clustering to perform. kmeans performs standard kmeans clustering while hkmeans uses the Hartigan method 
    ''' 
    
    silhouette_scores = []

    if clusttype == "kmeans":
        kmeans_kwargs = {
            "init": "random",
            "n_init": 10,
            "max_iter": 300,
            "random_state": 42,
        }

        for k in range(2, kmax):
            kmeans = KMeans(n_clusters = k, **kmeans_kwargs)
            kmeans.fit(pca_df)
            
            if len(np.unique(kmeans.labels_)) <= numsamps - 1: 
                score = silhouette_score(pca_df, kmeans.labels_)
                silhouette_scores.append(score)
                
    elif clusttype == "hkmeans":
        from sklearn.feature_extraction.text import TfidfVectorizer
        from sklearn import metrics
        from hkmeans import HKMeans
        
        for k in range(2, kmax):
            hkmeans = HKMeans(n_clusters = k, random_state = 128, n_init = 10,
                            n_jobs = 16, max_iter = 15)
            hkmeans.fit(pca_df)
            
            if len(np.unique(hkmeans.labels_)) <= numsamps - 1: 
                score = silhouette_score(pca_df, hkmeans.labels_)
                silhouette_scores.append(score)
    
    return(silhouette_scores)


def assign_to_cluster(pca_df, data_t, k, clusttype):
    '''
    Function that assigns the samples to clusters and returns a data frame containing the cluster assignment for each sample
    
    Parameters 
    ----------
        pca_df : panda data frame
            A df of PCA values that were calculated by pca_fit_transform
        data_t : panda data frame
            A transposed data frame, where samples are rows and variables (e.g. genes) are columns.
        k : int
            The value of k that corresponds to the maximum silhouette score
        clusttype : str [choices are kmeans or hkmeans]
            The type of clustering to perform. kmeans performs standard kmeans clustering while hkmeans uses the Hartigan method 
    '''     
    cluster_map = pd.DataFrame()
    cluster_map['sample'] = data_t.index.values

    if clusttype == "kmeans":
        kmeans_kwargs = {
            "init": "random",
            "n_init": 10,
            "max_iter": 300,
            "random_state": 42,
        }
        
        km = KMeans(n_clusters=k, **kmeans_kwargs).fit(pca_df)
        cluster_map['cluster'] = km.labels_

    elif clusttype == "hkmeans":
        from sklearn.feature_extraction.text import TfidfVectorizer
        from sklearn import metrics
        from hkmeans import HKMeans
        
        hkmeans = HKMeans(n_clusters = k, random_state = 128, n_init = 10,
                        n_jobs = 16, max_iter = 15)
        hkmeans.fit(pca_df)
        cluster_map['cluster'] = hkmeans.labels_

    print("\nBelow, the first column is the cluster number and the second column is the number of samples assigned to that cluster.")
    print(cluster_map['cluster'].value_counts())
    
    return(cluster_map)
