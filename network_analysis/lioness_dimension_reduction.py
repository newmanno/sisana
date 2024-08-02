import numpy as np
import pandas as pd
import argparse
import sys
import pickle
from pathlib import Path
from post import save_results
import time
import os 

if __name__ == '__main__':
    """
    Description:
        This code performs dimensionality reduction on the lioness dataframe. Users have the option of first performing PCA then running UMAP and HDBScan to find silhouette scores,
        or skipping the PCA step 
    """
    
    parser = argparse.ArgumentParser(description="Example command: python lioness_dimension_reduction.py -p <file.pickle> -m <metadata.tsv> -o ./output")
    reqArgGroup = parser.add_argument_group('Required arguments')  
    optArgGroup = parser.add_argument_group('Optional arguments')  
    
    # Required args
    reqArgGroup.add_argument("-p", "--picklefile", type=str, help="Path to pickle file created by lioness_to_pickle_df.py script", required=True)
    reqArgGroup.add_argument("-m", "--clinmeta", type=str, help="Clinical metadata for the patients, used for plotting", required=True)         
    reqArgGroup.add_argument("-m", "--meta", type=str, help="Clinical metadata for the patients, used for plotting", required=True)         
    reqArgGroup.add_argument("-b", "--colorby", type=str, help="Which metadata column to use for coloring plots", required=True)     
    reqArgGroup.add_argument("-o", "--outdir", type=str, help="Path to directory to output files to", required=True) 
    
    # Optional args    
    optArgGroup.add_argument("--pca", action='store_true', help="Flag; used if you wish to perform PCA prior to UMAP")
    optArgGroup.add_argument("--commonsamps", type = str, help="A two column csv file with a name of the ")

    args = parser.parse_args()
    
    try: 
        motfile = pd.read_csv(args.metadata, sep = "\t", engine = "pyarrow", header = None)
    except:
        raise Exception("There was an error reading in the metadata. Please make sure the metadata file is tab-delimited.")
        sys.exit()

    # Then load in the expression and metadata and filter the data for just the unique samples here,
    # then subset metadata for those samples as well

    common = "C:/Users/nolankn/Box/Kuijjer Lab/Kuijjer lab/Projects/Breast cancer subtypes/data/unique_samps_in_common_bw_data_and_metadata.csv"
    commondf = pd.read_csv(common)
    commonsamps = commondf["sampname"]
    commonshort = commondf.index

    # Filter the input data for the unique samples
    brca_data = pd.read_pickle(args.picklefile)
    brca_data_filtered = brca_data[brca_data.columns.intersection(commonshort)]
    
    # Scale the data (subtract the mean of each feature and divide by the stdev of each feature)
    # Example taken from https://scikit-learn.org/stable/modules/preprocessing.html
    from sklearn import preprocessing
    import numpy as np
    brca_data_filtered_scaled_np_array = brca_data_filtered.to_numpy()
    brca_data_filtered_scaled_np_array = preprocessing.StandardScaler().fit(brca_data_filtered_scaled_np_array)
    print(brca_data_filtered_scaled_np_array)

    # Read in and filter metadata
    metadata = "C:/Users/nolankn/Box/Kuijjer Lab/Kuijjer lab/Projects/Breast cancer subtypes/data/brca_samples_subtype_info.csv"
    meta = read.csv(metadata, header = T, row.names = 1)
    meta_filtered = meta[row.names(meta) %in% commonshort, ]

    pc <- prcomp(brca_data_filtered,
                center = TRUE,
                scale. = TRUE)
    rot = pc[["rotation"]]
    ggplot(rot[,1:2], aes(x=PC1, y=PC2)) + geom_point()
    autoplot(pc, data = dataClin, colour = 'BRCA_Subtype_PAM50')
    
    # Implementing PCA calculation: https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html
    # Implementing UMAP calculaiton: https://umap-learn.readthedocs.io/en/latest/