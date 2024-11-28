import os

def find_ex_paths():
    '''
    Description:
        This code finds the paths to the files in this folder so files can be copied to user's current working directory
     
    Parameters:
    -----------
        - None
    
    Returns:
    -----------
        - list of paths to be copied in main sisana script
    '''
    
    fnames = ["BRCA_TCGA_20_LumA_LumB_samps_5000_genes_exp.tsv",
              "BRCA_TCGA_20_LumA_LumB_samps_mapping.csv",
              "BRCA_TCGA_20_LumA_LumB_samps_survival_data.csv",
              "genes.txt",
              "motif_tcga_brca.tsv",
              "params.yml",
              "ppi_tcga_brca.tsv"]
    dir_path = os.path.dirname(os.path.realpath(__file__))
    files = [dir_path + "/" + x for x in fnames]
    return(files)