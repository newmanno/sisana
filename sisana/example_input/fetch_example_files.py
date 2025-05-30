import requests
import os

def fetch_files():
    """
    Description:
        This code fetches the example input files from Zenodo and saves them in a directory called ./example_inputs/.
        
    Parameters:
    -----------
        - None
    
    Returns:
    -----------
        - Nothing
    """
    
    os.makedirs('./example_inputs/', exist_ok=True)

    urls = ['https://zenodo.org/records/15552563/files/BRCA_TCGA_200_LumA_LumB_samps_mapping_w_header.csv',
            'https://zenodo.org/records/15552563/files/BRCA_TCGA_200_LumA_LumB_samps_survival_data.csv',
            'https://zenodo.org/records/15552563/files/BRCA_TCGA_20_LumA_LumB_samps_5000_genes_exp.tsv',
            'https://zenodo.org/records/15552563/files/c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt',
            'https://zenodo.org/records/15552563/files/c2.cp.reactome.v2023.2.Hs.symbols.gmt',
            'https://zenodo.org/records/15552563/files/genes_to_extract.txt',
            'https://zenodo.org/records/15552563/files/Hallmark.v2023.2.Hs.symbols.gmt',
            'https://zenodo.org/records/15552563/files/heatmap_genes.txt',
            'https://zenodo.org/records/15552563/files/heatmap_genes.txt',
            'https://zenodo.org/records/15552563/files/lioness_df_indegree_3_decimal_places_subset_200_LumALumB_samps.csv',
            'https://zenodo.org/records/15552563/files/params.yml',
            'https://zenodo.org/records/15552563/files/quantity_plot_genes.txt',
            'https://zenodo.org/records/15552563/files/volcano_plot_genes.txt',
            
            ### From the SPONGE Zenodo repo:
            'https://zenodo.org/records/13628785/files/ppi_prior_2024.tsv',
            'https://zenodo.org/records/13628785/files/motif_prior_names_2024.tsv']

    for link in urls:
        r = requests.get(link, stream=True)
        fname = os.path.basename(link)

        with open(f'./example_inputs/{fname}', 'wb') as f:
            for chunk in r.iter_content(chunk_size=16*1024):
                f.write(chunk)
        f.close()