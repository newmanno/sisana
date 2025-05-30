import requests
import os
from tqdm import tqdm
import sys

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

    curfile = 1

    for link in urls:
        filename = link.split("/")[-1]
        print(f"Downloading file {curfile} of {len(urls)}: {filename}")
        
        r = requests.get(link, stream=True)
        fname = os.path.basename(link)
        
        # Sizes in bytes.
        total_size = int(r.headers.get("content-length", 0))
        block_size = 1024       
        
        with tqdm(total=total_size, unit="B", unit_scale=True) as progress_bar:
            with open(f'./example_inputs/{fname}', 'wb') as f:
                for chunk in r.iter_content(block_size):
                    progress_bar.update(len(chunk))
                    f.write(chunk)
        print("\n")
        curfile += 1
        f.close()