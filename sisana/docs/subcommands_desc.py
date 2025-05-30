preprocess_desc = """
-----------------------------------
Step 1 - Preprocess
-----------------------------------

Preprocess is the first step to running SiSaNA. In this step, SiSaNA removes genes that 
are not expressed in at least the specified number of samples. 

Example: sisana preprocess params.yml

Input files (specified in the params.yml file): 
  - exp_file: A gene expression file with sample names as columns and gene names as rows

Output files: 
  - *_preprocessed.txt: The gene expression file, filtered only to contain genes that are 
  expressed in at least the user-defined number of samples. E.g. if the user sets "5" as 
  their value, then SiSaNA will remove genes not expressed in at least 5 samples.
"""

generate_desc = """
-----------------------------------
Step 2 - Generate
-----------------------------------

Generate is the second step in the SiSaNA pipeline. In this step, the networks
(both PANDA and per-sample LIONESS networks) are reconstructed. Then, the networks
get saved as a .pickle format in the ./tmp/ directory for quicker loading into 
downstream steps. SiSaNA will also calculate the indegrees and outdegrees for each
gene and transcription factor respectively, per sample.

Example: sisana generate params.yml

Input files (specified in the params.yml file):
  - exp: ./output/preprocess/BRCA_TCGA_20_LumA_LumB_samps_5000_genes_exp_preprocessed.txt # Path to the expression file
  - motif: ./example_inputs/motif_tcga_brca.tsv # Path to the motif prior file
  - ppi: ./example_inputs/ppi_tcga_brca.tsv # Path to the PPI prior file

Output files: 
  - panda_output.txt: saved to the user-specified output directory. Four column, tab-delimited
    file. Column 1: TF; Column 2: Gene; Column 3: Whether the interaction was in the prior motif;
    Column 4: Edge weight
  - lioness_df.pickle: saved to the ./tmp/ directory (not human-readable)
  - lioness.npy: saved to the user-specified output directory, which is a numpy formatted
    file of the lioness dataframe, where rows are TF-gene interactions and columns are samples
  - *_indegree.csv and *outdegree.csv: contain the indegrees for each gene
    per-sample and the outdegrees for each TF per-sample, saved to the user-specified output directory
    
"""

compare_desc = """
-----------------------------------
Step 3 - Compare
-----------------------------------

Compare is the third step in the SiSaNA pipeline. In this step, the samples can be compared 
in a number of different ways.

  ----- Option 1: means -----
  Compares the means between two sample groups. This performs a t-test or Mann-Whitney
  test (or their paired alternatives) between two groups. Users can compare differences
  in indegree/outdegree or expression at this step.
  
    Example: sisana compare means params.yml
  
    - Input files:
      - datafile: Data file (either expression, indegree, or outdegree) in tsv or csv format
      - mapfile: CSV mapping file, which maps sample name (column 1) to sample group (column 2). Assumes file has a header.
      
    - Output files:
      - A .rnk file that is in rnk format, ranked on test statistic, for use with "siana compare gsea"
      - A .txt file containing all the statistical outputs for the analysis
      
  ----- Option 2: survival -----
  Performs a survival analysis between two groups.

    Example: sisana compare survival params.yml

    - Input files:
      - Either a comma-separated or tab-separated metadata file mapping sample names to sample groups.
        Must also contain a column that contains the status of the individual. Must contain True/False 
        values only, where True = dead (event occurred) and False = alive.
    
    - Output files:
      - *survival_plot.png: A Kaplan-Meier plot with the two-sided p-value calculated with sksurv.compare.compare_survival
  
"""




gsea_desc = """
-----------------------------------
Step 4 - GSEA
-----------------------------------

  Perform gene set enrichment analysis (GSEA) to find enriched pathways between two groups

    Example: sisana gsea params.yml

    - Input files:
      - genefile: File (.rnk format, which is two column, tab delimited, no header) containing the 
        genes and test statistics to do enrichment on. The "sisana compare means" option creates
        a .rnk file that can be used for this, or users can supply an alternative file, such as one
        made using limma
      - gmtfile: Gene set file in gmt format
    
    - Output files:
      - Prerank_GSEA_*.txt: Tab-separated file containing the output of GSEA
      - GSEA_*_basic_enrichment_plot.png: GSEA basic enrichment plot, displaying the top 5 enriched terms
      - GSEA_*_basic_enrichment_dotplot.png: GSEA basic enrichment dotplot, displaying the top 10 enriched terms
      - Plus some additional files automatically created with the gseapy module

"""
      
visualize_desc = """
-----------------------------------
Step 4 - Visualize
-----------------------------------

Visualize is the fourth step in the SiSaNA pipeline. In this step, a variety of visualizations for the 
expression data and resulting networks can be created.

  ----- Option 1: volcano -----
  Creates a volcano pot of the input expression data
  
  Example: sisana visualize volcano params.yml
  
    - Input files:
      - datafile: TSV file that contains the fold change, p-value, FDR, and mean expression for each gene
    
    - Output files:
      - volcano_plot_*.png: A PNG volcano plot, where X-axis is the log2 fold change and y-axis is the 
        Benjamini-Hochberg adjusted p-value. 
        
  ----- Option 2: quantity -----
  Creates a boxplot or violin plot of the degrees/expression per-group
  
  Example: sisana visualize quantity params.yml
  
    - Input files:
      - datafile: Either a TSV or CSV file that contains the expression or indegrees of each gene per sample.
      - metadata: CSV metadata file mapping samples to groups (groups must match names of the groupnames arg), 
        must have a header of the format 'name,group'
      - genelist: .txt file containing simply a list of genes to plot, must match the name of genes in the datafile. 
        Recommended to not use more than 10 genes, otherwise visualize with a heatmap.    
   
    - Output files:
      - *violin_plot.png or *boxplot.png: The resulting violin plot or box plot in PNG format.

  ----- Option 3: heatmap -----
  Creates a heatmap on a per-group basis of the user-specified genes
  
  Example: sisana visualize heatmap params.yml
  
    - Input files:
      - datafile: Either a TSV or CSV file that contains the expression or indegrees of each gene per sample.
      - metadata: CSV metadata file mapping samples to groups (groups must match names of the groupnames arg), 
        must have a header of the format 'name,group'
      - genelist: .txt file containing simply a list of genes to plot, must match the name of genes in the datafile. 
        Recommended to not use more than 10 genes, otherwise use a heatmap.    
   
    - Output files:
      - *_heatmap.png: The heatmap itself, where rows are genes/TFs, columns are samples, and then the columns are
        grouped based on their input sample groups
        
"""

extract_desc = """
-----------------------------------
(Optional step) - Extract
-----------------------------------

This is an optional step that allows the user to extract edges connected to specific TFs or genes.
Note that after "sisana extract", you also need to supply either "genes" or "tfs" to tell SiSaNA
what type of data you are inputting.

Examples: sisana extract genes params.yml

          OR
          
          sisana extract tfs params.yml

  - Input files:
    - pickle: The lioness.pickle file created in the "generate" step
    - sampnames: The path to the tmp/samples.txt file, which is a text file that is generated in the preprocess step
    - symbols: The path to the file that contains the names of TFs or genes you wish to extract
  
  - Output files:
    - lioness_filtered_for_*.csv: Lioness networks filtered for the user-defined TFs/genes
    - lioness_available_*.csv: List of the names of genes or TSVs available to filter for. 
        - Note: This file is only created if an error occurs when trying to filter for genes/TFs.
      
"""