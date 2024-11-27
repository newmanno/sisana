preprocess_desc = """
-----------------------------------
Step 1 - Preprocess
-----------------------------------

Preprocess is the first step to running SiSaNA. In this step, SiSaNA takes the 
three provided input files given to params.yml (the protein-protein interaction file,
prior motif, and gene expression) and filters them to only contain the proteins/genes
that are present in all files. It also uses the "number" argument of the params.yml
file to filter out genes that are not expressed in at least the specified number of
samples. 

Example: sisana preprocess params.yml

Input files (specified in the params.yml file): 
  - exp_file: A protein-protein interaction file of known PPIs. Three column, tab-separated, 
    no header. The first and seccond columns are TF names and the third column is
    a binary value (0 or 1), indicating whether there is a known interaction between
    the two TFs.
  - motif_file: A tab-separated transcription factor binding motif file (no header) describing 
    whether the TF (column 1) binds to the motif of the gene (column 2). The third 
    column is a binary value (0 or 1) that indicates whether the TF is known to bind 
    to that motif.
  - exp_file: A tab-separated gene expression file with a header, where rows are gene names and
    columns are sample names.

Output files: 
  - *_filtered.txt files: (3 total) in the user-specified output directory. These are the 
    input files, filtered to only contain the proteins/genes across all three. The gene
    expression file also only contains genes that are expressed in the user-defined number of 
    samples. 
  - *filtering_statistics.txt: Contains information regarding the number of genes/proteins filtered
    out of each file
    
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
  - processed_paths: The path to the "processed_data_paths.yml" file created in the "preprocess"
    step. This is saved in the ./tmp/ directory.

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
      - mapfile: CSV mapping file, with a header in the format "name,group", which maps sample name to sample group
      
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
                
  ----- Option 3: gsea -----
  Performs gene set enrichment analysis (GSEA) to find enriched pathways between two groups

    Example: sisana compare gsea params.yml

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

Example: sisana extract genes params.yml

  - Input files:
    - pickle: The lioness.pickle file created in the "generate" step
    - namefile: The path to the file that contains the names of TFs or genes you wish to extract
  
  - Output files:
    - lioness_filtered_for_*.csv: Lioness networks filtered for the user-defined TFs/genes
    - lioness_available_*.csv: List of the names of genes or TSVs available to filter for. 
      Note: This file is only created if an error occurs when trying to filter for genes/TFs.
      
"""