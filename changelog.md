### 1.XX
- Changed the listed files in the log output files to all not have "./" at the beginning of their paths, just for consistency's sake 
- Added this changelog file

### 1.4.0
This new version of SiSaNA now generates log files as well, allowing the ability to find and reference the parameters you used for each analysis performed. The log files are automatically generated from the root project directory into a folder titled log_files. An example of the log file is given below.

### 1.3.0
SiSaNA now has the option to create clustermaps, allowing the user to visualize clusters of samples and genes/TFs. For this option, users also specify categorical metadata columns to color the samples on. 

### 1.2.0
With this newest version of SiSaNA, all analysis is now performed with the use of a params.yml file instead of specifying arguments directly via the command line.

### 1.1.0
Many features are now available in the new version. These include the following:

1. Filtering of all PPI, motif, and gene expression files, which is a prerequisite for running PANDA/LIONESS.
2. Filtering the output of Lioness for only edges found in the prior.
3. Calculation of the in-and out-degree of genes and TFs, respectively.
4. Reducing the number of decimal places in either the PANDA/LIONESS output or the calculated in-/out-degrees, which greatly saves on storage space.
5. Extraction of specific TFs/genes, which is useful for analyses such as limma (see Ritchie et al., 2015). An option to perform analyses with limma is not available as part of this software.
6. Comparison of groups identified in dimension reduction techniques such as UMAP or tSNE. These include comparisons of TFs/genes between two groups, survival analysis between groups, and gene set enrichment analysis (GSEA)
7. Visualization of these results via volcano plots, box plots, violin plots, and heatmaps is also possible.

### 1.0.0
Initial release. SiSaNA can reconstruct networks using PANDA/LIONESS as well as calculate in- and out-degree and compare degrees between groups.
