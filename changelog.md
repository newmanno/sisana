### 1.8.1
- Removed a file that was not used in the pipeline and was meant for testing purposes originally

### 1.8.0
- Fixed issue where the text for the quantity plot in the summarize file did not show the correct number of genes plotted

### 1.7.3
- Fixed issue where the header validation did not catch files without headers when using a metadata file in the format for a paired comparison
- Changed gsea color bar formatting so the title does not overlap the dot plot
 
### 1.7.2
- Added additional input validation for checking whether the user's file name matched the input file extension
- Added the requirement for user to specify metadata file type, as previously the user was forced to use a csv file in some steps
- Fixed problem with paired comparisons where code would not execute
- Added unit tests for all comparison types and reformatted so tt, mw, and paired_tt are now all modular functions that are called

### 1.7.1
- Pinned the versions of required modules to requirements.txt

### 1.7.0
- Changed the name of the "generate" command to "reconstruct", with backwards compatibility for generate still
- Fixed an issue where the error "ValueError: setting an array element with a sequence." would reportedly occur for some users if using newer versions of pandas/NumPy
- Added a "data_color" parameter for the heatmap so the user can choose the color map to use for the heatmap (based on the given color map from https://matplotlib.org/stable/users/explain/colors/colormaps.html)
- Added a "color" parameter for gsea that will now color the resulting dotplot based on the given color map from https://matplotlib.org/stable/users/explain/colors/colormaps.html
- Added a "colors" parameter for the volcano plot to specify the colors used in the plot
- Added a "colors" parameter for the survival plot to specify the colors for the lines
- Visualize and survival commands now output a pickled file containing the plot(s) created in that step, so advanced users can import them into python and modify them
- Made an exceptions directory in the root sisana folder for all custom exceptions, moved the previous custom exceptions from the analyze subdirectory to this new directory

### 1.6.1
- Fixed an issue in the previous release where the older example Zenodo files were installed. The visualization params in this file were incompatible with the new version.
- Set a required GSEApy version, as newer updates to GSEApy were causing different plot colors and labels
- Added version information for other packages in pyproject as well

### 1.6.0
- (Note: skipped 1.5.3 GitHub release to stay consistent with PyPI)
- Corrected the bug introduced in v1.5.0 where the panda output was sorted, then used as input to an unsorted df later on, resulting in incorrect results
- Fixed a bug where the visualization commands would fail due to incorrect arguments being supllied
- Added input validation for the params file
- Added version information to the log files
- Added unit tests
- Labeled required and optional parameters in the params.yml file. Some parameters, like the "genelist" parameters in the visualization parameters are required to exist in the params.yml file, but must be left blank if the user does not want to use them (e.g. it will look like "genelist: " instead of "genelist: file.txt")

### 1.5.2 [Note: BAD RELEASE, DO NOT USE]
- Fixed a bug where the incorrect genes are labeled on a volcano plot
- Fixed a bug where SiSaNA required the "preprocess" parameter to be defined in the params.yml file, even if the user was not running the preprocess step.
- Added functionality for the modeProcess argument from netZooPy
- Added functionality for specifiying whether to run PANDA on CPU or GPU
- Added a hotstart option for the lioness_to_pickle_df.py step, so that if the user has an issue in the "sisana generate" step after network reconstruction, they can jump back into the analysis without needing to re-reconstruct all the networks

### 1.5.1 [Note: BAD RELEASE, DO NOT USE]
- Added the example Jupyter notebook for clustering
- Fixed the spacing issue and the text size of the legend in the heatmap step. Note that this fix only works when using a single metadata group (e.g. subtype). Two or more metadata types are still not supported at this time.
- Added the example output results_summarized.html file to the GitHub repo.
- Fixed the example file downloads so it downloads from the Kuijjer Lab repo 

### 1.5.0 [Note: BAD RELEASE, DO NOT USE]
- Added option to run LIONESS on a subset of samples for batching of the "generate" step. PANDA is still ran on all samples, so the same background is used for each batch. 
- Removed some hard-coded values in volcano plot, so now all user-defined values should work properly
- Changed the labels on the volcano plot to make it more clear which groups were which
- When visualizing heatmaps previously, an error would occur if the metadata header was not in 
the correct format. This has now been fixed.
- The "sisana compare means" command has been changed to just "sisana compare". Likewise, the "sisana compare survival" command is now just "sisana survival"
- Added a "sisana summarize" command, which summarizes the analysis the user performs in an html format.

### 1.4.1
- Changed the listed files in the log output files to all not have "./" at the beginning of their paths, just for consistency's sake
- Added more descriptive descriptions of each parameter in the params.yml file
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
