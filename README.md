# SiSaNA
Single Sample Network Analysis

SiSaNA is used after creating both Panda and Lioness networks from the package netZooPy. SiSaNA takes the Lioness output, processes it to be analyzed downstream, and then calculates in- and out-degree for each of the reconstructed networks.

### Pre-processing of data
This step is actually performed prior to running Panda/Lioness, and it filters the expression matrix, PPI file, and prior motif to contain the same genes/TFs, which is necessary for running Panda/Lioness.


#### Usage
```
python filter_exp_min_samps.py -e <expression_file.tsv> -m <motif_file.txt> -p <ppi_file.txt> -n 10
```

#### Inputs
 - `-e`: Path to file containing the gene expression data. Row names must be genes, the expression file does not have a header, and the cells are separated by a tab
 - `m`: Path to motif file, which gets filtered to only contain genes that pass the minimum number of samples threshold
 - `p`: Path to ppi file, which gets filtered to only contain genes that pass the minimum number of samples threshold
 - `n`: Minimum number of samples a gene must be expressed in; expression data will be filtered for only genes that pass

#### Outputs
Three files, one for each of the three filtered input files. 


### Run PANDA
This step creates a Panda network from the filtered files. See documentation for netZooPy (https://github.com/netZoo/netZooPy/tree/master). An example command is given below.

#### Usage
```
python run_panda.py -e <expression_data_filtered>.txt -m <motif_data_filtered>.txt -p <ppi_data_filtered>.txt -r True -o <output_file>.txt
```

### Run LIONESS
Similar to the PANDA step, this step creates Lioness networks from the filtered files. See documentation for netZooPy (https://github.com/netZoo/netZooPy/tree/master). An example command is given below.

#### Usage
```
python run_lioness.py -e <expression_data_filtered>.txt -m <motif_data_filtered>.txt -p <ppi_data_filtered>.txt -g cpu -r single -c 4 -o ./output/ -f mat
```

