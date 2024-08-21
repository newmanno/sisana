class NotASubsetError(Exception):
    """
    Raise when the supplied gene list is not a subset of the genes in the dataset
    
    Attributes:
        user_gene_list -- List of genes the user wants to plot
        data_genes -- List of genes in the dataset
        message -- explanation of the error
    """

    def __init__(self, user_gene_list, data_genes, message=f"Error: The genes in the supplied list are not a subset of the genes in the dataset."):
        
        self.user_gene_list = user_gene_list
        self.data_genes = data_genes
        
        genes_missing = list(set(self.user_gene_list).difference(self.data_genes)) # Find genes not in the data provided
        print(genes_missing)
        print("\nError: The following genes in the supplied gene list are not present in the data provided:")
        [print(i) for i in genes_missing]
        print("\n")

        self.message = message 
        super().__init__(self.message)
        
        
        
        
        