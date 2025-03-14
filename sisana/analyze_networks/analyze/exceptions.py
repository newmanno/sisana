import pandas as pd
import numpy as np

class NotASubsetError(Exception):
    """
    Raise when the supplied list (genes/samples) is not a subset of the items in the data frame
    
    Attributes:
        user_gene_list : List of genes/samples the user inputs
        data_genes : List of genes/samples in the dataset
        dtype : data type that makes up the list, either "genes" or "samples" 
        message -- explanation of the error
    """

    def __init__(self, user_list, data_list, dtype, message="Error: The genes in the supplied list are not a subset of the genes in the dataset."):
        
        self.user_list = user_list
        self.data_list = np.unique(data_list)
        self.dtype = dtype
        
        items_missing = list(set(self.user_list).difference(self.data_list)) # Find genes/samples not in the data provided
        print(items_missing)
        print(f"\nError: The following {self.dtype} in the supplied list are not present in the data provided:")
        [print(i) for i in items_missing]
        print(f"\nYou have provided the following {self.dtype}:")
        [print(i) for i in self.user_list]
        print(f"\nThe following are the {self.dtype} found in the data provided:")
        [print(i) for i in self.data_list]
        print("\n")

        self.message = message 
        super().__init__(self.message)
        
class IncorrectHeaderError(Exception):
    """
    Raise when the user's metadata file header is not in the correct format of "name,group"
    
    Attributes:
        metadf : Data frame of the user's metadata file
        message -- explanation of the error
    """
    def __init__(self, metadf: pd.DataFrame, message="Error: The header of the metadata file must be in the format 'name,group'."):
        self.metadf = metadf
        print(f"\nThe header of your metadata file looks like the following:\n {self.metadf.head()} \n")
        self.message = message 
        super().__init__(self.message)

class WrongAmountOfColorsError(Exception):
    """
    Raise when the user does not supply the correct number of colors for their given sub-categories (e.g. only one color was given for the category "Sex" that contained "male" and "female" for subcategories)

    Attributes:
        category: str, Name of the category 
        unique_categories: list
        unique_subcategories: int, number of given color codes for the given category 
    """
    def __init__(self, category: str, num_unique_subcategories: int, num_unique_color_codes: str, message="Error: You have not supplied the correct number of colors for a category."):
        self.category = category
        self.num_unique_subcategories = num_unique_subcategories
        self.num_unique_color_codes = num_unique_color_codes
        print(f"\nError: The number of colors specified in your category_column_colors parameter in the params.yml file ({self.num_unique_color_codes}) does not match the number of unique subcategories found for the {category} category ({self.num_unique_subcategories}). Please fix and try running this script again.\n")
        self.message = message 
        super().__init__(self.message)