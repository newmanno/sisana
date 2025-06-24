import sys
import os

def create_log_file(subcommand: str, params_dict: dict, filenames: list): 
    """
    Description:
        This function creates log files for each command the user performs
        
    Parameters:
    -----------
        - subcommand: str, The name of the subcommand used
        - params_dict: dict, The dictionary of the parameters the user supplied
        - filenames: list, The list of file paths that were created with the subcommand
        
    Returns:
    -----------
        - str of the output file location 
    """
    os.makedirs("./log_files/", exist_ok=True)
    
    basename = f"{subcommand}_log.txt"
    file_outloc = os.path.join("./log_files/", basename)
    
    # Remove the "./" prefix to file names if found. This is for sake of 
    # clarity, since otherwise some file names had them and some did not,
    # just depending on how the user defined them in params file
    fixed_names = [n[2:] if n[:2] == "./" else n for n in filenames]
    
    with open(file_outloc, "w") as file:
        file.write("Directory analysis was performed in:\n")
        file.write(f"   - {os.getcwd()}\n")
        
        file.write("\nParameters used:\n")
        for k,v in params_dict.items():
            file.write(f"   - {k}: {v}\n")
        
        file.write("\nFile(s) generated:\n")
        for i in fixed_names: 
            file.write("   - " + i + "\n")

            
