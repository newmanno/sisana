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
    
    with open(file_outloc, "w") as file:
        file.write("Directory analysis was performed in:\n")
        file.write(f"   - {os.getcwd()}\n")
        
        file.write("\nParameters used:\n")
        for k,v in params_dict.items():
            file.write(f"   - {k}: {v}\n")
        
        file.write("\nFile(s) generated:\n")
        for i in filenames:
            file.write("   - " + i + "\n")

            
