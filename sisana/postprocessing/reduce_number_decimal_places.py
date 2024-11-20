
import numpy as np
import pandas as pd
import argparse
import sys
import pickle
from pathlib import Path
from post import save_results
import time
import os 


if __name__ == '__main__':
    """
    Description:
        The python code for run_lioness.py only allows for either 7 or 15 decimal places, which creates large files. This code takes the output of lioness and reduces the number of decimal places saved in the file.
        The previous lioness file will then need to be deleted manually if the user wishes to clear up space.
    """
    
    parser = argparse.ArgumentParser(description="Example command: python reduce_number_decimal_places.py -n lioness_df.pickle> -i pickle -o ./output/ -f csv -d 3")
    ArgGroup = parser.add_argument_group('Required arguments')  
    
    ArgGroup.add_argument("-n", "--filename", type=str, help="Path to file in csv or pickle format", required=True)
    ArgGroup.add_argument("-i", "--informat", type=str, choices = ["csv", "pickle"], help="Format of the input file", required=True)    
    ArgGroup.add_argument("-o", "--outdir", type=str, help="Path to directory to output file to", required=True) 
    ArgGroup.add_argument("-f", "--outformat", type=str, choices = ["csv", "pickle", "RData"], help="Format of file to output the network to", required=True) 
    ArgGroup.add_argument("-d", "--decimal", type=str, help="Number of decimal places to keep", required=True) 
    
    start = time.time()    

    args = parser.parse_args()
    
    pd.set_option("display.precision", 3)
    
    if args.informat == "csv":
        nwdf = pd.read_csv(args.filename)
    elif args.informat == "pickle":
        nwdf = pd.read_pickle(args.filename)

    base_file_name = Path(args.filename).stem
    
    print("Converting decimals, please wait...")
    round_nwdf = nwdf.round(decimals=int(args.decimal))
    
    print("Saving results, please wait...")
    
    if args.outformat == "csv" or args.outformat == "pickle":
        save_results(round_nwdf, args.outformat, base_file_name, f"{args.decimal}_decimal_places", args.outdir)
    else:
            
        import rpy2
        from rpy2 import robjects
        from rpy2.robjects import pandas2ri
        import rpy2.robjects as ro
        pandas2ri.activate()

        # write pandas dataframe to an .RData file
        def save_rdata_file(df, filename):
            #r_data = pandas2ri.py2ri(df)

            convert_time_start = time.time()  

            with (ro.default_converter + pandas2ri.converter).context():
                r_from_pd_df = ro.conversion.get_conversion().py2rpy(df)
            
            robjects.r.assign("lioness_df", r_from_pd_df)
            robjects.r("save(lioness_df, file='{}')".format(filename))  
        
            convert_time_end = time.time()
            convert_time = convert_time_end - convert_time_start
            print(f"Time taken for converting to RData file: {convert_time}")
            
        save_rdata_file(round_nwdf, os.path.join(args.outdir, f"{base_file_name}_{args.decimal}_decimal_places.RData"))
    
    end = time.time()  
    totaltime = end - start
    print(f"Time taken for whole script: {totaltime} seconds")
