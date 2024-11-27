import argparse
import os
import pandas as pd
import pickle
import numpy as np
from .analyze import file_to_list, map_samples
from sksurv.compare import compare_survival
from pathlib import Path


__author__ = 'Nolan Newman'
__contact__ = 'nolankn@uio.no'
    
def survival_analysis(metadata: str, filetype: str, sampgroup_colname: str, alivestatus_colname: str, days_colname: str, groups: list, outdir: str):
    """
    Description:
        This code performs a survival analysis between two user-defined groups and outputs
        both the survival plot and the statistics for the comparison(s)
        
    Parameters:
    -----------
        - metadata: str, Path to mapping file ("csv" or "txt") that maps samples to groups
        - filetype: str, Type of delimiter used for metadata
        - sampgroup_colname: str, Name of column containing sample group names
        - alivestatus_colname: str, Name of column that contains the status of the individual. Must contain True/False values only, where True = dead (event occurred) and False = alive.
        - days_colname: str, Name of column containing either the number of days an individual survived or the number of days to the last follow up.
        - groups: str, The names of the two groups (from the metadata file) to compare
        - outdir: str, The directory to save the output to
        
    Returns:
    -----------
        - Nothing
    """

    # parser = argparse.ArgumentParser(description="Example command: python survival_analysis.py -m map.csv -t txt -c subtypes -g group1 group2 -o ./output")
    # ArgGroup = parser.add_argument_group('Required arguments') 
    # ArgGroup.add_argument("-m", "--metadata", type=str, help="Path to mapping file (csv) that maps samples to groups", required=True) 
    # ArgGroup.add_argument("-t", "--filetype", choices = ["csv", "txt"], help="Type of delimiter used for --datafile", required=True)    
    # ArgGroup.add_argument("-c", "--colname", type=str, help="Name of column containing sample group names", required=True) 
    # ArgGroup.add_argument("-a", "--alivestatus", type=str, help="Name of column that contains the status of the individual. Must contain True/False values only, where True = dead (event occurred) and False = alive.", required=True) 
    # ArgGroup.add_argument("-d", "--days", type=str, help="Name of column containing either the number of days an individual survived or the number of days to the last follow up.", required=True) 
    # ArgGroup.add_argument("-g", "--compgroups", type=str, nargs=2, help="Name of groups in mapping file to compare", required=True) 
    # ArgGroup.add_argument("-o", "--outdir", type=str, help="Path to directory to output file to", required=True) 
    
    # args = parser.parse_args()

    # Create output directory if one does not already exist
    os.makedirs(outdir, exist_ok=True)

    # Get data and metadata
    if filetype == "csv":
        meta = pd.read_csv(metadata, engine = "python", index_col=[0], header = 0)
    elif filetype == "txt":
        meta = pd.read_csv(metadata, sep='\t', engine = "python", index_col=[0], header = 0)
            
    # Make suvival plot (code adapted from https://scikit-survival.readthedocs.io/en/stable/user_guide/00-introduction.html)   
    import matplotlib.pyplot as plt
    from sksurv.nonparametric import kaplan_meier_estimator

    for treatment_type in (groups[0], groups[1]):
        mask_treat = meta[sampgroup_colname] == treatment_type

        time_treatment, survival_prob_treatment, conf_int = kaplan_meier_estimator(
            meta[alivestatus_colname][mask_treat],
            meta[days_colname][mask_treat],
            conf_type="log-log",
        )

        plt.step(time_treatment, survival_prob_treatment, where="post", label=f"{treatment_type} (n = {mask_treat.sum()})")
        plt.fill_between(time_treatment, conf_int[0], conf_int[1], alpha=0.25, step="post")

    # Take just the survival status and the time columns and create a structured array for the survival analysis
    small_meta = meta[meta[sampgroup_colname].isin(groups)]
    meta_for_array = small_meta.loc[:,[alivestatus_colname, days_colname]]   
    meta_array = meta_for_array.to_records(index=False)
    group_list = np.array(small_meta[alivestatus_colname])

    print(meta_array)

    surv = compare_survival(meta_array, group_list, return_stats = True)
  
    plt.ylim(0, 1)
    plt.ylabel(r"est. probability of survival $\hat{S}(t)$")
    plt.xlabel("time $t$")
    plt.legend(loc="best")
    plt.text(.7, .7, f"pval: {surv[1]:.2E}", ha='center', va='top', transform=plt.gca().transAxes)
    
    print(f"\nSurvival analysis results:")
    print(f"chi-square test-statistic: {surv[0]:.2E}")
    print(f"p-value: {surv[1]:.2E}")
    
    print(f"\nSurvival analysis statistics:")
    print(surv[2])
     
    print(f"\nSurvival analysis covariance matrix:")
    print(surv[3])

    outdir_path = Path(outdir)
    plt.savefig(f"{outdir_path}/{groups[0]}_v_{groups[1]}_survival_plot.png")
    print(f"\nFile saved: {outdir_path}/{groups[0]}_v_{groups[1]}_survival_plot.png")     
  


        
