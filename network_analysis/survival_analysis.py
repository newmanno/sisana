import argparse
import os
import pandas as pd
import pickle
import numpy as np
from analyze import file_to_list, map_samples

__author__ = 'Nolan Newman'
__contact__ = 'nolankn@uio.no'
    
if __name__ == '__main__':
    """
    Description:
        This code performs a survival analysis between two user-defined groups and outputs
        both the survival plot and the statistics for the comparison(s)
    """

    parser = argparse.ArgumentParser(description="Example command: python survival_analysis.py -m map.csv -t txt -c subtypes -g group1 group2 -o ./output")
    ArgGroup = parser.add_argument_group('Required arguments') 
    ArgGroup.add_argument("-m", "--metadata", type=str, help="Path to mapping file (csv) that maps samples to groups", required=True) 
    ArgGroup.add_argument("-t", "--filetype", choices = ["csv", "txt"], help="Type of delimiter used for --datafile", required=True)    
    ArgGroup.add_argument("-c", "--colname", type=str, help="Name of column containing sample group names", required=True) 
    ArgGroup.add_argument("-g", "--compgroups", type=str, nargs=2, help="Name of groups in mapping file to compare", required=True) 
    ArgGroup.add_argument("-o", "--outdir", type=str, help="Path to directory to output file to", required=True) 
    
    args = parser.parse_args()

    # Get data and metadata
    if args.filetype == "csv":
        meta = pd.read_csv(args.metadata, engine = "pyarrow", index_col=[0], header = 0)
    elif args.filetype == "txt":
        meta = pd.read_csv(args.metadata, sep='\t', engine = "pyarrow", index_col=[0], header = 0)
        
    print(meta)
    
    # Assign samples from mapping file to groups
    # sampdict = map_samples(meta, groups[0], groups[1])
    
    # Make suvival plot (code adapted from https://scikit-survival.readthedocs.io/en/stable/user_guide/00-introduction.html)   
    import matplotlib.pyplot as plt
    from sksurv.nonparametric import kaplan_meier_estimator

    for treatment_type in (args.compgroups[0], args.compgroups[1]):
        mask_treat = meta["Basal_sub_subtype"] == treatment_type
        #print(mask_treat)
        time_treatment, survival_prob_treatment, conf_int = kaplan_meier_estimator(
            meta["Alive_status"][mask_treat],
            meta["days_to_death_or_follow_up"][mask_treat],
            conf_type="log-log",
        )

        plt.step(time_treatment, survival_prob_treatment, where="post", label=f"{treatment_type} (n = {mask_treat.sum()})")
        plt.fill_between(time_treatment, conf_int[0], conf_int[1], alpha=0.25, step="post")

    from sksurv.compare import compare_survival

    # Take just the survival status and the time columns and create a structured array for the survival analysis
    small_meta = meta[meta["Basal_sub_subtype"].isin(args.compgroups)]
    meta_for_array = small_meta.loc[:,['Alive_status','days_to_death_or_follow_up']]   
    meta_array = meta_for_array.to_records(index=False)
    group_list = np.array(small_meta["Basal_sub_subtype"])

    surv = compare_survival(meta_array, group_list, return_stats = True)
    
    print(f"\nSurvival analysis results:")
    print(f"chi-square test-statistic: {surv[0]:.2E}")
    print(f"p-value: {surv[1]:.2E}")
    
    print(f"\nSurvival analysis statistics:")
    print(surv[2])
     
    print(f"\nSurvival analysis covariance matrix:")
    print(surv[3])

    plt.ylim(0, 1)
    plt.ylabel(r"est. probability of survival $\hat{S}(t)$")
    plt.xlabel("time $t$")
    plt.legend(loc="best")
    plt.text(.7, .7, f"pval: {surv[1]:.2E}", ha='center', va='top', transform=plt.gca().transAxes)
    #plt.text(4, 0.9, f"pval: {surv[1]:.2E}", ha='center')
    
    plt.savefig(args.outdir + "survival_plot.png")
    print(f"\nFile saved: {args.outdir}survival_plot.png")     
            



        