import numpy as np
import pandas as pd
import os
import argparse
import random

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Example command: python permutation_test.py -f TFs.txt -o ./output/")
    ArgGroup = parser.add_argument_group('Required arguments')  
    ArgGroup.add_argument("-f", "--inputfile", type=str, help="Path to file containing just a list of TFs/genes to select from, no header", required=True) 
    ArgGroup.add_argument("-n", "--number", type=int, help="Number of times to randomly select values", required=True) 
    ArgGroup.add_argument("-o", "--outputdir", type=str, help="Path to ouput folder", required=True)
    
    args = parser.parse_args()

    if args.number < 100:
        raise Exception("Error: Please enter a value of 100 or greater for --number.")  

    # Read in data to a list
    permlist = []
    with open(args.inputfile) as f:
        for line in f:
            permlist.append(line.strip())

    # Remove empty values in list (which stem from empty lines in file)
    permlist = [x for x in permlist if x]

    # Perform random selection from list
    selected_vals = []
    i = 0
    while i < args.number:
        selected_vals.append(random.choice(permlist))
        i += 1

    # Count how many times each value appears
    uniqvals = np.unique(permlist)
    print(uniqvals)

    counts = {}
    for val in uniqvals:
        counts[val] = selected_vals.count(val)

    all_stats = {}
    for k,v in counts.items():
        stats = {}
        stats['Number_occurrences'] = int(v)
        stats['Percent_of_total'] = f"{v/args.number * 100:.2f}%"
        stats['p-value'] = round(v/args.number, 2)
        stats['Occurrences/expected'] = round(v/(args.number/len(uniqvals)), 2)

    #    print(f"{k}: Number occurrences - {v}    Percent of total - {v/args.number * 100}%    p-value - {1 - v/args.number:.2f}") 
    
        all_stats[k] = stats

    stats_table = pd.DataFrame.from_dict(all_stats, orient='columns').T
    print(stats_table)
