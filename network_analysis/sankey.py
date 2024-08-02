import argparse
import plotly.graph_objects as go
import pandas as pd
import csv
import os 
import matplotlib.pyplot as plt

if __name__ == '__main__':
    """
    Description:
        This code plots the expression or degrees of genes used in the Lioness pipeline
    """

    parser = argparse.ArgumentParser(description="Example command: python heatmap.py -d <indegree.csv> -t csv -f -m <metadata.csv> -g <genelist.txt> -n group1 group2 group3 -o ./output/")
    ArgGroup = parser.add_argument_group('Required arguments')  

    ArgGroup.add_argument("-i", "--input", type=str, help="Path to file containing the two classifications in the first two columns, and the weights in the third column", required=True)
    ArgGroup.add_argument("-o", "--outdir", help="Directory to save resulting plot to", required=True)        

    args = parser.parse_args()

    def create_sankey_dict(mapfile):
        '''
        Function that creates a dictionary of the structure needed for making the Sankey plot

            Arguments:
                - mapfile: pandas dataframe of the mapping file from the user
        '''
        
        with open(mapfile) as samp_file:
            samp_file = csv.reader(samp_file, delimiter = ',')
                        
            source_list = []
            target_list = []
            weights = []
            
            # Populate lists with each column
            for row in samp_file:
                print(row)
            
                source_list.append(row[0])
                target_list.append(row[1])
                weights.append(row[2])
            
            sank = dict(source = source_list, target = target_list, value = weights)

        return(source_list, sank)

    sank_dict = create_sankey_dict(args.input)
    
    for k,v in sank_dict[1].items():
        print(k,v)
    
    # Create figure
    fig = go.Figure(data=[go.Sankey(
        node = dict(
        pad = 15,
        thickness = 20,
        line = dict(color = "black", width = 0.5),
        label = sank_dict[0]
        
        ),
        link = sank_dict[1]
    )])

    fig.update_layout(title_text="Basic Sankey Diagram", font_size=10,width=600, height=400)
    fig.show()

    outname = os.path.join(args.outdir, f"sankey.png")
    fig.write_image(outname)

    print(f"\nFile created: {outname}")
