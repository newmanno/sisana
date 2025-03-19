import yaml
import argparse
from netZooPy.panda.panda import Panda
from netZooPy.lioness.lioness import Lioness
from sisana.preprocessing import preprocess_data
from sisana.postprocessing import convert_lion_to_pickle, extract_tfs_genes
from sisana.analyze_networks import calculate_panda_degree, calculate_lioness_degree, compare_bw_groups, survival_analysis, perform_gsea, plot_volcano, plot_expression_degree, plot_heatmap, plot_clustermap
from sisana.example_input import find_example_paths
import sisana.docs
import os 
import pandas as pd
import sys
import numpy as np
from pathlib import Path

def cli():
    """
    SiSaNA command line interface
    """

    DESCRIPTION = """
    SiSaNA - Single Sample Network Analysis
    A command line interface tool used to generate and analyze 
    PANDA and LIONESS networks. It works through subcommands. 
    The command 'sisana generate -p params.yaml', for example,
    will reconstruct a PANDA or LIONESS network, using the parameters 
    set in the params.yaml file.
    Developed by Nolan Newman (nolan.newman@ncmm.uio.no).
    """
    EPILOG = """
    Code available under MIT license:
    https://github.com/kuijjerlab/sisana
    """
    
    # create the top-level parser
    parser = argparse.ArgumentParser(prog='sisana.py', description=DESCRIPTION, epilog=EPILOG)    
    parser.add_argument('-e', '--example', action='store_true', help='Flag; Copies the example input files into a directory called "./example_inputs"')    
    parser.add_argument('-s', '--setAndForget', action='store_true', help='Flag; Will attempt to run ALL STEPS of SiSaNA at once. Warning: This requires a very well-formatted params file and should not be used by first-time users. Most users will want to run each of the steps individually."')    

    # Add subcommands
    subparsers = parser.add_subparsers(title='Subcommands', dest='command')
    pre = subparsers.add_parser('preprocess', help='Filters expression data for parameters (e.g. genes) that are only present in at least m samples. Also filters each input file so they have the same genes and TFs across each', epilog=sisana.docs.preprocess_desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    gen = subparsers.add_parser('generate', help='Generates PANDA and LIONESS networks', epilog=sisana.docs.generate_desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    ext = subparsers.add_parser('extract', help='Extract edges connected to specified TFs/genes', epilog=sisana.docs.extract_desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    comp = subparsers.add_parser('compare', help='Compare networks between sample groups', epilog=sisana.docs.compare_desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    gsea = subparsers.add_parser('gsea', help='Perform gene set enrichment analysis between sample groups', epilog=sisana.docs.gsea_desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    vis = subparsers.add_parser('visualize', help='Visualize the calculated degrees of each sample group', epilog=sisana.docs.visualize_desc, formatter_class=argparse.RawDescriptionHelpFormatter)

    # options for preprocess subcommand
    # pre.add_argument("-t", "--template", action="store_true", help='Flag for whether to show the path to the template file')
    pre.add_argument("params", type=str, help='Path to yaml file containing the parameters to use')
        
    # options for generate subcommand    
    gen.add_argument("params", type=str, help='Path to yaml file containing the parameters to use')

    # options for compare subcommand
    ext.add_argument("extractchoice", type=str, choices = ["genes", "tfs"], help="Do you want to extract specific gene or TF edges?")   
    ext.add_argument("params", type=str, help='Path to yaml file containing the parameters to use')

    # options for compare subcommand
    comp.add_argument("compchoice", type=str, choices = ["means", "survival"], help="The type of comparison to do")   
    comp.add_argument("params", type=str, help='Path to yaml file containing the parameters to use')

    # options for gsea subcommand    
    gsea.add_argument("params", type=str, help='Path to yaml file containing the parameters to use')

    # options for visualize subcommand
    vis.add_argument("plotchoice", type=str, choices = ["all", "quantity", "heatmap", "volcano", "clustermap"], nargs='?', default="all", help="The type of plot to create")   
    vis.add_argument("params", type=str, help='Path to yaml file containing the parameters to use')

    args = parser.parse_args()
      
    # If user wants example files, retrieve them from their installed paths
    if args.example:
        print("Copying example files. Please wait...")

        import shutil
        import glob
        os.makedirs('./example_inputs/', exist_ok=True)
        all_ex_files = find_example_paths()

        for fname in all_ex_files:
            shutil.copy2(fname, './example_inputs')
        print("Example input files have been created in ./example_inputs/")
        sys.exit(0)
        
    params = yaml.load(open(args.params), Loader=yaml.FullLoader)

    # Create output for temp files if one does not already exist
    os.makedirs('./tmp/', exist_ok=True)

    ########################################################
    # 1) Preprocess the data
    ########################################################

    preprocess_params = params['preprocess']
    if args.command == 'preprocess':
        
        # Save the order of the sample names to their own file, then export the data frame without a header, since that is what is required for CLI version of PANDA
        expdf = pd.read_csv(preprocess_params['exp_file'], sep='\t', index_col=0)
        name_list = list(expdf.columns.values)
        
        with open('./tmp/samples.txt', 'w') as f:
            for line in name_list:
                f.write(f"{line}\n")
        
        expdf.to_csv('./tmp/exp_no_header.csv', header=False, sep="\t")
        
        preprocess_data('./tmp/exp_no_header.csv', 
                        preprocess_params['motif_file'],
                        preprocess_params['ppi_file'],
                        preprocess_params['number'],
                        preprocess_params['outdir'])    
        
    ########################################################
    # 2) Run PANDA, using the parameters from the yaml file
    ########################################################

    elif args.command == 'generate':
        
        if params['generate']['method'].lower() == 'panda' or params['generate']['method'].lower() == 'lioness':

            # gen_params = yaml.load(open(args.params), Loader=yaml.FullLoader)
            # data_paths = yaml.load(open(params['generate']['processed_paths']), Loader=yaml.FullLoader)    
            data_paths = yaml.load(open('./tmp/processed_data_paths.yml'), Loader=yaml.FullLoader)
            
            # Create output dir if one does not already exist
            panda_output_location = params['generate']['pandafilepath']

            pandapath = Path(panda_output_location)
            if str(pandapath)[-4:] != ".txt":
                print(str(pandapath)[-4:])
                raise Exception("Error: Panda output file must have a .txt extension. Please edit your pandafilepath variable in your params file")
            os.makedirs(pandapath.parent, exist_ok=True)

            panda_obj = Panda(expression_file=data_paths['exp'], 
                            motif_file=data_paths['motif'], 
                            ppi_file=data_paths['ppi'], 
                            save_tmp=False, 
                            remove_missing=False, 
                            keep_expression_matrix=True, 
                            save_memory=False)

            panda_obj.save_panda_results(panda_output_location)     
            
            print("Now calculating PANDA degrees...")
            calculate_panda_degree(inputfile=panda_output_location)
               
        if params['generate']['method'].lower() == 'lioness':
            lioness_output_location = params['generate']['lionessfilepath']
            if lioness_output_location[-4:] != ".npy":
                raise Exception("Error: Lioness output file must have a .npy extension. Please edit your lionessfilepath variable in your params file.")
            
            lionesspath = Path(lioness_output_location)

            lion = Lioness(panda_obj, 
                           computing=params['generate']['compute'], 
                           precision="double",
                           ncores=params['generate']['ncores'], 
                           save_dir=lionesspath.parent, 
                           save_fmt="npy")
            
            # Rename the default name of the lioness output file, which is not an option of the current Lioness NetZooPy cli
            os.rename(os.path.join(lionesspath.parent, "lioness.npy"), lioness_output_location)

            #lion_loc = params['generate']['outdir'] + "lioness.npy"
            lion_loc = lioness_output_location
            liondf = pd.DataFrame(np.load(lion_loc))            
                
            # To make the edges positive values for log2FC calculation later on, first need to transform 
            # edges by doing ln(e^w + 1), then calculate degrees. Then you can do the log2FC of degrees
            # in next step
            # 
            # This transformation is described in the paper "Regulatory Network of PD1 Signaling Is Associated 
            # with Prognosis in Glioblastoma Multiforme"
            # print("Now transforming edges...")

            # print("Datafile before transformation")
            # print(liondf.head(n=20))
            
            # lion_transformed = liondf.apply(np.vectorize(transform_edge_to_positive_val))
            
            # print("Datafile after transformation")
            # print(lion_transformed.head(n=20))        
            
            # print("LIONESS network with transformed edge values saved to " + os.path.join(params['generate']['outdir'], "lioness_transformed_edges.npy"))
                        
            pickle_path = './tmp/lioness.pickle'
            convert_lion_to_pickle(panda_output_location,
                                liondf,
                                "npy", 
                                './tmp/samples.txt',  
                                pickle_path)
            
            print("\nLIONESS networks created. Now calculating LIONESS degrees...")
            calculate_lioness_degree(inputfile=pickle_path,
                            datatype="pickle")
            print("LIONESS degrees have now been calculated.")

            # Move degree files from .tmp to user's output location
            Path("./tmp/lioness_indegree.csv").rename(f"{Path(lioness_output_location).parent}/lioness_indegree.csv")
            Path("./tmp/lioness_outdegree.csv").rename(f"{Path(lioness_output_location).parent}/lioness_outdegree.csv")
                
        print(f"\nPANDA network saved to {panda_output_location}")
        print(f"PANDA degrees saved to:") 
        print(f"{str(panda_output_location)[:-4]}_outdegree.csv")
        print(f"{str(panda_output_location)[:-4]}_indegree.csv")
        
        if params['generate']['method'].lower() == 'lioness':        
            print("\nLIONESS network saved to " + lioness_output_location)
            print(f"LIONESS degrees saved to:")
            print(f"{os.path.splitext(lioness_output_location)[0]}_outdegree.csv")
            print(f"{os.path.splitext(lioness_output_location)[0]}_indegree.csv")

    ########################################################
    # 3) Compare between sample groups
    ########################################################

    elif args.command == 'compare':
        
        if args.compchoice == "means":        
            compare_bw_groups(datafile=params["compare"]["means"]["datafile"], 
                              mapfile=params["compare"]["means"]["mapfile"], 
                              datatype=params["compare"]["means"]["datatype"], 
                              groups=params["compare"]["means"]["groups"],
                              testtype=params["compare"]["means"]["testtype"], 
                              filetype=params["compare"]["means"]["filetype"], 
                              outdir=params["compare"]["means"]["outdir"])
        
        if args.compchoice == "survival":     
            try:
                survival_analysis(metadata=params["compare"]["survival"]["metadata"],
                                filetype=params["compare"]["survival"]["filetype"], 
                                sampgroup_colname=params["compare"]["survival"]["sampgroup_colname"],
                                alivestatus_colname=params["compare"]["survival"]["alivestatus_colname"],
                                days_colname=params["compare"]["survival"]["days_colname"],
                                groups=params["compare"]["survival"]["groups"],
                                outdir=params["compare"]["survival"]["outdir"],
                                appendname=params["compare"]["survival"]["appendname"])
            except:
                survival_analysis(metadata=params["compare"]["survival"]["metadata"],
                                filetype=params["compare"]["survival"]["filetype"], 
                                sampgroup_colname=params["compare"]["survival"]["sampgroup_colname"],
                                alivestatus_colname=params["compare"]["survival"]["alivestatus_colname"],
                                days_colname=params["compare"]["survival"]["days_colname"],
                                groups=params["compare"]["survival"]["groups"],
                                outdir=params["compare"]["survival"]["outdir"])

    ########################################################
    # 4) Perform gene set enrichment analysis
    ########################################################   
        
    elif args.command == 'gsea':    
        perform_gsea(genefile=params["gsea"]["genefile"], 
                        gmtfile=params["gsea"]["gmtfile"], 
                        geneset=params["gsea"]["geneset"], 
                        outdir=params["gsea"]["outdir"])
    
    ########################################################
    # 5) Visualize results
    ########################################################       

    elif args.command == "visualize":

        # if args.plotchoice == "all":
            # # plot top genes in each visualization mode
            # plot_volcano(statsfile=params["visualize"]["volcano"]["statsfile"],
            #              diffcol=params["visualize"]["volcano"]["diffcol"],
            #              adjpcol=params["visualize"]["volcano"]["adjpcol"],
            #              adjpvalthreshold=params["visualize"]["volcano"]["adjpvalthreshold"],
            #              outdir=params["visualize"]["volcano"]["outdir"],
            #              top=True,
            #              numlabels=25)   

            # plot_expression_degree(datafile=params["visualize"]["quantity"]["datafile"],
            #             filetype=params["visualize"]["quantity"]["filetype"], 
            #             statsfile=params["visualize"]["quantity"]["statsfile"],                         
            #             metadata=params["visualize"]["quantity"]["metadata"],
            #             genelist=params["visualize"]["quantity"]["genelist"],
            #             plottype=params["visualize"]["quantity"]["plottype"],
            #             groups=params["visualize"]["quantity"]["groups"],
            #             colors=params["visualize"]["quantity"]["colors"],
            #             prefix=params["visualize"]["quantity"]["prefix"],
            #             yaxisname=params["visualize"]["quantity"]["yaxisname"],
            #             outdir=params["visualize"]["quantity"]["outdir"],
            #             top=True)   
            
            # plot_heatmap(datafile=params["visualize"]["heatmap"]["datafile"],
            #             filetype=params["visualize"]["heatmap"]["filetype"], 
            #             statsfile=params["visualize"]["heatmap"]["statsfile"],
            #             metadata=params["visualize"]["heatmap"]["metadata"],
            #             genelist=params["visualize"]["heatmap"]["genelist"],
            #             hierarchicalcluster=params["visualize"]["heatmap"]["hierarchicalcluster"],
            #             groups=params["visualize"]["heatmap"]["groups"],
            #             prefix=params["visualize"]["heatmap"]["prefix"],
            #             plotnames=params["visualize"]["heatmap"]["plotnames"],
            #             outdir=params["visualize"]["heatmap"]["outdir"],
            #             top=True)                           

        if args.plotchoice == "volcano":    
            plot_volcano(statsfile=params["visualize"]["volcano"]["statsfile"],
                         diffcol=params["visualize"]["volcano"]["diffcol"],
                         adjpcol=params["visualize"]["volcano"]["adjpcol"],
                         adjpvalthreshold=params["visualize"]["volcano"]["adjpvalthreshold"],
                         xaxisthreshold=params["visualize"]["volcano"]["xaxisthreshold"],
                         difftype=params["visualize"]["volcano"]["difftype"],
                         genelist=params["visualize"]["volcano"]["genelist"],
                         outdir=params["visualize"]["volcano"]["outdir"],
                         top=False)      
    
        if args.plotchoice == "quantity":   
            if params["visualize"]["quantity"]["genelist"] != None:
                plot_expression_degree(datafile=params["visualize"]["quantity"]["datafile"],
                            filetype=params["visualize"]["quantity"]["filetype"], 
                            statsfile=params["visualize"]["quantity"]["statsfile"], 
                            metadata=params["visualize"]["quantity"]["metadata"],
                            plottype=params["visualize"]["quantity"]["plottype"],
                            groups=params["visualize"]["quantity"]["groups"],
                            colors=params["visualize"]["quantity"]["colors"],
                            prefix=params["visualize"]["quantity"]["prefix"],
                            yaxisname=params["visualize"]["quantity"]["yaxisname"],
                            outdir=params["visualize"]["quantity"]["outdir"],
                            genelist=params["visualize"]["quantity"]["genelist"],
                            top=False)   
            else:
                plot_expression_degree(datafile=params["visualize"]["quantity"]["datafile"],
                            filetype=params["visualize"]["quantity"]["filetype"], 
                            statsfile=params["visualize"]["quantity"]["statsfile"], 
                            metadata=params["visualize"]["quantity"]["metadata"],
                            plottype=params["visualize"]["quantity"]["plottype"],
                            groups=params["visualize"]["quantity"]["groups"],
                            colors=params["visualize"]["quantity"]["colors"],
                            prefix=params["visualize"]["quantity"]["prefix"],
                            yaxisname=params["visualize"]["quantity"]["yaxisname"],
                            outdir=params["visualize"]["quantity"]["outdir"],
                            genelist=params["visualize"]["quantity"]["genelist"],
                            numgenes=params["visualize"]["quantity"]["numgenes"],
                            top=True)                   
                
        # if args.plotchoice == "heatmap":    
        #     plot_heatmap(datafile=params["visualize"]["heatmap"]["datafile"],
        #                 filetype=params["visualize"]["heatmap"]["filetype"], 
        #                 statsfile=params["visualize"]["heatmap"]["statsfile"],
        #                 metadata=params["visualize"]["heatmap"]["metadata"],
        #                 genelist=params["visualize"]["heatmap"]["genelist"],
        #                 groups=params["visualize"]["heatmap"]["groups"],
        #                 prefix=params["visualize"]["heatmap"]["prefix"],
        #                 plotnames=params["visualize"]["heatmap"]["plotnames"],
        #                 outdir=params["visualize"]["heatmap"]["outdir"],
        #                 top=False)  
            
        if args.plotchoice == "clustermap":    
            plot_clustermap(datafile=params["visualize"]["clustermap"]["datafile"],
                        filetype=params["visualize"]["clustermap"]["filetype"], 
                        metadata=params["visualize"]["clustermap"]["metadata"],
                        genelist=params["visualize"]["clustermap"]["genelist"],
                        column_cluster=params["visualize"]["clustermap"]["column_cluster"],
                        row_cluster=params["visualize"]["clustermap"]["row_cluster"],
                        prefix=params["visualize"]["clustermap"]["prefix"],
                        outdir=params["visualize"]["clustermap"]["outdir"],
                        plot_gene_names=params["visualize"]["clustermap"]["plot_gene_names"],
                        plot_sample_names=params["visualize"]["clustermap"]["plot_sample_names"],
                        category_label_columns=params["visualize"]["clustermap"]["category_label_columns"],
                        category_column_colors=params["visualize"]["clustermap"]["category_column_colors"],                       
                        top=False)   
            
    ########################################################
    # 6) Optional, extract edges that connect to specific TFs/genes
    ########################################################

    elif args.command == 'extract':
        
        extract_tfs_genes(pickle=params["extract"]["pickle"], 
                         datatype=args.extractchoice, 
                         namefile=params["extract"]["namefile"], 
                         outdir=params["extract"]["outdir"])
            
