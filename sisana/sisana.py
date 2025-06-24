import yaml
import argparse
from netZooPy.panda.panda import Panda
from netZooPy.lioness.lioness import Lioness
from sisana.preprocessing import preprocess_data
from sisana.postprocessing import convert_lion_to_pickle, extract_tfs_genes
from sisana.analyze_networks import calculate_panda_degree, calculate_lioness_degree, compare_bw_groups, survival_analysis, perform_gsea, plot_volcano, plot_expression_degree, plot_heatmap, plot_clustermap
from sisana.example_input import find_example_paths, fetch_files
import sisana.docs
from sisana.docs import create_log_file
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

    # options for extract subcommand
    ext.add_argument("extractchoice", type=str, choices = ["genes", "tfs"], help="Do you want to extract specific gene or TF edges?")   
    ext.add_argument("params", type=str, help='Path to yaml file containing the parameters to use')

    # options for compare subcommand
    comp.add_argument("compchoice", type=str, choices = ["means", "survival"], help="The type of comparison to do")   
    comp.add_argument("params", type=str, help='Path to yaml file containing the parameters to use')

    # options for gsea subcommand    
    gsea.add_argument("params", type=str, help='Path to yaml file containing the parameters to use')

    # options for visualize subcommand
    vis.add_argument("plotchoice", type=str, choices = ["all", "quantity", "heatmap", "volcano"], nargs='?', default="all", help="The type of plot to create")   
    vis.add_argument("params", type=str, help='Path to yaml file containing the parameters to use')

    args = parser.parse_args()
      
    # If user wants example files, retrieve them from their installed paths
    if args.example:
        
        print("Downloading example input files from Zenodo. Please wait...")
        fetch_files()
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
        
        # # Save the order of the sample names to their own file, then export the data frame without a header, since that is what is required for CLI version of PANDA
        # expdf = pd.read_csv(preprocess_params['exp_file'], sep='\t', index_col=0)
        # name_list = list(expdf.columns.values)
        
        # with open('./tmp/samples.txt', 'w') as f:
        #     for line in name_list:
        #         f.write(f"{line}\n")
        
        # Remove genes that are not expressed in at least the user-defined minimum ("number")
        f = preprocess_data(preprocess_params['exp_file'], 
                        preprocess_params['number'],
                        preprocess_params['outdir'])    
        
        file_name = [f]
        
        create_log_file("preprocess", 
                        preprocess_params, 
                        file_name)
        
    ########################################################
    # 2) Run PANDA, using the parameters from the yaml file
    ########################################################

    if args.command == 'generate':
        
        generate_params = params['generate']

        if generate_params['method'].lower() == 'panda' or generate_params['method'].lower() == 'lioness':

            # data_paths = yaml.load(open('./tmp/processed_data_paths.yml'), Loader=yaml.FullLoader)
            
            # Create output dir if one does not already exist
            panda_output_location = generate_params['pandafilepath']

            pandapath = Path(panda_output_location)
            if str(pandapath)[-4:] != ".txt":
                raise Exception("Error: Panda output file must have a .txt extension. Please edit your pandafilepath variable in your params file")
            os.makedirs(pandapath.parent, exist_ok=True)
            
            expdf = pd.read_csv(generate_params['exp'], sep='\t', index_col=0)

            panda_obj = Panda(expression_file=generate_params['exp'], 
                motif_file=generate_params['motif'], 
                ppi_file=generate_params['ppi'], 
                save_tmp=False, 
                remove_missing=False, 
                keep_expression_matrix=True, 
                save_memory=False,
                modeProcess="intersection",
                with_header=True)

            panda_obj.save_panda_results(panda_output_location, old_compatible=False)     
            
            print("Now calculating PANDA degrees...")
            calculate_panda_degree(inputfile=panda_output_location)
               
        if generate_params['method'].lower() == 'lioness':
            lioness_output_location = generate_params['lionessfilepath']
            if lioness_output_location[-4:] != ".npy":
                raise Exception("Error: Lioness output file must have a .npy extension. Please edit your lionessfilepath variable in your params file.")
            
            lionesspath = Path(lioness_output_location)

            lion = Lioness(panda_obj, 
                           computing=generate_params['compute'], 
                           precision="double",
                           ncores=generate_params['ncores'], 
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
        
        if generate_params['method'].lower() == 'lioness':        
            print("\nLIONESS network saved to " + lioness_output_location)
            print(f"LIONESS degrees saved to:")
            print(f"{Path(lioness_output_location).parent}/lioness_indegree.csv")
            print(f"{Path(lioness_output_location).parent}/lioness_outdegree.csv")
            
        outfiles = [panda_output_location,
                    f"{str(panda_output_location)[:-4]}_outdegree.csv",
                    f"{str(panda_output_location)[:-4]}_indegree.csv",
                    lion_loc,
                    f"{Path(lioness_output_location).parent}/lioness_indegree.csv",
                    f"{Path(lioness_output_location).parent}/lioness_outdegree.csv"]
            
        create_log_file("generate", 
                        generate_params, 
                        outfiles)
        
    ########################################################
    # 3) Compare between sample groups
    ########################################################

    if args.command == 'compare':
        
        if args.compchoice == "means":     
            compare_means_params = params['compare']['means']
   
            outfiles = compare_bw_groups(datafile=compare_means_params["datafile"], 
                                        mapfile=compare_means_params["mapfile"], 
                                        datatype=compare_means_params["datatype"], 
                                        groups=compare_means_params["groups"],
                                        testtype=compare_means_params["testtype"], 
                                        filetype=compare_means_params["filetype"],
                                        rankby_col=compare_means_params["rankby"],
                                        outdir=compare_means_params["outdir"])
            
            create_log_file("compare_means", 
                compare_means_params, 
                outfiles)
        
        if args.compchoice == "survival":     
            compare_survival_params = params['compare']['survival']

            try:
                outfiles = survival_analysis(metadata=compare_survival_params["metadata"],
                                filetype=compare_survival_params["filetype"], 
                                sampgroup_colname=compare_survival_params["sampgroup_colname"],
                                alivestatus_colname=compare_survival_params["alivestatus_colname"],
                                days_colname=compare_survival_params["days_colname"],
                                groups=compare_survival_params["groups"],
                                outdir=compare_survival_params["outdir"],
                                appendname=compare_survival_params["appendname"])
            except:
                outfiles = survival_analysis(metadata=compare_survival_params["metadata"],
                                filetype=compare_survival_params["filetype"], 
                                sampgroup_colname=compare_survival_params["sampgroup_colname"],
                                alivestatus_colname=compare_survival_params["alivestatus_colname"],
                                days_colname=compare_survival_params["days_colname"],
                                groups=compare_survival_params["groups"],
                                outdir=compare_survival_params["outdir"])
                
            create_log_file("compare_survival", 
                compare_survival_params, 
                [outfiles])

    ########################################################
    # 4) Perform gene set enrichment analysis
    ########################################################   
        
    if args.command == 'gsea':    
        gsea_params = params["gsea"]
        
        outfiles = perform_gsea(genefile=gsea_params["genefile"], 
                        gmtfile=gsea_params["gmtfile"], 
                        geneset=gsea_params["geneset"], 
                        outdir=gsea_params["outdir"])
        
        create_log_file("gsea", 
            gsea_params, 
            outfiles)
    
    ########################################################
    # 5) Visualize results
    ########################################################       

    if args.command == "visualize":                  

        if args.plotchoice == "volcano":    
            volcano_params = params["visualize"]["volcano"]

            outfiles = plot_volcano(statsfile=volcano_params["statsfile"],
                         diffcol=volcano_params["diffcol"],
                         adjpcol=volcano_params["adjpcol"],
                         adjpvalthreshold=volcano_params["adjpvalthreshold"],
                         xaxisthreshold=volcano_params["xaxisthreshold"],
                         difftype=volcano_params["difftype"],
                         genelist=volcano_params["genelist"],
                         outdir=volcano_params["outdir"],
                         top=False)      
            
            create_log_file("volcano_plot", 
                volcano_params, 
                [outfiles])
    
        if args.plotchoice == "quantity":   
            quantity_params = params["visualize"]["quantity"]
            
            if quantity_params["genelist"] != None:
                outfiles = plot_expression_degree(datafile=quantity_params["datafile"],
                            filetype=quantity_params["filetype"], 
                            statsfile=quantity_params["statsfile"], 
                            metadata=quantity_params["metadata"],
                            plottype=quantity_params["plottype"],
                            groups=quantity_params["groups"],
                            colors=quantity_params["colors"],
                            prefix=quantity_params["prefix"],
                            yaxisname=quantity_params["yaxisname"],
                            outdir=quantity_params["outdir"],
                            genelist=quantity_params["genelist"],
                            top=False)   
            else:
                outfiles = plot_expression_degree(datafile=quantity_params["datafile"],
                            filetype=quantity_params["filetype"], 
                            statsfile=quantity_params["statsfile"], 
                            metadata=quantity_params["metadata"],
                            plottype=quantity_params["plottype"],
                            groups=quantity_params["groups"],
                            colors=quantity_params["colors"],
                            prefix=quantity_params["prefix"],
                            yaxisname=quantity_params["yaxisname"],
                            outdir=quantity_params["outdir"],
                            genelist=quantity_params["genelist"],
                            numgenes=quantity_params["numgenes"],
                            top=True)   
                
            create_log_file("quantity_plot", 
                quantity_params, 
                [outfiles])               
                
        # For now, the plot_heatmap option is being deprecated for use of the plot_clustermap option instead,
        # as the clustermap option allows for more user control and clustering of patients/parameters
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
            
        if args.plotchoice == "heatmap":    
            heatmap_params = params["visualize"]["heatmap"]

            outfiles = plot_clustermap(datafile=heatmap_params["datafile"],
                        filetype=heatmap_params["filetype"], 
                        metadata=heatmap_params["metadata"],
                        genelist=heatmap_params["genelist"],
                        column_cluster=heatmap_params["column_cluster"],
                        row_cluster=heatmap_params["row_cluster"],
                        prefix=heatmap_params["prefix"],
                        outdir=heatmap_params["outdir"],
                        plot_gene_names=heatmap_params["plot_gene_names"],
                        plot_sample_names=heatmap_params["plot_sample_names"],
                        category_label_columns=heatmap_params["category_label_columns"],
                        category_column_colors=heatmap_params["category_column_colors"],                       
                        top=False)   
            
            create_log_file("heatmap", 
                heatmap_params, 
                outfiles)  
            
    ########################################################
    # 6) Optional, extract edges that connect to specific TFs/genes
    ########################################################

    if args.command == 'extract':
        extract_params = params["extract"]

        outfiles = extract_tfs_genes(pickle=extract_params["pickle"], 
                         datatype=args.extractchoice, 
                         sampnames=extract_params["sampnames"],
                         symbols=extract_params["symbols"], 
                         outdir=extract_params["outdir"])
        
        create_log_file("extract", 
                extract_params, 
                [outfiles])  
            
