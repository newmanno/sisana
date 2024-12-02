import yaml
import argparse
from netZooPy.panda.panda import Panda
from netZooPy.lioness.lioness import Lioness
from sisana.preprocessing import preprocess_data
from sisana.postprocessing import convert_lion_to_pickle, extract_tfs_genes
from sisana.analyze_networks import calculate_degree, compare_bw_groups, survival_analysis, perform_gsea, plot_volcano, plot_expression_degree, plot_heatmap
from sisana.example_input import find_ex_paths
import sisana.docs
import os 
import pandas as pd
import sys

def cli():
    """
    SiSaNA command line interface
    """

    DESCRIPTION = """
    SiSaNA - Single Sample Network Analysis
    A sommand line interface tool used to generate and analyze 
    PANDA and LIONESS networks. It works through subcommands. 
    The command 'sisana generate -p params.yaml', for example,
    will reconstruct a panda or lioness network, using the parameters 
    set in the params.yaml file.
    Developed by Nolan Newman (nolan.newman@ncmm.uio.no).
    """
    EPILOG = """
    Code available under GPL-3.0 license:
    https://github.com/newmanno/sisana
    """
    
    # create the top-level parser
    parser = argparse.ArgumentParser(prog='sisana.py', description=DESCRIPTION, epilog=EPILOG)    
    parser.add_argument('-e', '--example', action='store_true', help='Copies the example input files into a directory called "./example_inputs"')    

    # Add subcommands
    subparsers = parser.add_subparsers(title='Subcommands', dest='command')
    pre = subparsers.add_parser('preprocess', help='Filters expression data for parameters (e.g. genes) that are only present in at least m samples. Also filters each input file so they have the same genes and TFs across each', epilog=sisana.docs.preprocess_desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    gen = subparsers.add_parser('generate', help='Generates PANDA and LIONESS networks', epilog=sisana.docs.generate_desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    ext = subparsers.add_parser('extract', help='Extract edges connected to specified TFs/genes', epilog=sisana.docs.extract_desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    comp = subparsers.add_parser('compare', help='Compare networks between sample groups', epilog=sisana.docs.compare_desc, formatter_class=argparse.RawDescriptionHelpFormatter)
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
    comp.add_argument("compchoice", type=str, choices = ["means", "survival", "gsea"], help="The type of comparison to do")   
    comp.add_argument("params", type=str, help='Path to yaml file containing the parameters to use')

    # options for visualize subcommand
    vis.add_argument("plotchoice", type=str, choices = ["volcano", "quantity", "heatmap"], help="The type of plot to create")   
    vis.add_argument("params", type=str, help='Path to yaml file containing the parameters to use')

    args = parser.parse_args()
        
    # If user wants example files, retrieve them from their installed paths
    if args.example:
        print("Copying example files. Please wait...")

        import shutil
        import glob
        os.makedirs('./example_inputs/', exist_ok=True)
        all_ex_files = find_ex_paths()

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
        
        # gen_params = yaml.load(open(args.params), Loader=yaml.FullLoader)
        # data_paths = yaml.load(open(params['generate']['processed_paths']), Loader=yaml.FullLoader)    
        data_paths = yaml.load(open('./tmp/processed_data_paths.yml'), Loader=yaml.FullLoader)
        
        # Create output dir if one does not already exist
        os.makedirs(params['generate']['outdir'], exist_ok=True)

        panda_obj = Panda(expression_file=data_paths['exp'], 
                        motif_file=data_paths['motif'], 
                        ppi_file=data_paths['ppi'], 
                        save_tmp=False, 
                        remove_missing=False, 
                        keep_expression_matrix=True, 
                        save_memory=False)

        panda_obj.save_panda_results(params['generate']['outdir'] + "panda_output.txt")
        
        print("PANDA network saved to " + params['generate']['outdir'])
        
        if params['generate']['method'] == 'lioness':
            lion = Lioness(panda_obj, 
                        computing=params['generate']['compute'], 
                        precision="double",
                        ncores=params['generate']['ncores'], 
                        save_dir=params['generate']['outdir'], 
                        save_fmt="npy")

        pickle_path = './tmp/lioness.pickle'
        convert_lion_to_pickle(os.path.join(params['generate']['outdir'], 'panda_output.txt'),
                            params['generate']['outdir'] + 'lioness.npy',
                            "npy", 
                            './tmp/samples.txt',  
                            pickle_path)

        calculate_degree(inputfile=pickle_path,
                         datatype="pickle",
                         outdir=params['generate']['outdir'])

        print("PANDA network saved to " + params['generate']['outdir'])
        print("LIONESS network saved to " + os.path.join(params['generate']['outdir'], "lioness.npy"))
        print("Pickled LIONESS network (for faster loading into Python) saved to " + pickle_path)

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
            survival_analysis(metadata=params["compare"]["survival"]["metadata"],
                            filetype=params["compare"]["survival"]["filetype"], 
                            sampgroup_colname=params["compare"]["survival"]["sampgroup_colname"],
                            alivestatus_colname=params["compare"]["survival"]["alivestatus_colname"],
                            days_colname=params["compare"]["survival"]["days_colname"],
                            groups=params["compare"]["survival"]["groups"],
                            outdir=params["compare"]["survival"]["outdir"])
        
        if args.compchoice == "gsea":    
            perform_gsea(genefile=params["compare"]["gsea"]["genefile"], 
                        gmtfile=params["compare"]["gsea"]["gmtfile"], 
                        geneset=params["compare"]["gsea"]["geneset"], 
                        outdir=params["compare"]["gsea"]["outdir"])
    
    ########################################################
    # 4) Visualize results
    ########################################################       

    elif args.command == "visualize":

        if args.plotchoice == "volcano":    
            plot_volcano(datafile=params["visualize"]["volcano"]["datafile"],
                         fccol=params["visualize"]["volcano"]["FCcol"],
                         adjpcol=params["visualize"]["volcano"]["adjpcol"],
                         fcthreshold=params["visualize"]["volcano"]["fcthreshold"], 
                         adjpvalthreshold=params["visualize"]["volcano"]["adjpvalthreshold"],
                         numlabels=params["visualize"]["volcano"]["numlabels"],
                         outdir=params["visualize"]["volcano"]["outdir"])      
    
        if args.plotchoice == "quantity":    
            plot_expression_degree(datafile=params["visualize"]["quantity"]["datafile"],
                        filetype=params["visualize"]["quantity"]["filetype"], 
                        metadata=params["visualize"]["quantity"]["metadata"],
                        genelist=params["visualize"]["quantity"]["genelist"],
                        plottype=params["visualize"]["quantity"]["plottype"],
                        groups=params["visualize"]["quantity"]["groups"],
                        colors=params["visualize"]["quantity"]["colors"],
                        prefix=params["visualize"]["quantity"]["prefix"],
                        yaxisname=params["visualize"]["quantity"]["yaxisname"],
                        outdir=params["visualize"]["quantity"]["outdir"])   
            
        if args.plotchoice == "heatmap":    
            plot_heatmap(datafile=params["visualize"]["heatmap"]["datafile"],
                        filetype=params["visualize"]["heatmap"]["filetype"], 
                        metadata=params["visualize"]["heatmap"]["metadata"],
                        genelist=params["visualize"]["heatmap"]["genelist"],
                        hierarchicalcluster=params["visualize"]["heatmap"]["hierarchicalcluster"],
                        groups=params["visualize"]["heatmap"]["groups"],
                        prefix=params["visualize"]["heatmap"]["prefix"],
                        plotnames=params["visualize"]["heatmap"]["plotnames"],
                        outdir=params["visualize"]["heatmap"]["outdir"])   
            
    ########################################################
    # 5) Optional, extract edges that connect to specific TFs/genes
    ########################################################

    elif args.command == 'extract':
        
        extract_tfs_genes(pickle=params["extract"]["pickle"], 
                            datatype=args.extractchoice, 
                            namefile=params["extract"]["namefile"], 
                            outdir=params["extract"]["outdir"])
            
