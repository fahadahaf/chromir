import torch

from argparse import ArgumentParser
from torch.backends import cudnn

#local imports
from experiment import motif_analysis, run_experiment
from utils import annotate_motifs, get_params_dict



###########################################################################################################################
#--------------------------------------------------Argument Parser--------------------------------------------------------#
###########################################################################################################################
def parseArgs():
    """Parse command line arguments
    
    Returns
    -------
    a : argparse.ArgumentParser
    
    """
    parser = ArgumentParser(description='Main chromIR script to run experiments.')
    parser.add_argument('-v', '--verbose',dest='verbose', action='store_true', 
                        default=False, help="verbose output [default is quiet running]")
    parser.add_argument('-o','--outDir',dest='directory',type=str,
                        action='store',help="output directory", default='')
    parser.add_argument('-m','--mode', dest='mode',type=str,
                        action='store',help="Mode of operation: train or test.", default='train')     
    parser.add_argument('--deskload', dest='deskLoad',
                        action='store_true',default=False,
                        help="Load dataset from desk. If false, the data is converted into tensors and kept in main memory (not recommended for large datasets).")  
    parser.add_argument('-w','--numworkers',dest='numWorkers',type=int,
                        action='store',help="Number of workers used in data loader. For loading from the desk, use more than 1 for faster fetching.", default=1)        
    parser.add_argument('--splitperc',dest='splitperc',type=float, action='store',
                        help="Pecentages of test, and validation data splits, eg. 10 for 10 percent data used for testing and validation.", default=10)
    parser.add_argument('--motifanalysis', dest='motifAnalysis',
                        action='store_true',default=False,
                        help="Analyze CNN filters for motifs and search them against known TF database.")
    parser.add_argument('--scorecutoff',dest='scoreCutoff',type=float,
                        action='store',default=0.65,
                        help="In case of binary labels, the positive probability cutoff to use.")
    parser.add_argument('--tomtompath',dest='tomtomPath',
                        type=str,action='store',default=None,
                        help="Provide path to where TomTom (from MEME suite) is located.") 
    parser.add_argument('--database',dest='tfDatabase',type=str,action='store',
                        help="Search CNN motifs against known TF database. Default is Human CISBP TFs.", default=None)
    parser.add_argument('--annotate',dest='annotateTomTom',type=str,action='store',
                        default=None, help="Annotate tomtom motifs. The value of this variable should be path to the database file used for annotation. Default is None.")
    parser.add_argument('-s','--store', dest='storeCNN',
                        action='store_true',default=False,
                        help="Store per batch CNN outpout matrices. If false, the are kept in the main memory.")
    parser.add_argument('--tomtomdist', dest='tomtomDist',type=str,
                        action='store',default = 'pearson',
                        help="TomTom distance parameter (pearson, kullback, ed etc). Default is pearson. See TomTom help from MEME suite.")
    parser.add_argument('--tomtompval', dest='tomtomPval',type=float,
                        action='store',default = 0.05,
                        help="Adjusted p-value cutoff from TomTom. Default is 0.05.")
    parser.add_argument('--wvpath',dest='wvPath',
                        type=str,action='store',default=None,
                        help="Path to where the word2vec trained embeddings are located. Default is None.")
    parser.add_argument('--nettype',dest='netType',
                        type=str,action='store',default='basset',
                        help="Model type to use: either basset or attention. Default is basset.")  
    parser.add_argument('--embd', dest='useEmbeddings',
                        action='store_true',default=False,
                        help="Whether to use word2vec embeddings instead of one-hot encoded input. Default is False. Make sure to provide the appropriate hyperparameters file.")                       					
    parser.add_argument('inputprefix', type=str,
                        help="Input file prefix for the bed/text file and the corresponding fasta file (sequences).")
    parser.add_argument('hparamfile',type=str,
                        help='Name of the hyperparameters file to be used.')
    
    args = parser.parse_args()
    return args
###########################################################################################################################
#---------------------------------------------------------End-------------------------------------------------------------#
###########################################################################################################################



def main():
    #CUDA for pytorch
    use_cuda = torch.cuda.is_available()
    device = torch.device(torch.cuda.current_device() if use_cuda else "cpu")
    cudnn.benchmark = True
    arg_space = parseArgs()
    #create params dictionary
    params_dict = get_params_dict(arg_space.hparamfile)
    test_resBlob, CNNWeights = run_experiment(device, arg_space, params_dict)
    if arg_space.motifAnalysis:
        motif_dir_pos, num_examples_pos = motif_analysis(test_resBlob, CNNWeights, arg_space)
        motif_dir_neg, num_examples_neg = motif_analysis(test_resBlob, CNNWeights,  arg_space, for_negative=True)
    
    if arg_space.annotateTomTom != None:
        if arg_space.verbose:
            print("Annotating motifs...")
        annotate_motifs(arg_space.annotateTomTom, motif_dir_pos)
        annotate_motifs(arg_space.annotateTomTom, motif_dir_neg)

        
 


if __name__ == "__main__":
    main()
