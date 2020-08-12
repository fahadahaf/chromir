import gensim
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim

from argparse import ArgumentParser
from fastprogress import progress_bar
from gensim.models import Word2Vec
from gensim.models.word2vec import LineSentence
from random import randint
from sklearn import metrics
from torch.backends import cudnn
from torch.utils import data
from torch.utils.data import Dataset, DataLoader
from torch.autograd import Variable
from torch.utils.data.sampler import SubsetRandomSampler
from torch.autograd import Function # import Function to create custom activations
from torch.nn.parameter import Parameter # import Parameter to create custom activations with learnable parameters
from torch import optim # import optimizers for demonstrations

#local imports
from datasets import DatasetLoadAll, DatasetLazyLoad, DatasetEmbd
from models import Basset, AttentionNet
from utils import get_params_dict



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
                        default=None, help="Annotate tomtom motifs. The options are: 1. path to annotation file, 2. No (not to annotate the output) 3. None (default where human CISBP annotations are used)")
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
    parser.add_argument('inputprefix', type=str,
                        help="Input file prefix for the bed/text file and the corresponding fasta file (sequences).")
    parser.add_argument('hparamfile',type=str,
                        help='Name of the hyperparameters file to be used.')
    
    args = parser.parse_args()
    return args
###########################################################################################################################
#---------------------------------------------------------End-------------------------------------------------------------#
###########################################################################################################################


###########################################################################################################################
#--------------------------------------------Train and Evaluate Functions-------------------------------------------------#
###########################################################################################################################
def trainRegular(model, device, iterator, optimizer, criterion):
    model.train()
    running_loss = 0.0
    train_auc = []
    for batch_idx, (headers, seqs, data, target) in enumerate(iterator):
        #pdb.set_trace()
        data, target = data.to(device,dtype=torch.float), target.to(device,dtype=torch.long)
        optimizer.zero_grad()
        outputs,_ = model(data)
        loss = criterion(outputs, target)
        #loss = F.binary_cross_entropy(outputs, target)
        labels = target.cpu().numpy()
        
        softmax = torch.nn.Softmax(dim=1)
        pred = softmax(outputs)
        pred = pred.cpu().detach().numpy()
        #print(pred)
        try:
            train_auc.append(metrics.roc_auc_score(labels, pred[:,1]))
        except:
            train_auc.append(0.0)

        loss.backward()
        optimizer.step()
        running_loss += loss.item()
        #return outputs
    return running_loss/len(train_loader),train_auc


def evaluateRegular(net, device, iterator, criterion, out_dirc=None, getCNN=False, storeCNNout = False, getSeqs = False):
    running_loss = 0.0
    valid_auc = [] 
    roc = np.asarray([[],[]]).T
    per_batch_labelPreds = {}
    per_batch_CNNoutput = {}
    per_batch_testSeqs = {}
    per_batch_info = {}

    net.eval()
    CNNlayer = net.layer1[0:3]
    CNNlayer.eval()
    
    with torch.no_grad():
        for batch_idx, (headers, seqs, data, target) in enumerate(iterator):
            data, target = data.to(device,dtype=torch.float), target.to(device,dtype=torch.long)
            # Model computations
            outputs = net(data)
            loss = criterion(outputs, target)
            softmax = torch.nn.Softmax(dim=1)
            labels=target.cpu().numpy()
            pred = softmax(outputs)
            pred=pred.cpu().detach().numpy()
            label_pred = np.column_stack((labels,pred[:,1]))
            per_batch_labelPreds[batch_idx] = label_pred
            roc = np.row_stack((roc,label_pred))

            try:
                valid_auc.append(metrics.roc_auc_score(labels, pred[:,1]))
            except:
                valid_auc.append(0.0)
            running_loss += loss.item()
        
            outputCNN = CNNlayer(data).cpu().detach().numpy()
            if getCNN == True:
                outputCNN = CNNlayer(data)
                if storeCNNout == True:
                    if not os.path.exists(out_dirc):
                        os.makedirs(out_dirc)	
                    with open(out_dirc+'/CNNout_batch-'+str(batch_idx)+'.pckl','wb') as f:
	                    pickle.dump(outputCNN.cpu().detach().numpy(),f)
                    per_batch_CNNoutput[batch_idx] = out_dirc+'/CNNout_batch-'+str(batch_idx)+'.pckl'
                else:
                    per_batch_CNNoutput[batch_idx] = outputCNN.cpu().detach().numpy()
            
            if getSeqs == True:
                per_batch_testSeqs[batch_idx] = np.column_stack((headers,seqs))
            
    labels = roc[:,0]
    preds = roc[:,1]
    valid_auc = metrics.roc_auc_score(labels,preds)
        
    return running_loss/len(iterator),valid_auc,roc,per_batch_labelPreds,per_batch_CNNoutput,per_batch_testSeqs 
###########################################################################################################################
#---------------------------------------------------------End-------------------------------------------------------------#
###########################################################################################################################
def get_indices(dataset_size, test_split, output_dir, shuffle_data=True, seed_val=100):
    indices = list(range(dataset_size))
    split_val = int(np.floor(test_split*dataset_size))
    if shuffle_data:
        np.random.seed(seed_val)
        np.random.shuffle(indices)
    train_indices, test_indices, valid_indices = indices[2*split_val:], indices[:split_val], indices[split_val:2*split_val]
    #--save indices for later use, when testing for example---#
    np.savetxt(output_dir+'/valid_indices.txt',valid_indices,fmt='%s')
    np.savetxt(output_dir+'/test_indices.txt',test_indices,fmt='%s')
    np.savetxt(output_dir+'/train_indices.txt',train_indices,fmt='%s')
    return train_indices, test_indices, valid_indices


def run_experiment(device, arg_space, output_dir, params):
    input_prefix = arg_space.inputprefix
    w2v_path = arg_space.wvPath #'Word2Vec_Models/'
    net_type = arg_space.netType

    num_labels = params['num_classes']
    get_CNNout = params['get_CNNout']
    get_sequences = params['get_seqs']
    batch_size = params['batch_size']
    max_epochs = params['num_epochs']
    use_embds = params['use_embeddings']
    prefix = 'modelRes' #Using generic, not sure if we need it as an argument or part of the params dict

    test_split = arg_space.splitperc/100
    if arg_space.verbose:
        print("test/validation split val: %.2f"%test_split)

    if use_embds == False:
        if arg_space.deskLoad:
            final_dataset = DatasetLazyLoad(input_prefix)
        else:
            final_dataset = DatasetLoadAll(input_prefix)

        train_indices, test_indices, valid_indices = get_indices(len(final_dataset), test_split, output_dir, shuffle_data=True, seed_val=100)
        modelvw = None
    else:
        data_all = DatasetLoadAll(input_prefix, for_embeddings=True)
        train_indices, test_indices, valid_indices = get_indices(len(final_dataset), test_split, output_dir, shuffle_data=True, seed_val=100)
        final_dataset = pd.merge(data_all.df_seq_final, data_all.df, on='header')[['sequence',7]]
        final_dataset = DatasetEmbd(final_dataset.values.tolist(),modelwv,kmer_len,stride)
        w2v_filename = 'Word2Vec_Model_kmerLen'+str(params['embd_kmersize'])+'_win'+str(params['embd_window'])+'_embSize'+str(params['embd_size'])
        modelwv = Word2Vec.load(w2v_path+w2v_filename)

    train_sampler = SubsetRandomSampler(train_indices)
    test_sampler = SubsetRandomSampler(test_indices)
    valid_sampler = SubsetRandomSampler(valid_indices)
    train_loader = DataLoader(data_all, batch_size = batchSize, sampler = train_sampler)
    test_loader = DataLoader(data_all, batch_size = batchSize, sampler = test_sampler)
    valid_loader = DataLoader(data_all, batch_size = batchSize, sampler = valid_sampler)

    if net_type == 'basset':
        net = Basset(wvmodel=modelwv).to(device)
    else:
        net = AttentionNet(wvmodel=modelwv).to(device)
    
    criterion = nn.CrossEntropyLoss(reduction='mean')
    optimizer = optim.Adam(net.parameters())

    ##-------Main train/test loop----------##
    if arg_space.mode == 'train':
        best_valid_loss = np.inf
        best_valid_auc = np.inf
        for epoch in progress_bar(range(1, max_epochs + 1)):
            res_train = trainRegular(net, device, train_loader, optimizer, criterion)
            res_valid = evaluateRegular(net, device, valid_loader, criterion)
            res_train_auc = np.asarray(res_train[1]).mean()
            res_train_loss = res_train[0]
            res_valid_auc = np.asarray(res_valid[1]).mean()
            res_valid_loss = res_valid[0]
            if res_valid_loss < best_valid_loss:
                best_valid_loss = res_valid_loss
                best_valid_auc = res_valid_auc
                if arg_space.verbose:
                    print("Best Validation Loss: %.3f and AUC: %.2f"%(best_valid_loss, best_valid_auc), "\n")
                torch.save({'epoch': epoch,
                        'model_state_dict': net.state_dict(),
                        'optimizer_state_dict':optimizer.state_dict(),
                        'loss':res_valid_loss
                        },output_dir+'/'+prefix+'_model')
    else:
        try:    
            checkpoint = torch.load(output_dir+'/model')
            net.load_state_dict(checkpoint['model_state_dict'])
            optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
            epoch = checkpoint['epoch']
            loss = checkpoint['loss']
        except:
            print("No pre-trained model found! Please run with --mode set to train.")
        res_test = evaluateRegular(net, device, test_loader, criterion, output_dir+"/Stored_Values",
                                   getCNN=get_CNNout, storeCNN=argParse.storeCNN, getSeqs=get_sequences)
        test_loss = res_test[0]
        labels = res_test[2][:,0]
        preds = res_test[2][:,1]
        auc_test = metrics.roc_auc_score(labels, preds)
        if arg_space.verbose:
            print("Test Loss: %.3f and AUC: %.2f"%(test_loss, auc_test), "\n")
        auprc_test = metrics.average_precision_score(labels,preds)

        some_res = [['Best_Valid_Loss','Best_Valid_AUC','Test_Loss','Test_AUC', 'Test_AUPRC']]
        some_res.append([best_valid_loss,best_valid_auc,test_loss,auc_test,auprc_test])

        fpr,tpr,thresholds = metrics.roc_curve(labels,preds)
        precision,recall,thresholdsPR = metrics.precision_recall_curve(labels,preds)
        roc_dict = {'fpr':fpr, 'tpr':tpr, 'thresholds':thresholds}
        prc_dict = {'precision':precision, 'recall':recall, 'thresholds':thresholdsPR}

        with open(output_dir+'/'+prefix+'_roc.pckl','wb') as f:
	        pickle.dump(roc_dict,f)
        with open(output_dir+'/'+prefix+'_prc.pckl','wb') as f:
	        pickle.dump(prc_dict,f)
        np.savetxt(output_dir+'/'+prefix+'_results.txt',some_res,fmt='%s',delimiter='\t')

        return res_test
        


def main():
    #CUDA for pytorch
    use_cuda = torch.cuda.is_available()
    device = torch.device(torch.cuda.current_device() if use_cuda else "cpu")
    cudnn.benchmark = True

    arg_space = parseArgs()
    output_dir = arg_space.directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    #save arguments to keep record
    with open(output_dir+'/arguments.txt','w') as f:
	    f.writelines(str(argSpace))

    #create params dictionary
    params_dict = get_params_dict(argSpace.hparamfile)

    test_resBlob = run_experiment(device, arg_space, output_dir, params_dict)


if __name__ == "__main__":
    main()
