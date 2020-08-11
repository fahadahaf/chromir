import numpy as np
from utils import get_params_dict


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


def evaluateRegular(net, iterator, criterion, out_dirc, getCNN=False, storeCNNout = False, getSeqs = False):
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
                    output_dir = out_dirc
                    if not os.path.exists(output_dir):
                        os.makedirs(output_dir)	
                    with open(output_dir+'/CNNout_batch-'+str(batch_idx)+'.pckl','wb') as f:
	                    pickle.dump(outputCNN.cpu().detach().numpy(),f)
                    per_batch_CNNoutput[batch_idx] = output_dir+'/CNNout_batch-'+str(batch_idx)+'.pckl'
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




data = np.loadtxt('models/BASSET-noEmbds_hyperParams.txt',dtype=str,delimiter='\t')
get_params_dict(data)