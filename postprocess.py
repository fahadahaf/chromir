import glob
import matplotlib.pyplot as plt
import pickle
import numpy as np
import sys


def some_processing_func():
    #final_list = ['Basset-Embds','Basset-noEmbds','MH-PosEnc','CNN-RNN-MH-noPosEnc','RNN-MH-noPosEnc']
    final_list = {'Basset-noEmbds':['Basset', 'darkorange'],
                'Basset-Embds':['Basset-E','saddlebrown'],
                'MH-PosEnc':['MSA','crimson'],
                'CNN-MH-noPosEnc':['CNN-MSA','rebeccapurple'],
                'CNN-RNN-MH-noPosEnc':['CNN-RNN-MSA','limegreen'],
                'RNN-MH-noPosEnc':['RNN-MSA','dodgerblue']}

    all_pickles = glob.glob('Finalzed_Results_TrainTest/*roc.pckl') 

    file_dict = {}
    for entry in all_pickles:
        pck_name = entry.split('/')[1].split('_')[0]
        file_dict[pck_name] = entry

    base_FPR = [round(i*0.01,2) for i in range(0,101)]

    plt.plot([0,1],[0,1],'k--')

    ###########for GKMSM###############
    gkm_file = '/s/jawar/p/nobackup/altsplice1/fahad/Gapped_SVM_Analysis/Encode_IR_DHSs-idiffIR/TrainTest_Split/ROC_data_FPR-TPR.txt'
    gkm_data = np.loadtxt(gkm_file,dtype=float)

    TPR = gkm_data[:,1]
    FPR = gkm_data[:,0]

    exp_TPR = []
        
    for i in range(0, len(base_FPR)):
        index = np.argwhere(FPR >= base_FPR[i])[0][0]
        exp_TPR.append(TPR[index])

    test_auc = np.loadtxt('/s/jawar/p/nobackup/altsplice1/fahad/Gapped_SVM_Analysis/Encode_IR_DHSs-idiffIR/TrainTest_Split/final_results_TrainTest.txt',dtype=str,skiprows=1)[5]
    test_auc = round(float(test_auc),2)
    plt.plot(base_FPR,exp_TPR,lw=1,label='Gapped SVM (AUC = '+str(test_auc)+')',color='dimgrey')
    ###################################

    for key in file_dict:
        if key not in final_list:
            continue
        print("Running for: ",key)
        pckl_name = file_dict[key]
        with open(pckl_name,'rb') as f:
            roc_dict = pickle.load(f)
        
        TPR = roc_dict['tpr']
        FPR = roc_dict['fpr']
        
        exp_TPR = []
        
        for i in range(0, len(base_FPR)):
            index = np.argwhere(FPR >= base_FPR[i])[0][0]
            exp_TPR.append(TPR[index])
        
        res_data = np.loadtxt('Finalzed_Results_TrainTest/'+key+'_results.txt',delimiter='\t',skiprows=1)
        
        test_auc = round(res_data[-2],2)
        value = final_list[key][0]
        clr = final_list[key][1]
        plt.plot(base_FPR,exp_TPR,lw=1,label=value+' (AUC = '+str(test_auc)+')',color=clr)
            
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.title('ROC curves')
    plt.legend(loc=4, fontsize=9)#'best')
    plt.grid(which='major',axis='both',linestyle='--', linewidth=1)
    plt.savefig('Finalzed_Results_TrainTest/ROC_curves_selected.pdf')
    plt.savefig('Finalzed_Results_TrainTest/ROC_curves_selected.png')
    plt.clf()
