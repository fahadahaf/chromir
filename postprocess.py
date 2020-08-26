import glob
import matplotlib.pyplot as plt
import pickle
import numpy as np
import os
import sys

from argparse import ArgumentParser

from utils import get_params_dict


def parseArgs():
    """Parse command line arguments
    
    Returns
    -------
    a : argparse.ArgumentParser
    
    """
    parser = ArgumentParser(description='Post process the ROC and PRC data to generate the corresponding plots.')
    parser.add_argument('-v', '--verbose',dest='verbose', action='store_true', 
                        default=False, help="verbose output [default is quiet running]")
    parser.add_argument('-o','--outDir',dest='out_dir',type=str,
                        action='store',help="output directory. Default: results/ directory (will be created if doesn't exists).", default='results')
    parser.add_argument('-t','--type', dest='type',type=str,
                        action='store',help="Plot type: either ROC or PRC. Default: ROC", default='ROC')
    parser.add_argument('--suffix', dest='suffix',type=str,
                        action='store',help="A unique suffix to add to plot name. Default '' (empty string)", default='')
    parser.add_argument('--curve20',dest='useCurve20', action='store_true', 
                        default=False, help="Plot ROC/PRC cuve at maxed at 0.2 on X-axis (zoom-in version). Default: False")                        					
    parser.add_argument('infofile',type=str,
                        help='The text file containing names and locations of each experiment for which the ROC/PRC curve will be generated.')
    
    args = parser.parse_args()
    return args


def roc_prc_curve(arg_space, exp_dict):
    suffix = '_'+arg_space.suffix if len(arg_space.suffix) > 0 else arg_space.suffix
    curve20 = '_curve20' if arg_space.useCurve20 else ''
    #some colors to be used for individual curves.
    colors = ['darkorange', 'saddlebrown', 'crimson', 'rebeccapurple', 'limegreen', 'teal', 'dimgray']
    out_dir = arg_space.out_dir.strip('/')+'/'

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    pckl_text = ''
    xval,yval = '',''
    areaType = ''
    if arg_space.type == 'ROC':
        areaType = 'AUC'
        pckl_text = 'roc'
        xval,yval = 'fpr','tpr'
        plt.plot([0,1],[0,1],'k--')
    elif arg_space.type == 'PRC':
        areaType = 'AUPRC'
        pckl_text = 'prc'
        xval,yval = 'recall','precision'
        plt.plot([0,1],[0.5,0.5],'k--')
    else:
        print('invalid argument! --type can only have one of the following values: ROC or PRC')
        return

    count = 0
    for key in exp_dict:
        if arg_space.verbose:
            print('Running for: %s', key)
        label = key
        with open(exp_dict[key]+'/modelRes_%s.pckl'%pckl_text, 'rb') as f:
            pckl = pickle.load(f)
        stats = np.loadtxt(exp_dict[key]+'/modelRes_results.txt',delimiter='\t',skiprows=1)
        Xval = pckl[xval]
        Yval =  pckl[yval]
        if arg_space.type == 'ROC':
            test_stat = round(stats[-2],2)
        else:
            test_stat = round(stats[-1],2)
        clr = colors[count]
        plt.plot(Xval, Yval, lw=1, label='%s (%s = %.2f)'%(label,areaType,test_stat), color=clr)
        count += 1

    plt.grid(which='major',axis='both',linestyle='--', linewidth=1)        
    if arg_space.useCurve20:
        plt.xlim(0, 0.2)
        if arg_space.type == 'ROC':
            plt.ylim(0, 0.6)
            plt.xlabel('False positive rate',fontsize=10.5)
            plt.ylabel('True positive rate',fontsize=10.5)
            plt.legend(loc=4, fontsize=10.5)
        else:
            plt.ylim(0.5, 1)
            plt.xlabel('Recall',fontsize=10.5)
            plt.ylabel('Precision',fontsize=10.5)
            plt.legend(loc=1, fontsize=10.5)
        #plt.title('Precision-Recall curves')
    else:
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        if arg_space.type == 'ROC':
            plt.xlabel('False positive rate',fontsize=10.5)
            plt.ylabel('True positive rate',fontsize=10.5)
            plt.legend(loc=4, fontsize=10.5)
        else:
            plt.xlabel('Recall',fontsize=10.5)
            plt.ylabel('Precision',fontsize=10.5)
            plt.legend(loc=3, fontsize=10.5)
        #plt.title('Precision-Recall curves')
    plt.savefig(out_dir+'%s_curves_selected%s%s.pdf'%(pckl_text.upper(),curve20,suffix))
    plt.savefig(out_dir+'%s_curves_selected%s%s.png'%(pckl_text.upper(),curve20,suffix))
    plt.clf()


def main():
    arg_space = parseArgs()
    #create params dictionary
    params_dict = get_params_dict(arg_space.infofile)
    #print(params_dict)
    roc_prc_curve(arg_space, params_dict)


if __name__ == "__main__":
    main()

