import numpy as np

def get_params_dict(params_path):
    param_data = np.loadtxt(params_path, dtype=str, delimiter='|')
    params = {}
    for entry in param_data:
        if entry[1] == 'False':
            params[entry[0]] = False
        elif entry[1] == 'True':
            params[entry[0]] = True
        else:
            try:
                params[entry[0]] = int(entry[1])
            except:
                params[entry[0]] = entry[1]    
    
    return params

def calculate_padding(inputLength, filterSize):
    padding = inputLength - (inputLength - filterSize + 1)
    return int(padding/2) #appended to both sides the half of what is needed

def annotate_motifs(annotate_arg, motif_dir):
    ###-----------------Adding TF details to TomTom results----------------###
        try:
            tomtom_res = np.loadtxt(motif_dir+'/tomtom/tomtom.tsv',dtype=str,delimiter='\t')
        except:
            print("Error! motif file not found. Make sure to do motif analysis first.")
            return
        if annotate_arg == 'default':
            database = np.loadtxt('../../Basset_Splicing_IR-iDiffIR/Analysis_For_none_network-typeB_lotus_posThresh-0.60/MEME_analysis/Homo_sapiens_2019_01_14_4_17_pm/TF_Information_all_motifs.txt',dtype=str,delimiter='\t')
        else:
            database = np.loadtxt(annotate_arg, dtype=str, delimiter='\t')
        final = []                                     
        for entry in tomtom_res[1:]:
            motifID = entry[1]                         
            res = np.argwhere(database[:,3]==motifID)
            TFs = ','.join(database[res.flatten(),6])
            final.append(entry.tolist()+[TFs])
        np.savetxt(motif_dir+'/tomtom/tomtom_annotated.tsv', final, delimiter='\t', fmt='%s')