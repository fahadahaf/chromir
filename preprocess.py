
import numpy as np
import random
import sys

def check_overlap(element,DHS):
        sE = element[0]
        eE = element[1]
        sD = DHS[0]
        eD = DHS[1]
        return (sE <= eD) and (eE >= sD)

def process_IR(data, intragenic, save_results=False, verbose=True, gene_cols=[0,-2]):

    count = 0
    final_list = []
    final_list_v2 = []
    final_list_v3 = []
    final_list_v4 = []

    for i in range(0,len(intragenic)):
        entry = intragenic[i]
        for element in data:
            if entry[gene_cols[0]] == element[gene_cols[1]]:
                event = element[2].strip('(').strip(')').split(',')
                event = [int(event[0]),int(event[1])]
                
                preExon = element[1].strip('(').strip(')').split(',')
                postExon = element[3].strip('(').strip(')').split(',')
                eventWhole = [int(preExon[0]),int(postExon[1])]
                #event = [int(element[1]),int(element[2])]
                DHS = [int(entry[1]),int(entry[2])]
                #res = check_overlap(event,DHS)
                
                if check_overlap(event,DHS): #should overlap the retained intron ____
                    if [i]+entry.tolist() not in final_list:
                        final_list.append([i]+entry.tolist())
                                    
                elif check_overlap(eventWhole,DHS): #should be overlapping the whole event ===____==
                    if [i]+entry.tolist() not in final_list_v2:
                        final_list_v2.append([i]+entry.tolist())
        count += 1
        if verbose:
            if count % 100 == 0:
                print(count,len(final_list),len(final_list_v2))

    if save_results:                    
        np.savetxt('final_overlapping_DHSs_IR_iDiffIR.txt',final_list,fmt='%s',delimiter='\t') #events overlapping retained introns
        np.savetxt('final_overlapping_DHSs_IR_iDiffIR_V2.txt',final_list_v2,fmt='%s',delimiter='\t') #events overlapping the event (flanking exons) but not the retained intron
    
    return final_list, final_list_v2

def remove_overlap_duplicates(data, verbose=True):

    final = [data[0,1:].tolist()]
    final_indices = [data[0,0]]
    count = 0
    for entry in data:
        count += 1
        if verbose:
            if count%1000 == 0:
                print(f"Removing duplicates, done: {count}")
        if entry[1:].tolist() in final:
            continue
        else:
            final.append(entry[1:].tolist())
            final_indices.append(entry[0])

    final_data = np.column_stack((final_indices,final))
    #outfname = fname.strip('.txt')+'_uniqNew.txt'
    #np.savetxt(outfname,final_data,fmt='%s')
    return final_data

def create_dataset(IRdata, intragenicDHSs, load_from_file=False, save_data=True, remove_duplicates=True):
    if load_from_file:
        try:
            pos_data_V1 = np.loadtxt('final_overlapping_DHSs_IR_iDiffIR.txt',dtype=str) #these events (true positive) are those which overlap a DHS with the retained intron
            pos_data_V2 = np.loadtxt('final_overlapping_DHSs_IR_iDiffIR_V2.txt',dtype=str) #these events have a DHS not in retained intron but in the flanking exons so we still don't want any DHS overlapping these.
        except FileNotFoundError:
            print("DHS and IR overlapping files not found! set load_from_file=False to re-run the overlapping analysis.")
    else:
        pos_data_V1, pos_data_V2 = process_IR(IRdata, intragenicDHSs, save_results=False)

    if remove_duplicates:
        pos_data_V1 = remove_overlap_duplicates(pos_data_V1)
        pos_data_V2 = remove_overlap_duplicates(pos_data_V2)

    rest_data = intragenicDHSs


    #pos_data and pos_data_V2 are mutually exclusive so combine them when getting making sure no negative examples are selected mixed with positives
    pos_data = np.concatenate((pos_data_V1,pos_data_V2))

    pos_indices = pos_data[:,0].astype(int)
    rest_noPos = np.delete(rest_data,pos_indices,axis=0)



    pos_final = pos_data_V1[:,1:] #we don't care about V2 since its not going to be part of the final positive set
    neg_indices = np.random.choice([i for i in range(0,len(rest_noPos))],size=len(pos_final)*2,replace=False)

    #pos_final = pos_data_V1[:,1:] #we don't care about V2 since its not going to be part of the final positive set
    neg_final = rest_noPos[neg_indices]



    ############This step is to make sure we don't have information leakage between positive and negative examples (examples appearing in both cases)############
    #pos_final_test = pos_final[:,:-1].tolist() 
    pos_final_test = pos_final.tolist() 

    pos_set = []
    count = 0
    for entry in neg_final:
        if entry.tolist() in pos_final_test:
            count += 1
            pos_set.append(entry.tolist())

    pos_final_last = []

    for entry in pos_final:
        if entry.tolist() not in pos_set:
            pos_final_last.append(entry.tolist())

    pos_final_last = np.asarray(pos_final_last)
    ##########################################################################

    pos_labels = np.asarray(['1' for i in range(0,len(pos_final_last))])
    neg_labels = np.asarray(['0' for i in range(0,len(neg_final))])

    pos_final_last = np.column_stack((pos_final_last,pos_labels))
    neg_final = np.column_stack((neg_final,neg_labels))

    data_final = np.concatenate((pos_final_last,neg_final),axis=0)

    #np.savetxt('Labelled_Data_'+dtype+'.txt',data_final,fmt='%s',delimiter='\t')
    if save_data:
        fname = 'data/Labelled_Data_IR_iDiffIR_corrected.txt'
        np.savetxt(fname,data_final,fmt='%s',delimiter='\t')
        return fname, data_final
    else:
        return data_final

def main():
    IRdatafile = sys.argv[1] #'/s/jawar/p/nobackup/altsplice1/fahad/Human_hg19/Annotations/MISOHelp/IR_events_iDiffIR.txt'
    intragenicDHSfile = sys.argv[2] #'../encode_roadmap_inGenesOnly_orig.bed'
    IRdata = np.loadtxt(IRdatafile, dtype=str, skiprows=1, delimiter='\t')
    intragenicDHSs = np.loadtxt(intragenicDHSfile, dtype=str)

    result = create_dataset(IRdata, intragenicDHSs, load_from_file=False, save_data=True, remove_duplicates=True)

    if isinstance(result, str):
        fname = result
    else:
        fname = result[0]
    
    print(f'Done! IR and non-IR DHSs are at: {fname}')
    


if __name__ == "__main__":
    main()