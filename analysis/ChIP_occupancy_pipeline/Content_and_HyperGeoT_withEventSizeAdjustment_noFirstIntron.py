from scipy.stats import hypergeom

import numpy as np
import pandas as pd
import sys
import re

TF_name = sys.argv[1]

def getIntronLens(data):
	lens = []
	if len(data[0])>5:
		intron_index = 1
	else:
		intron_index = 2
	for i in range(0,len(data)):
		entry = data[i]                  
		start = int(entry[intron_index].split(',')[0].split('(')[1])
		end = int(entry[intron_index].split(',')[1].split(')')[0])
		lens.append((end-start)+1)
	return np.asarray(lens)

                                                                                                                                                                                                 
                                                                                          
df = pd.read_csv('/s/jawar/p/nobackup/altsplice1/fahad/Human_hg19/Annotations/MISOHelp/all_introns_ensemble.txt',sep='\t')                                                                                          
df = df[df['isFirst']==True].reset_index(drop=True)
df['Coordinates'] = df['Coordinates'].apply(lambda x: str((int(re.split('\(|\,|\)|',x)[1]), int(re.split('\(|\,|\)|',x)[2])-1)))  	
df['Key'] = df.apply(lambda x: x['Gene']+'_'+x['Coordinates'], axis=1)

def filter_first_intron(data, df, pop=False):
	print('filtering for first Intron...')
	indices = []
	for i in range(0, data.shape[0]):
			gene = data[i][0]
			intron = data[i][1] if pop == False else data[i][2]
			if gene+'_'+intron in df['Key'].values:
				indices.append(i)
	return np.delete(data, indices, axis=0)

exp_limit = ['1.0','5.0','10.0','20.0']

IR = []
IR_Pop = []

IE = []
IE_Pop = []

IE_All = []
IE_All_Pop = []

IR_Content = []
IE_Content = []
HypergeoPMF = []
Final_Res = []

for i in range(0,len(exp_limit)):
	#########For IR#########
	fname = 'HotSpot_Data/'+TF_name+'_in_IR_'+exp_limit[i]+'_expressed_StartInGene_intronCoverage.txt'
	dataIR = np.loadtxt(fname,dtype=str,delimiter='\t', skiprows=1)
	dataIR = filter_first_intron(dataIR, df, pop=False)
	IR.append(len(dataIR))
	#########For IR_Pop###########
	fname = 'HotSpot_Data/Population_Stats/Union_IR_Events_'+exp_limit[i]+'_expressed_population_IRCovered.txt'
	dataIRpop = np.loadtxt(fname,dtype=str,delimiter='\t', skiprows=1)
	dataIRpop = filter_first_intron(dataIRpop, df, pop=True)
	IR_Pop.append(len(dataIRpop))
	
	IRpopLens = getIntronLens(dataIRpop)
	maxLen = np.mean(IRpopLens) + (1*np.std(IRpopLens)) #mean + 2 St dev (previously was mean + 1 st dev)
	
	
	#########For IE##########
	fname = 'HotSpot_Data/'+TF_name+'_in_IE_Neg_'+exp_limit[i]+'_expressed_StartInGene_intronCoverage.txt'
	dataIE = np.loadtxt(fname,dtype=str,delimiter='\t', skiprows=1)
	dataIE = filter_first_intron(dataIE, df, pop=False)
	#IE.append(len(dataIE))
	
	IEIntronLens = getIntronLens(dataIE)
	indicesIE = np.where(IEIntronLens < maxLen)[0]
	
	dataIE = dataIE[indicesIE]
	IE.append(len(dataIE))
	
	
	#########For IE_Pop###########
	fname = 'HotSpot_Data/Population_Stats/Union_IR_Events_Neg_'+exp_limit[i]+'_expressed_population_IRCovered.txt'
	dataIEpop = np.loadtxt(fname,dtype=str,delimiter='\t', skiprows=1)
	dataIEpop = filter_first_intron(dataIEpop, df, pop=True)
	#IE_Pop.append(len(dataIEpop))
	
	IEpopIntronLens = getIntronLens(dataIEpop)
	indicesIEpop = np.where(IEpopIntronLens < maxLen)[0]
	
	dataIEpop = dataIEpop[indicesIEpop]
	IE_Pop.append(len(dataIEpop))
	
	#print len(dataIR)/float(len(dataIRpop)),len(indicesIE)/float(len(indicesIEpop))
	print len(dataIR)/float(len(dataIRpop)),len(dataIE)/float(len(dataIEpop))
	
	#########For IE_All##########
	#fname = 'HotSpot_Data/EGR1_in_IE_Neg_ALL_'+exp_limit[i]+'_expressed.txt'
	#data = np.loadtxt(fname,dtype=str)
	#IE_All.append(len(data))
	IE_All.append(IR[i]+IE[i])
	#########For IE_All_Pop###########
	#fname = 'HotSpot_Data/Union_IR_Events_Neg_ALL_'+exp_limit[i]+'_expressed_population.txt'
	#data = np.loadtxt(fname,dtype=str)
	#IE_All_Pop.append(len(data))
	IE_All_Pop.append(IR_Pop[i]+IE_Pop[i])
	
	##############IR Content################
	IR_Content.append(float(IR[i])/float(IR_Pop[i]))
	
	##############IE Content################
	IE_Content.append(float(IE[i])/float(IE_Pop[i]))
	
	##############Hypergeometric Test###################
	#HypergeoPMF.append(hypergeom.pmf(IR[i],IE_All_Pop[i],IR_Pop[i],IE_All[i]))
	HypergeoPMF.append(hypergeom.sf(IR[i]-1, IE_All_Pop[i], IR_Pop[i], IE_All[i])) #see medium article understanding and implementing the hypergeo test in python
	
	Final_Res.append([ IR[i],IR_Pop[i],round(IR_Content[i]*100,2), IE[i],IE_Pop[i],round(IE_Content[i]*100,2), HypergeoPMF[i]])

np.savetxt(TF_name+'-some_res_EventLenAdj_StartInGene_noFirstIntron.txt',Final_Res,fmt='%s',delimiter='\t')
