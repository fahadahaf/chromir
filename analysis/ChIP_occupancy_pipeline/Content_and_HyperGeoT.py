from scipy.stats import hypergeom

import numpy as np
import sys

TF_name = sys.argv[1]

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
	data = np.loadtxt(fname,dtype=str,delimiter='\t')
	IR.append(len(data))
	#########For IR_Pop###########
	fname = 'HotSpot_Data/Population_Stats/Union_IR_Events_'+exp_limit[i]+'_expressed_population_IRCovered.txt'
	data = np.loadtxt(fname,dtype=str,delimiter='\t')
	IR_Pop.append(len(data))
	
	#########For IE##########
	fname = 'HotSpot_Data/'+TF_name+'_in_IE_Neg_'+exp_limit[i]+'_expressed_StartInGene_intronCoverage.txt'
	data = np.loadtxt(fname,dtype=str,delimiter='\t')
	IE.append(len(data))
	#########For IE_Pop###########
	fname = 'HotSpot_Data/Population_Stats/Union_IR_Events_Neg_'+exp_limit[i]+'_expressed_population_IRCovered.txt'
	data = np.loadtxt(fname,dtype=str,delimiter='\t')
	IE_Pop.append(len(data))
	
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
	HypergeoPMF.append(hypergeom.sf(IR[i]-1, IE_All_Pop[i], IR_Pop[i], IE_All[i]))

	Final_Res.append([ IR[i],IR_Pop[i],round(IR_Content[i]*100,2), IE[i],IE_Pop[i],round(IE_Content[i]*100,2), HypergeoPMF[i]])

np.savetxt(TF_name+'-some_res_StartInGene.txt',Final_Res,fmt='%s',delimiter='\t')
