import os,sys
#sys.path.append('/s/chopin/c/proj/protfun/arch/x86_64/lib/python')


import math
import numpy as np
import pickle

import os, sys, numpy, pysam
#from iDiffIR.IntronModel import *
#from iDiffIR.Plot import *
#from iDiffIR.Stat import *
#import matplotlib
#matplotlib.use('agg')
from iDiffIR.IntronModel import *
from SpliceGrapher.formats.fasta import *
from argparse import ArgumentParser, ArgumentTypeError
from SpliceGrapher.shared.utils      import *
from SpliceGrapher.formats.loader import loadGeneModels
import math

import numpy as np
import itertools
import pickle


def get_coordinates(strObj):
    start,end = strObj.split(',')
    start = int(start.strip('('))
    end = int(end.strip(')'))
    return start,end

def determine_length(strObj):
	start,end = get_coordinates(strObj)
	dist = max(start,end)-min(start,end)
	return dist

def is_within(point,part):
	if point >= part[0] and point <= part[1]:
		return True
	else:
		return False

def individual_overalap_check(hotspot,part):
	ret_statement = 'none';
	if is_within(hotspot[0],part) and is_within(hotspot[1],part):
		ret_statement = 'all'
	elif is_within(hotspot[0],part) and is_within(hotspot[1],part)==False:
		ret_statement = 'start'
	elif is_within(hotspot[0],part)==False and is_within(hotspot[1],part):
		ret_statement = 'end'
	elif is_within(hotspot[0],part)==False and is_within(hotspot[1],part)==False and hotspot[0]<part[0] and hotspot[1]>part[1]:
		ret_statement = 'accross'
	return ret_statement
	
def overlap_check(hotspot,pre_exon,intron,post_exon):
	hotspot = get_coordinates(hotspot)
	pre_exon = get_coordinates(pre_exon)
	intron = get_coordinates(intron)
	post_exon = get_coordinates(post_exon)
	even_coordinates = (pre_exon[0],post_exon[1])
	f_pre = individual_overalap_check(hotspot,pre_exon)
	f_in = individual_overalap_check(hotspot,intron)
	f_post = individual_overalap_check(hotspot,post_exon)
	#print hotspot,pre_exon,intron,post_exon
	#print f_pre,f_in,f_post
	#return 'No' if f_pre=='none' else 'Yes','No' if f_in=='none' else 'Yes','No' if f_post=='none' else 'Yes'
	return f_pre,f_in,f_post


def validate_example(gene,hotspot,pre_exon,event,post_exon):
	
	found = False	
	
	pre,intron,post = overlap_check(hotspot,pre_exon,event,post_exon)
	gene_firstExon = (gene.exons[0][0]+gene.minpos,gene.exons[0][1]+gene.minpos)
	gene_lastExon = (gene.exons[-1][0]+gene.minpos,gene.exons[-1][1]+gene.minpos)
				#gene_firstExon = str(gene_firstExon)
				#gene_lastExon = str(gene_lastExon)
				
	#pre_exon = get_coordinates(instance[12])
	#post_exon = get_coordinates(instance[13])
	
	pre_exon = get_coordinates(pre_exon)
	post_exon = get_coordinates(post_exon)
				
	if gene_firstExon==gene_lastExon: #if there is only one exon
		return found
				
					
				
	if pre_exon == gene_firstExon and post_exon != gene_lastExon and (pre=='start' or pre=='all'):
		#count_x1 += 1
		#print pre,intron,post,DHS_control,pre_exon,event,post_exon,gene.minpos,gene.maxpos,len(gene.introns),gene_firstExon,gene_lastExon,'First'
		found = True
		
	elif post_exon == gene_lastExon and pre_exon != gene_firstExon and (post=='end' or post=='all'):
		#count_x2 += 1
		#print pre,intron,post,DHS_control,pre_exon,event,post_exon,gene.minpos,gene.maxpos,len(gene.introns),gene_firstExon,gene_lastExon,'Second'
		found = True
				
	elif pre_exon == gene_firstExon and post_exon == gene_lastExon and (pre=='start' or pre=='all') and (post=='end' or post=='all'):
		#count_x3 += 1
		#print pre,intron,post,DHS_control,pre_exon,event,post_exon,gene.minpos,gene.maxpos,len(gene.introns),gene_firstExon,gene_lastExon,'Third'
		found = True
					
	elif pre_exon != gene_firstExon and post_exon != gene_lastExon:# and found > 0 and found < len(gene.introns)-1:
		#count_x4 += 1
		#print pre,intron,post,DHS_control,pre_exon,event,post_exon,gene.minpos,gene.maxpos,len(gene.introns),gene_firstExon,gene_lastExon,'Fourth'
		found = True
	
	return found

######################################################################################################################################################

#Usage expressed_gene_filtering.py exp_limit file_name dir_name


exp_limit = float(sys.argv[1])
inp_f = sys.argv[2]
dirr = sys.argv[3]
Target = np.loadtxt(inp_f,dtype=str,delimiter='\t')
Exp_File = np.loadtxt('exonic_expression_list.txt',dtype=str)
f_name = inp_f.split('.txt')[0]+'_'+str(exp_limit)+'_expressed'

Exp_File_Dict = {}
for entry in Exp_File[1:]:
	Exp_File_Dict[entry[0]]=entry[1]

Final_List = [Target[0].tolist()]


for entry in Target[1:]:
	if entry[0] in Exp_File_Dict.keys():
		exp = float(Exp_File_Dict[entry[0]])
		if exp >= exp_limit:
			Final_List.append(entry.tolist())

if not os.path.exists(dirr):
    os.makedirs(dirr)

np.savetxt(dirr+'/'+f_name+'.txt',Final_List,fmt='%s',delimiter='\t')	

###############For noBoundaryModified###################
if 'geneRecords' not in locals():
	print 'geneRecord not found'
	with open('Gene_Records_idiffIR.pckl','r') as f:
		geneRecords = pickle.load(f)

#count_nob = 0
Final_List_nob = [Final_List[0]]
data = np.asarray(Final_List)

track_index = []
for gene in geneRecords:
	indices = np.argwhere(data[:,0]==gene.gid)
	if len(indices) > 0:
		for index in indices:
			event = data[index][0][1] #1 is the event in the big table of events
			#event = (event[0]-gene.minpos,event[1]-gene.minpos+1)
			pre_exon = data[index][0][10]
			post_exon = data[index][0][11]
			#DHS_control = 0
			#DHS_heat = 0
			DHS = data[index][0][9]
			
			result = validate_example(gene,DHS,pre_exon,event,post_exon)
			if result == True:
				if index not in track_index:
					track_index.append(index)
					Final_List_nob.append(data[index][0].tolist())
			
np.savetxt(dirr+'/'+f_name+'_noBoundaryModified'+'.txt',Final_List_nob,fmt='%s',delimiter='\t')	
