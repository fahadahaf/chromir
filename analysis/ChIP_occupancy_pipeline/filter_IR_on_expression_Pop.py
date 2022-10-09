import os, sys, numpy
#from iDiffIR.IntronModel import *
#from iDiffIR.Plot import *
from iDiffIR.Stat import *
#import matplotlib
#matplotlib.use('agg')
from SpliceGrapher.formats.fasta import *
from argparse import ArgumentParser, ArgumentTypeError
from SpliceGrapher.shared.utils      import *
from SpliceGrapher.formats.loader import loadGeneModels
from SpliceGrapher.formats.GeneModel   import *
#from SpliceGrapher.formats.loadGTF     import *
from SpliceGrapher.formats.fasta       import FastaRecord
from SpliceGrapher.SpliceGraph         import *
from SpliceGrapher.formats.GTFLoader   import *
from SpliceGrapher.formats.loader import loadGeneModels
from SpliceGrapher.formats.FastaLoader import FastaLoader
import math
import matplotlib.pyplot as plt
from itertools import groupby
#plt.switch_backend('Qt4Agg') #somehow Mike's code mess up backend so use this to view plots

import numpy as np
import pickle
import re

#This script is for the population of IR events (all IR), the other one is when IR and DHSs overlap.


###############################################################################################
print ("Now loading gene dictionary...")
with open ('ALL_GENE_INFO_DICT.pckl','rb') as f:
    gene_dict = pickle.load(f)
print ("Loading complete!")

control = '/s/jawar/h/nobackup/fahad/Human_Chromatin/RNAseq/K562/wgEncodeCaltechRnaSeqK562R2x75Il200Aligns_Merged.bam'
###############################################################################################

def determine_length(strObj):
	start,end = get_coordinates(strObj)
	dist = max(start,end)-min(start,end)+1
	return dist

def get_coordinates(strObj):
    start,end = strObj.split(',')
    start = int(start.strip('('))
    end = int(end.strip(')'))
    return start,end

def check_coverage(lst,perc):
	length = len(lst)
	count = numpy.count_nonzero(lst)
	perc_calc = (float(count)/length) * 100
	#print perc_calc,count,length
	if perc_calc >= perc:
		return True,perc_calc
	else:
		return False,perc_calc

def coverage_stats(bamfile,chromosome,start,end,cov_perc):
	exp_array = getDepthsFromBam(bamfile,chromosome,start,end)[0] #[0] is to get the actual array
	if len(exp_array)==0:
		return -1,-1
	result,perc = check_coverage(exp_array,cov_perc)
	return result,perc
	
	
fname = sys.argv[1]
data = np.loadtxt(fname,dtype=str,delimiter='\t')
cov_perc = 100 # per cent of coverage required (X% of locations should have at least one or more reads)
fname = fname.split('.txt')[0]

percentages = [] #to see distribution of percentages of IR
count_control = 0 #counting track of IR events in control qualifying the minimum coverage
remove_track = [] #to remove straightforward events (just control or heat and not mixed)
for i in range(1,len(data)):
	gene = data[i][0]
	event = data[i][2]
	chromosome = gene_dict[gene][-1] #-1 is the chromosome name (string)
	start,end = get_coordinates(event)
	
	result_control,perc_control = coverage_stats(control,chromosome,start,end,cov_perc)
	
	if result_control==-1:
		continue

	
	percentages.append(perc_control)
	
	if result_control == True:
		count_control+=1
	
	else:
		remove_track.append(i)
	
	
	
	
	
	


data_final = np.delete(data,remove_track,axis=0)

np.savetxt(fname+'_IRCovered.txt',data_final,fmt='%s',delimiter='\t')

print len(data_final)
print 'control',count_control

		
		
		
		
		
		
