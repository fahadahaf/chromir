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
	#print count,length
	perc_calc = (float(count)/length) * 100
	#print perc_calc,count,length
	if perc_calc <= perc:
		#print count,length,lst
		return True,perc_calc
	else:
		return False,perc_calc

def coverage_stats(bamfile,chromosome,start,end,cov_perc):
	exp_array = getDepthsFromBam(bamfile,chromosome,start,end)[0] #[0] is to get the actual array
	if len(exp_array)==0:
		return -1,-1
	result,perc = check_coverage(exp_array,cov_perc)
	return result,perc

def intron_coverage_level_check(intron,preExon,postExon,bamfile,chromosome,factor):
	start_in,end_in = get_coordinates(intron)
	int_mean = getDepthsFromBam(bamfile,chromosome,start_in,end_in)[0].mean()
	
	start_pre,end_pre = get_coordinates(preExon)
	pre_cov = getDepthsFromBam(bamfile,chromosome,start_pre,end_pre)[0]
	
	start_post,end_post = get_coordinates(postExon)
	post_cov = getDepthsFromBam(bamfile,chromosome,start_post,end_post)[0]
	
	flanking_mean = np.concatenate((pre_cov,post_cov)).mean()
	#print flanking_mean,int_mean,factor
	if int_mean > flanking_mean/factor:
		return True
	else:
		return False

#usage filter_IR_on_expression.py file_name
fname = sys.argv[1]
data = np.loadtxt(fname,dtype=str,delimiter='\t')
cov_perc = 30 # per cent of coverage required (X% of locations should have at least one or more reads)
reads_thresh = 1 #float(fname.split('_')[4])/2 #minimum avg. reads that IR should have (at least half of the no. of reads in overall depth of the gene)
fname = fname.split('.txt')[0]

percentages = [] #to see distribution of percentages of IR
count_control = 0 #counting track of IR events in control qualifying the minimum coverage
remove_track = [] #to remove straightforward events (just control or heat and not mixed)
for i in range(1,len(data)):
	gene = data[i][0]
	event = data[i][1]
	pre_exon = data[i][10]
	post_exon = data[i][11]
	if gene not in gene_dict.keys():
		continue #to avoid the key not found error
	chromosome = gene_dict[gene][-1] #-1 is the chromosome name (string)
	start,end = get_coordinates(event)
	
#######################For Just Control Events#####################	
	result,perc = coverage_stats(control,chromosome,start,end,cov_perc)
	if result==-1:
		continue
	percentages.append(perc)
	if result == True:
		count_control+=1
		########################Below check is for avg. coverage of intron in comparison to its flanking exons. I am not using it for now##############
		#result_2 = intron_coverage_level_check(event,pre_exon,post_exon,control,chromosome,reads_thresh)
		#if result_2 == True: #second check for IR mean
		#	count_control+=1
		#else:
		#	remove_track.append(i)
		###############################################################################################################################################
	else:
		remove_track.append(i)
###################################################################

data_final = np.delete(data,remove_track,axis=0)
print len(data_final)
np.savetxt(fname+'_intronCoverage.txt',data_final,fmt='%s',delimiter='\t')
		
print 'control: ', count_control,len(data)
		
###############Testing the final data###################
#percent_track = []
#for entry in data_final[1:]:
	#gene = entry[0]
	#event = entry[1]
	#chromosome = gene_dict[gene][-1] #-1 is the chromosome name (string)
	#start,end = get_coordinates(event)
	#if entry[2]=='Yes':
		#result,perc = coverage_stats(control,chromosome,start,end,cov_perc)
		#percent_track.append(perc)
		#if perc<100.0:
			#print entry
	#if entry[3]=='Yes':
		#result,perc = coverage_stats(heat,chromosome,start,end,cov_perc)
		#percent_track.append(perc)
		#if perc<100.0:
			#print entry
		
		
		
		
		
		
		
		
