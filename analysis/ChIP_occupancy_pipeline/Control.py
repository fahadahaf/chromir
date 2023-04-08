import sys
#sys.path.append('/s/chopin/c/proj/protfun/arch/x86_64/lib/python')

from iDiffIR.IntronModel import *
import math

from SpliceGrapher.formats.fasta import *
import os,sys,numpy
from SpliceGrapher.shared.utils import *
from SpliceGrapher.formats.loader import loadGeneModels


##################################### When DHS is bigger than pre-intron-post region#################
#        ---  ---- ----    pre intron post
#    ______              Yes No No
#  
#    ___________         No Yes No
# 
#    ________________    No No Yes
#
#    ___________________ No No No
#
#         _________________ Yes No No
#
#               ___________ No Yes No
#
#                    ______ No No Yes
#####################################################################################################
#################################################
import numpy as np

#countx = 0
#county = 0

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
	return 'No' if f_pre=='none' else 'Yes','No' if f_in=='none' else 'Yes','No' if f_post=='none' else 'Yes'     
        

def resolve_parts(pre_start,pre_end,intron_start,intron_end,post_start,post_end,DHS_start,DHS_end,pre_exon,intron,post_exon):
	if DHS_end < pre_end and DHS_start < pre_start and DHS_end > pre_start:
		return 'Yes','No','No'
	elif DHS_end < intron_end and DHS_start < pre_start and DHS_end > intron_start:
		return 'Yes','Yes','No'
	elif DHS_end < post_end and DHS_start < pre_start and DHS_end > post_start:
		return 'Yes','Yes','Yes'
	elif DHS_start > pre_start and DHS_end > post_end and DHS_start < pre_end:
		return 'Yes','Yes','Yes'
	elif DHS_start > intron_start and DHS_end > post_end and DHS_start < intron_end:
		return 'No','Yes','Yes'
	elif DHS_start > post_start and DHS_end > post_end and DHS_start < post_end:
		return 'No','No','Yes'
	elif DHS_start < pre_start and DHS_end > post_end:
		return 'Yes','Yes','Yes'  
	else:
		return pre_exon,intron,post_exon
	       
        

print "Loading gene model..."

geneModel = loadGeneModels('../ensGene.gtf',verbose=True)

Diff_IR = np.loadtxt('Union_IR_Events.txt',dtype=str,delimiter='\t')

TF_name = sys.argv[1] #MAZ, EGR1 etc
celltype = sys.argv[2] #K562 in this case

#DHS_List = np.loadtxt('GSE34318_RAW/wdBuds_Pooled/F-Seq/filtered_DHSs.csv',dtype=str,delimiter=',')
DHS_List = np.loadtxt('../'+TF_name+'_mutliCellLines_peaks.bed',dtype=str,delimiter='\t')

tracker = [['GeneID','Event','is_within_gene','start_in_gene','end_in_gene','gene_Coord','in_5-Exon','in_Intron','in_3-Exon','TF_Coordinates','Pre_Exon','Post_Exon','TF_p-val','TF_Sig-val','TF_peak-val','gene_Strand']]

limit = len(Diff_IR)

for i in range(1, limit): #0 is header so starting from 1
	geneID = Diff_IR[i][0]
	gene = geneModel.getGeneByName(geneID)
	chrNo = gene.chromosome.lower()
	strand = Diff_IR[i][-1]
	
	pre_exon = Diff_IR[i][1]
	intron = Diff_IR[i][2]
	post_exon = Diff_IR[i][3]
	
	for j in range(0,len(DHS_List)):
		if celltype not in DHS_List[j][-1]: #TF CHIP peak should be in the current celltype
			continue
			
		chrDHS = DHS_List[j][0]#.lower()
		dhs_start = DHS_List[j][1]
		dhs_end = DHS_List[j][2]
		DHS = '('+dhs_start+', '+dhs_end+')'
		geneCoord = '('+str(gene.minpos)+', '+str(gene.maxpos)+')'
		if chrDHS == chrNo:
			gene_check = individual_overalap_check(get_coordinates(DHS),get_coordinates(geneCoord))
			if gene_check=='none':
				continue
			pre,intr,post = overlap_check(DHS,pre_exon,intron,post_exon)
			if pre=='No' and intr=='No' and post=='No':
				continue
		
			within_gene = 'No'
			start_in_gene = 'No'
			end_in_gene = 'No'
		
			if gene_check=='all':
				within_gene = 'Yes'
				start_in_gene = 'Yes'
				end_in_gene = 'Yes'
			elif gene_check=='start':
				start_in_gene = 'Yes'
			elif gene_check=='end':
				end_in_gene = 'Yes'
			elif gene_check=='across':
				within_gene = 'No*'
				start_in_gene = 'No*'
				end_in_gene = 'No*'
			if strand == '+':
				tracker.append([geneID,intron,within_gene,start_in_gene,end_in_gene,geneCoord,pre,intr,post,DHS,pre_exon,post_exon,'-','-','-',strand])
			else:
				tracker.append([geneID,intron,within_gene,start_in_gene,end_in_gene,geneCoord,post,intr,pre,DHS,pre_exon,post_exon,'-','-','-',strand])
		
		#print gene_check,geneCoord,pre,intr,post,DHS,pre_exon,intron,post_exon
		
	
	#print 'yes'
	if i%100==0:
		print i,chrDHS,chrNo

np.savetxt(TF_name+'_in_IR.txt',tracker,fmt='%s',delimiter='\t')











#bestA = 2
#count = 0
#countx = 0
#county = 0

#count_len = 0
#len_DHS = 150

#tracker = [['GeneID','idiffIR_pval','idiffIR_adj-pval','idiffIR_fold-change','is_within_Gene','start_in_gene','end_in_gene','in_pre_exon','in_intron','in_post_exon','hotspot_pval','hotspot_sigVal','gene_coordinates','idiffIR_pre-exon','idiffIR_intron','idiffIR_post-exon','hotspot_coordinates','gene_strand']]

#for i in range(1,limit):
    #geneID = Diff_IR[i][0]
    ##if int(Diff_IR[i][-2]) < bestA:
    ##    continue
    ##else:
    #chrNo = 'chr'+geneID[2]
    #gene = geneModel.getGene(chrNo,geneID)
        
    #for j in range(1,len(DHS_List)):   
        #within_gene = 'No'
        #start_gene = 'No'
        #end_gene = 'No'
        #pre_exon = 'No'
        #in_intron = 'No'
        #post_exon = 'No'  
        #if DHS_List[j][-2] == '-': #'-' before, make '+' for testing all DHSs (differential)
            #continue
        #else:
            #chrNo_IR = DHS_List[j][0]
            #if chrNo_IR == chrNo:
                ##print 'yes'
                #start = int(DHS_List[j][1])
                #end = int(DHS_List[j][2])
                
                #ans = is_within_gene(gene,start,end)
                #if ans=='all':
                    #count += 1
                    #within_gene = 'Yes'
                #elif ans=='start':
					#countx += 1
					#start_gene = 'Yes'
                #elif ans=='end':
					#county += 1
					#end_gene = 'Yes'
                #if ans!='no':
					#start_pre,end_pre = get_coordinates(Diff_IR[i][1])
					#pre_exon = check_part(start_pre,end_pre,start,end)
					#start_int,end_int = get_coordinates(Diff_IR[i][2])
					#in_intron = check_part(start_int,end_int,start,end)
					#start_post,end_post = get_coordinates(Diff_IR[i][3])
					#post_exon = check_part(start_post,end_post,start,end)
					#summ = determine_length(Diff_IR[i][1])+determine_length(Diff_IR[i][2])+determine_length(Diff_IR[i][3]) 
					#gene_coord = '('+str(gene.minpos)+','+str(gene.maxpos)+')'
					#hotspot_coord = '('+str(DHS_List[j][1])+','+str(DHS_List[j][2])+')'
					#if summ <= len_DHS: #Length of intron and flanking exons in less than the length of the DHS
						#count_len += 1
						#print 'when DHS site is bigger'
						#print count_len,pre_exon,in_intron,post_exon
						#print Diff_IR[i][1],Diff_IR[i][2],Diff_IR[i][3]
						#print hotspot_coord
						#pre_exon,in_intron,post_exon = resolve_parts(start_pre,end_pre,start_int,end_int,start_post,end_post,start,end,pre_exon,in_intron,post_exon)
						#print pre_exon,in_intron,post_exon
					
					#if pre_exon == 'Yes' and in_intron=='No' and post_exon=='Yes':
						#print 'when Yes,No,Yes'
						#in_intron = 'Yes'
						#print Diff_IR[i][1],Diff_IR[i][2],Diff_IR[i][3]
						#print hotspot_coord
					
					#tracker.append([Diff_IR[i][0],'-','-','-',within_gene,start_gene,end_gene,pre_exon,in_intron,post_exon,DHS_List[j][7],DHS_List[j][6],gene_coord,Diff_IR[i][1],Diff_IR[i][2],Diff_IR[i][3],hotspot_coord,gene.strand])
					##print count,ans
					##print Diff_IR[i][0:3],Diff_IR[i][4:7]
					##print DHS_List[j]
					##print '\n'

#np.savetxt('diff_only_tracker-Control_allbestA.txt',tracker,fmt='%s',delimiter='\t')
                
                
                
                
                
