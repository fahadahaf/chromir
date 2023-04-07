import sys
from SpliceGrapher.formats.loader import loadGeneModels
import numpy as np

def check_overlap(element,DHS):
        sE = element[0]
        eE = element[1]
        sD = DHS[0]
        eD = DHS[1]
        return (sE <= eD) and (eE >= sD)

def isWithin(element,gene):
        if gene.minpos<element[0] and gene.maxpos>element[1]:
                return True
        else:
                return False



bed_file = sys.argv[1] #DHS BED file
gtf_file = sys.argv[2] #gene model file
gmodel = loadGeneModels(gtf_file, verbose=True)
out_file = sys.argv[3] #output file name
#bed_data = np.loadtxt(bed_file,dtype=str)
allGenes = gmodel.getAllGenes()
count = 0
count_found = 0
with open(bed_file,'r') as f, open(out_file,'w') as g:
        for line in f:
                chrom,start,end = line.strip().split('\t')[:3]
                start = int(start)
                end = int(end)
                for gene in allGenes:
                        if chrom == gene.chromosome:
                                result = isWithin((start,end),gene)
                                if result == True:
                                        line = line.strip().split('\t')
                                        line = line[:len(line)-1]+[gene.id]+[line[-1]]
                                        line = '\t'.join(line)+'\n'
                                        g.write(line)
                                        #print chrom,gene.chromosome,gene.id,start,end,gene.minpos,gene.maxpos                                                                                                     
                                        count_found += 1
                count += 1
                if count%1000 == 0:
                        print "Done: ",count, count_found
                        #break

allGenes = gmodel.getAllGenes()