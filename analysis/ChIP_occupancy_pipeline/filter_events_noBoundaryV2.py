import numpy as np
import sys

dfile = sys.argv[1]

data = np.loadtxt(dfile,delimiter='\t',dtype=str)

dataFinal = [data[0].tolist()]

for entry in data[1:]:
	if entry[3] == 'No' and entry[-1] == '+': #entry[3] is "start_in_gene", so we are dropping events with a CHIP peak overlapping them but not starting in gene.
		continue
	elif entry[4] == 'No' and entry[-1] == '-': #this is for negative stranded genes
		continue
	dataFinal.append(entry.tolist())
	
doutfile = dfile.split('intronCoverage')[0] + 'StartInGene_intronCoverage.txt'

np.savetxt(doutfile,dataFinal,fmt='%s',delimiter='\t')
