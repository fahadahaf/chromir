from iDiffIR.IntronModel import *
from iDiffIR.Stat import *
from SpliceGrapher.formats.fasta import *
from SpliceGrapher.shared.utils      import *
from SpliceGrapher.formats.loader import loadGeneModels
geneModel = loadGeneModels('ensGene.gtf',verbose = True)
geneRecords = makeModels(geneModel,'extract_IR_test',verbose=True,exonic=False,procs=8)
