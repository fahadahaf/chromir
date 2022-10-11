# ChIP TF occupancy analysis
The ChIP TF occupancy analysis are done in *Python 2.7*, this is because we used SpliceGrapher and iDiffIR libraries which are written in the older Python version.

## Dependency
The following python 2.7 packages are required:  
[iDiffIR (version 0.0.1)](https://combi.cs.colostate.edu/idiffir/installation.html)  
[SpliceGrapher (version 0.2.5)](https://sourceforge.net/projects/splicegrapher/)  
[biopython (version 1.74)](https://biopython.org)   
[matplotlib (vresion 1.5.3)](https://matplotlib.org)  
[numpy (version 1.14.3)](www.numpy.org)   
[scipy (version 0.18.1)](www.scipy.org)  

##Example
To generate ChIP occupancy stats for FOXK2, run the following shell script:
```
sh run_all_nonIR_FOXK2.sh
```

For other TFs, this script needs to be modified accordingly. 
