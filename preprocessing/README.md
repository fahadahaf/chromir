# Data Preparation
## Human DHSs data
To get the DHS peaks, we followed pipeline mentioned in the Basset paper by Kelley et al.:
- Open the following notebook, part of Basset tutorial: https://github.com/davek44/Basset/blob/master/tutorials/prepare_compendium.ipynb
- Run the first command (that uses the preprocess_features.py script) to get all the Encode DHSs  

## Human Intragenic DHSs
To get intragenic only DHSs, run the provided script:
```
python get_intragenic_DHSs.py Basset_all_DHSs.bed Human_hg19_annotations.gtf Basset_all_DHSs_intragenic.bed
```
Arguments to the script are:  
- All Encode human DHSs we got in first step (.bed file)
- Human hg19 annotations (.gtf file)
- Output filename with intragenic only DHSs

## Human IR events data
We used Splicegrapher to get list of IR events from human genome. Since Splicegrapher is written in Python 2.7, you can use any other tool that will give you list of IR events, using either just human genome annotations, or both annotations and evidence from RNA-seq data.  
Note that the IR events file should have the following column format:  
```
Gene_ID, Pre_Exon, Intron, Post_Exon, Strand, Chromosome, First_Intron
```
- Gene_ID: ID of the gene in which the IR event is reported. [required]
- Pre_Exon: Genomic coordinates (start, end) of the exon preceeding the intron. Note that these are absolute genomic coordinates; the Strand field determines if this is a 3' exon or a 5' exon. [required]
- Intron: Genomic coordinates (start, end) of the intron. [required]
- Post_Exon: Genomic coordinates (start, end) of the exon following the intron. Note that these are absolute genomic coordinates; the Strand field determines if this is a 5' exon or a 3' exon. [required]
- Strand: Strand of the gene, "+" or "-". [required]
- Chromosome: Chromosome of the gene. [required]
- First_Intron: Whether the intron is first intron of the gene, values 1 or 0. [optional]