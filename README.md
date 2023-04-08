# chromIR
Code repository for our work on the co-transcriptional regulation of Intron Retention using chromatin accessibility data. 

## Manuscript
Fahad Ullah, Saira Jabeen, Maayan Salton, Anireddy S.N. Reddy, Asa Ben-Hur; Evidence for the role of transcription factors in the co-transcriptional regulation of intron retention, Genome Biology 2023; [https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-02885-1](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-02885-1)

## Dependency
*chromir* is written in python 3. The following python packages are required:  
[fastprogress (version 0.1.21)](https://github.com/fastai/fastprogress)  
[matplotlib (vresion 3.1.3)](https://matplotlib.org)  
[numpy (version 1.17.2)](www.numpy.org)   
[pandas (version 0.25.1)](www.pandas.pydata.org)  
[pytorch (version 1.2.0)](https://pytorch.org)  
[scikit-learn (vresion 0.24)](https://scikit-learn.org/stable/)  
[scipy (version 1.4.1)](www.scipy.org)  
[statsmodels (version 0.9.0)](http://www.statsmodels.org/stable/index.html)  
[gensim (version 3.8.1)](https://radimrehurek.com/gensim/)

and for motif analysis:  
[MEME suite](http://meme-suite.org/doc/download.html)  
[WebLogo](https://weblogo.berkeley.edu)

## Usage
```
usage: main.py [-h] [-v] [-o DIRECTORY] [-m MODE] [--deskload] [-w NUMWORKERS]
               [--splitperc SPLITPERC] [--motifanalysis]
               [--scorecutoff SCORECUTOFF] [--tomtompath TOMTOMPATH]
               [--database TFDATABASE] [--annotate ANNOTATETOMTOM] [-s]
               [--tomtomdist TOMTOMDIST] [--tomtompval TOMTOMPVAL]
               [--wvpath WVPATH] [--nettype NETTYPE] [--embd]
               inputprefix hparamfile

Main chromIR script to run experiments.

positional arguments:
  inputprefix           Input file prefix for the bed/text file and the
                        corresponding fasta file (sequences).
  hparamfile            Name of the hyperparameters file to be used.

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         verbose output [default is quiet running]
  -o DIRECTORY, --outDir DIRECTORY
                        output directory
  -m MODE, --mode MODE  Mode of operation: train or test.
  --deskload            Load dataset from desk. If false, the data is
                        converted into tensors and kept in main memory (not
                        recommended for large datasets).
  -w NUMWORKERS, --numworkers NUMWORKERS
                        Number of workers used in data loader. For loading
                        from the desk, use more than 1 for faster fetching.
  --splitperc SPLITPERC
                        Pecentages of test, and validation data splits, eg. 10
                        for 10 percent data used for testing and validation.
  --motifanalysis       Analyze CNN filters for motifs and search them against
                        known TF database.
  --scorecutoff SCORECUTOFF
                        In case of binary labels, the positive probability
                        cutoff to use.
  --tomtompath TOMTOMPATH
                        Provide path to where TomTom (from MEME suite) is
                        located.
  --database TFDATABASE
                        Search CNN motifs against known TF database. Default
                        is Human CISBP TFs.
  --annotate ANNOTATETOMTOM
                        Annotate tomtom motifs. The value of this variable
                        should be path to the database file used for
                        annotation. Default is None.
  -s, --store           Store per batch CNN outpout matrices. If false, the
                        are kept in the main memory.
  --tomtomdist TOMTOMDIST
                        TomTom distance parameter (pearson, kullback, ed etc).
                        Default is pearson. See TomTom help from MEME suite.
  --tomtompval TOMTOMPVAL
                        Adjusted p-value cutoff from TomTom. Default is 0.05.
  --wvpath WVPATH       Path to where the word2vec trained embeddings are
                        located. Default is None.
  --nettype NETTYPE     Model type to use: either basset or attention. Default
                        is basset.
  --embd                Whether to use word2vec embeddings instead of one-hot
                        encoded input. Default is False. Make sure to provide
                        the appropriate hyperparameters file.
```

## Tutorial
You can skip **1** and **2** since the data has been provided.
1. To generate the data, first run the preprocessing script:
```
python preprocess.py IR_events_file.txt encode_intragenic_DHSs.bed
```
This will create the final dataset file `Labelled_Data_IR_iDiffIR_corrected.txt` in **data/** directory.  
**Note:** Check out the **preprocessing/** directory for more details on the IR events file format, encode DHSs data, and how to get the intragenic DHSs.

2. To generate the fasta file with DNA sequences for the corresponding DHSs, use [bedtools](https://bedtools.readthedocs.io/en/latest/): 
```
bedtools getfasta -fi hg19_human_reference.fa -bed data/Labelled_Data_IR_iDiffIR_corrected.txt -s -fo data/Labelled_Data_IR_iDiffIR_corrected.fa
```

3. Running *chromir* with motif analysis. Note that the output directory is stored under **results/**. 
```
python main.py -v -o Basset_noEmbd -m train --motifanalysis --nettype basset data/Labelled_Data_IR_iDiffIR_corrected data/models/Basset-noEmbds_hyperParams.txt
```
