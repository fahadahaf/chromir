# chromIR
Code repository for our work on the co-transcriptional regulation of Intron Retention using chromatin accessibility data. 

## Manuscript
Fahad Ullah, Maayan Salton, Anireddy S.N. Reddy, Asa Ben-Hur; Evidence for the role of transcription factors in the co-transcriptional regulation of intron retention, bioRxiv 2021; [https://www.biorxiv.org/content/10.1101/2021.11.18.469150v1](https://www.biorxiv.org/content/10.1101/2021.11.18.469150v1)

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
