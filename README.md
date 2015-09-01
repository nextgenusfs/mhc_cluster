###USEARCH MHC Clustering Script:###

___
####Installation:####
You can do a git clone to copy this repository:

`git clone https://github.com/nextgenusfs/mhc_cluster.git`

And then you will need to add to your $PATH or always include the entire path to the scripts at runtime.

`export PATH="/location/of/packages/mhc_cluster:$PATH"`

You will also need to install USEARCH8 - get it here: http://www.drive5.com/usearch/download.html.  One way to make the program executable and move into your path:

```
#make executable
sudo chmod +x /path/to/usearch8.0.1623_i86osx32
```

```
#create softlink to folder in $PATH
sudo ln -s /path/to/usearch8.0.1623_i86osx32 /usr/local/bin/usearch8
```

You will also need to have HMMer3 installed, simpliest way on MacOSX is with Homebrew (http://brew.sh).  And then you can tap the homebrew/science rep (https://github.com/Homebrew/homebrew-science)

```
brew tap homebrew/science
brew install hmmer
```

####Processing Ion Torrent Data:####

From the Ion Torrent Server, analyze the data using the `--disable-all-filters` BaseCaller argument. This will leave the adapters/key/barcode sequence intact. The data need to be exported as a FASTQ file, or alternatively use a 3rd party tool to convert the BAM output file to FASTQ (i.e. `bedtools bamtofastq -i <BAM> -fq <FASTQ>`). You can then de-multiplex the data as follows:

`mhc-process_reads.py --barcodes 1,5,24 --rev_comp data.fastq > data.demux.fq`

This will find Ion barcodes (1, 5, and 24) and relabel header with that information (barcodelabel=BC_5;). By default,     it will look for all 96 Ion Xpress barcodes, specifiy the barcodes you used by a comma separated list. Next the script will find and trim both the forward and reverse primers (default is MHC_DRB region: drbF & drbR). This will then save to STDOUT only the sequences that have a valid barcode sequence and contain both forward/reverse primers, and finally it will reverse-complement the sequence so it is in the proper orientation for the next step (if you used the drbF/drbR primers).  These options can be customized using: `--fwd_primer`, `--rev_primer`, etc. Type `-h` for all the available options.

####Clustering your data into OTUs:####

The next step is to quality filter the data (remove low quality reads), remove contaminating reads by filtering with a DRB exon 2 HMM model (HMMER3), trimming the reads to a set length (default=average read length), and then utilize the UPARSE OTU clustering algorithm via USEARCH8.  This is done as follows:

`mhc-OTU_cluster.py -f data.demux.fq -o output`

This script is a wrapper for hmmscan and USEARCH8, resulting in an OTU table, multi-fasta OTU file, etc. for the end result.  Here the reads are filtered with a max expected errors of 1.0, filtered with an HMM DRB2 DNA model, trim/padded to average length of those reads that pass the filter, and then run through the UPARSE pipeline (dereplicate, sort, cluster, map to OTUs, build OTU table).
    
####Dependencies####
python, biopython, USEARCH8 (accessible in PATH; alternatively you can pass in the variable `-u /path/to/usearch8` to scripts requiring USEARCH8).

