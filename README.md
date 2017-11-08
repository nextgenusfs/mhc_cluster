### MHC DRB amplicon Script:

___
#### Installation:
You can do a git clone to copy this repository:

`git clone https://github.com/nextgenusfs/mhc_cluster.git`

And then you will need to add to your $PATH or always include the entire path to the scripts at runtime.

`export PATH="/location/of/packages/mhc_cluster:$PATH"`

You will also need to install USEARCH9  - get it here: http://www.drive5.com/usearch/download.html.  One way to make the program executable and move into your path:

```
#make executable
sudo chmod +x /path/to/usearch9.2.64_i86osx32
```

```
#create softlink to folder in $PATH
sudo ln -s /path/to/usearch9.2.64_i86osx32 /usr/local/bin/usearch9
```
You will also need VSEARCH, which you can download with homebrew/bioconda or directly from GitHub repo.
```
#install with homebrew
brew tap homebrew/science
brew install vsearch
```

You will also need to have HMMer3 installed, simpliest way on MacOSX is with Homebrew (http://brew.sh).  And then you can tap the homebrew/science rep (https://github.com/Homebrew/homebrew-science)

```
brew tap homebrew/science
brew install hmmer
```
   
#### Dependencies
python, biopython, VSEARCH, AMPtk, USEARCH9 (accessible in PATH; alternatively you can pass in the variable `-u /path/to/usearch8` to scripts requiring USEARCH8).

