# RYcode (and simple_concat.pl)
Script for RY/SW/KM coding of phylogenetic data (a simple concatenation script, first
used in the same study as this RY-coding script, is described below).

Used for analyses in:
Braun EL & Kimball RT 2021. Data types and the phylogeny of Neoaves. Birds, 2, 1-22. 
https://doi.org/10.3390/birds2010001

This script uses a non-interleaved relaxed phylip file of DNA sequences as input and
it produces either a binary relaxed phylip file or a nexus file with the data encoded
as R and Y.

The program is straightforward to use; it will provide instructions when called without
command line options, as follow:

./recodeRY.pl 

or

perl recode.pl

```
Usage:
  $ ./recodeRY.pl <infile> <outfile> <code> <binary>
  infile  = relaxed phylip format DNA data
  outfile = relaxed phylip format 01 (or RY) data
  code    = RY (A or G = 0; C or T = 1)
            SW (G or C = 0; A or T = 1)
            KM (G or T = 0; A or C = 1)
  binary  = yes (0/1) or no (R/Y)

Setting 'binary' to 'no' will always code the datas as R and Y with the
convention 0->R and 1->Y
(i.e., code=SW and binary=no will code G and C as R and A and T as Y; this
is done to make it possible to use of PAUP* to examine base composition)
exiting...
```

The 0/1 data can be analyzed in programs like IQ-TREE and RAxML to generate estimates
of phylogeny. The RY nexus files can be analyzed in PAUP (note that the use of RY in
this context is a kludge that allows use PAUP to calculate nucleotide frequencies).

This repository also includes the simple concatenation code used for the Braun and
Kimball (2020) project.

--------------------------------------------------------------------------------
simple_concat.pl

A simple perl program to concatenate data from multiple loci in a set of separate
alignment files, all in relaxed phylip format (sequential single-line format). As
with the recodeRY.pl script, calling simple_concat.pl with no command line options
prints information about usage:

Usage:
  $ simple_concat.pl <authority> <filelist> <outfile>
  authority = list of taxa
  filelist  = list of phylip format files
  outfile   = outfile prefix
     -- concatenated nexus file: <outfile>.nex
     -- nexus sets block:        <outfile>.partitions.txt
     -- RAxML paratition file:   <outfile>.raxparts.txt
exiting...

Authority file:

The authority file is a list of the taxa to include in the final dataset, with one
taxon per line. If a taxon in the authority file is absent from an individual gene
file it will be added to the block for that locus and coded as missing data (?). If a
taxon is present in an individual gene file it will be omitted from the concatenated
output file

Filelist:

The filelist is a list of all individual files to concatenate. The files should be
in sequential relaxed phylip format, with the sequences on a single line. Taxon names
should not include whitespace. An example of this format follows:

5    42
Turkey                  AAGCTNGGGCATTTCAGGGTGAGCCCGGGCAATACAGGGTAT
Salmo_schiefermuelleri  AAGCCTTGGCAGTGCAGGGTGAGCCGTGGCCGGGCACGGTAT
H_sapiens               ACCGGTTGGCCGTTCAGGGTACAGGTTGGCCGTTCAGGGTAA
Chimp                   AAACCCTTGCCGTTACGCTTAAACCGAGGCCGGGACACTCAT
Gorilla                 AAACCCTTGCCGGTACGCTTAAACCATTGCCGGTACGCTTAA
  
PAUP* will output data in this format using the following command:

  export file=FILENAME.phy format=RelPHYLIP charsperline=all replace;
  
Outfile prefix:
  
The concatenated file will be called OUTFILE.nex and it will be an interleaved nexus
file with a data block that is readable by PAUP* (and at least some other nexus readers). 
Each individual locus file is presented separately with comments identifying the locus. 
A sets block with the boundaries for each locus is also included in the file. A separate
nexus file with the locus boundaries is also available, as is a RAxML format partition
file. 
  
If OUTFILE.nex does not import properly when you read the file in another program, check the 
partition names. This is often a source of problems. For example, some phylogenetic programs may
find a partition name like gene1.phy problematic but be fine with gene1; in that example
you could correct the file using sed 's/.phy//g' OUTFILE.nex > temp (and then renaming the
temp file OUTFILE.nex). If the problem appears to be the nexus data block, execute the file
in PAUP* and then export in a usable format.
