# RYcode
Script for RY/SW/KM coding of phylogenetic data

Used for analyses in:
Braun EL & Kimball RT 2021. Data types and the phylogeny of Neoaves. Birds, 2, 1-22. 
https://doi.org/10.3390/birds2010001

This script uses a non-interleaved relaxed phylip file of DNA sequences as input and
it produces either a binary relaxed phylip file or a nexus file with the data encoded
as R and Y.

The program is straightforward to use; it will provide instructions when called without
command line options:

./recodeRY.pl 

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
