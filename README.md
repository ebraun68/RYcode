# RYcode
Script for RY/SW/KM coding of phylogenetic data

Used for analyses in:
Braun and Kimball (2020) Data types and the phylogeny of Neoaves. Submitted.

This script uses a non-interleaved relaxed phylip file of DNA sequences as input and
it produces either a binary relaxed phylip file or a nexus file with the data encoded
as R and Y.

The program is straightforward to use; it will provide instruction when called without
command line options:

./recodeRY.pl 
#Usage:
#  $ ./recodeRY.pl <infile> <outfile> <code> <binary>
#  infile  = relaxed phylip format DNA data
#  outfile = relaxed phylip format 01 (or RY) data
#  code    = RY (A or G = 0; C or T = 1)
#            SW (G or C = 0; A or T = 1)
#            KM (G or T = 0; A or C = 1)
#  binary  = yes (0/1) or no (R/Y)

Setting 'binary' to 'no' will always code the datas as R and Y with the
convention 0->R and 1->Y
(i.e., code=SW and binary=no will code G and C as R and A and T as Y; this
is done to make it possible to use of PAUP* to examine base composition)
exiting...

This repository also includes the simple concatenation code used for the Braun and
Kimball (2020) project.
