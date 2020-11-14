#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;

############################################################################
# Initialize variables
############################################################################
my($progname) = $0;
my($iter);
my($jter);
my($kter);
my($lter);
my($mter);
my($nter);

my($tempch);
my($tempvar);
my @temparray;

if ( @ARGV != 4 ) {
	print "Usage:\n  \$ $progname <infile> <outfile> <code> <binary>\n";
	print "  infile  = relaxed phylip format DNA data\n";
	print "  outfile = relaxed phylip format 01 (or RY) data\n";
	print "  code    = RY (A or G = 0; C or T = 1)\n";
	print "            SW (G or C = 0; A or T = 1)\n";
	print "            KM (G or T = 0; A or C = 1)\n";
	print "  binary  = yes (0/1) or no (R/Y)\n";
	print "\nSetting 'binary' to 'no' will always code the datas as R and Y with the\n";
	print "convention 0->R and 1->Y\n";
	print "(i.e., code=SW and binary=no will code G and C as R and A and T as Y; this\n";
	print "is done to make it possible to use of PAUP* to examine base composition)\n";
	print "exiting...\n";
	exit;
}

my($infile)=$ARGV[0];
my($outfile)=$ARGV[1];
my($code)=$ARGV[2];
my($usebinary)=$ARGV[3];

############################################################################
# Read the input files
############################################################################
print "Reading data from file $infile";
print "... ";

open (my $INF, $infile) or die "Could not open file $infile for input.\n";
my @seqinfo = <$INF>; # Read the DNA alignment file
close($INF) or die "Could not close file $infile\n";
my($seqinfonum) = $#seqinfo + 1;

print "done\n";

my($ntax);
my($nchar);

chomp($seqinfo[0]);
($ntax,$nchar) = split(/\s+/, $seqinfo[0]);



############################################################################
# Set up the coding
############################################################################
my($zero) = "0";
my($one)  = "1";
if ( lc($usebinary) eq "no" ) {
	$zero = "R";
	$one  = "Y";
}

############################################################################
# Convert to a file with the appropriate binary coding
############################################################################
my($name);
my($seq);
my($base);

# Find longest name
my($maxlen) = 10;
for ($iter=0; $iter<$ntax; $iter++) {
	chomp($seqinfo[$iter+1]);
	($name,$seq) = split(/\s+/, $seqinfo[$iter+1]);
	if ( length($name) > $maxlen ) { $maxlen = length($name); }
}

print "Writing ";
if ( lc($usebinary) eq "no" ) { print "RY symbol"; }
else { print "binary"; }
print " data to $outfile";

open (my $OUTF, ">$outfile") or die "Could not open file $outfile for output.\n";

if ( lc($usebinary) eq "no" ) {
	print $OUTF "#NEXUS\n\n";
	print $OUTF "Begin data;\n";
	print $OUTF "\tDimensions ntax=$ntax nchar=$nchar;\n";
	print $OUTF "\tFormat datatype=nucleotide gap=- missing=? matchchar=.;\n";
	print $OUTF "Matrix\n";
}
else { print $OUTF "$ntax $nchar\n"; }

for ($iter=0; $iter<$ntax; $iter++) {
	if ( $iter % 10 == 0 ) { print "."; }
#	chomp($seqinfo[$iter+1]);
	($name,$seq) = split(/\s+/, $seqinfo[$iter+1]);
	print $OUTF "$name  ";
	$tempvar = $maxlen - length($name) + 1;
	for ($jter=0; $jter<$tempvar; $jter++) {
		print $OUTF " ";
	}
	for ($jter=0; $jter<$nchar; $jter++) {
		$base = substr($seq,$jter,1);
		if ( $base eq "-" || $base eq "?" ) { print $OUTF "$base"; }
		elsif ( $code eq "SW" || $code eq "sw" ) {
			$tempch = "?";
			if ( uc($base) eq "G" || uc($base) eq "C" ) { $tempch = "$zero"; }
			elsif ( uc($base) eq "A" || uc($base) eq "T" ) { $tempch = "$one"; }
			elsif ( uc($base) eq "S" ) { $tempch = "$zero"; }
			elsif ( uc($base) eq "W" ) { $tempch = "$one"; }
			print $OUTF "$tempch";
		}
		elsif ( $code eq "KM" || $code eq "km" ) {
			$tempch = "?";
			if ( uc($base) eq "G" || uc($base) eq "T" ) { $tempch = "$zero"; }
			elsif ( uc($base) eq "A" || uc($base) eq "C" ) { $tempch = "$one"; }
			elsif ( uc($base) eq "K" ) { $tempch = "$zero"; }
			elsif ( uc($base) eq "M" ) { $tempch = "$one"; }
			print $OUTF "$tempch";
		}
		else { # NOTE - program will default to RY coding
			$tempch = "?";
			if ( uc($base) eq "A" || uc($base) eq "G" ) { $tempch = "$zero"; }
			elsif ( uc($base) eq "C" || uc($base) eq "T" ) { $tempch = "$one"; }
			elsif ( uc($base) eq "R" ) { $tempch = "$zero"; }
			elsif ( uc($base) eq "Y" ) { $tempch = "$one"; }
			print $OUTF "$tempch";
		}
	}
	print $OUTF "\n";
}

if ( lc($usebinary) eq "no" ) {
	print $OUTF "\t;\nEnd;\n\n";
}

close($OUTF) or die "Could not close file $outfile\n";

print " done\n\n";

