#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;

############################################################################
# recodeRY.pl version 1.1
#
# Takes a relaxed phylip DNA alignment as input and converts that alignment
# to one of three formats:
#
#   1. Binary (0/1) coding (relaxed phylip format outfile)
#         - can be used for RY, SW, or KM coding
#              RY (A or G = 0; C or T = 1)
#              SW (G or C = 0; A or T = 1)
#              KM (G or T = 0; A or C = 1)
#         - typical command:
#              perl recodeRY.pl inputfile.phy outputfile.phy RY binary
#
#      NOTE: binary coding is the default so any value passed as the fourth
#            item on the command line will result in binary coding
#
#   2. RY coding (all positions) (nexus format nucleotide outfile)
#         - can be used for RY, SW, or KM coding
#              RY (A or G = R; C or T = Y)
#              SW (G or C = R; A or T = Y)
#              KM (G or T = R; A or C = Y)
#         - typical command:
#              perl recodeRY.pl inputfile.phy outputfile.phy RY RY
#
#      NOTE: combining 'SW RY' or 'KM RY' as the third and fourth items 
#            passed on the command line results in the use of incorrect
#            IUPAC codes for the nucleotides. This option is included to
#            allow the use of PAUP* to calculate base frequencies.
#
#   3. Third position RY coding (nexus format nucleotide outfile)
#         - limited to RY coding (SW or KM ignored, but must be passed)
#         - every third base is changed to R or Y (or ?)
#         - typical command:
#              perl recodeRY.pl inputfile.phy outputfile.phy RY third
#
#      NOTE: the input file is assumed to be completely coding data with
#            the first position of the alignment corresponding to a first
#            codon positions. Any non-coding sequence or sequences that
#            would disrupt the reading frame (e.g., frameshifts) should be
#            removed from the input alignment before using this program.
############################################################################

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
	print "Usage:\n  \$ $progname <infile> <outfile> <code> <mode>\n";
	print "  infile  = relaxed phylip format DNA data\n";
	print "  outfile = relaxed phylip format 01 data or nexus format RY data\n";
	print "  code    = RY (A or G = 0; C or T = 1)\n";
	print "            SW (G or C = 0; A or T = 1)\n";
	print "            KM (G or T = 0; A or C = 1)\n";
	print "  mode    = binary (0/1), RY (R/Y), or third (every 3rd base R/Y)\n";
	print "\nSetting 'mode' to 'RY' will always code the datas as R and Y with the\n";
	print "following convention:\n     0->R and 1->Y\n";
	print "(i.e., code=SW and mode=RY will code G and C as R and A and T as Y; this is\n";
	print "done to make it possible to use PAUP* to examine base composition)\n";
	print "\nSetting 'mode' to 'third' will recode every third position. The alignment\n";
	print "is assumed to be exclusively coding, in frame, and start on a first position.\n";
	print "The third base positions will be coded as R/Y, not 0/1.\n";
	print "\nexiting...\n";
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
if ( lc($usebinary) eq "ry" || lc($usebinary) eq "third" ) {
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
if ( lc($usebinary) eq "ry" ) { print "RY symbol"; }
elsif ( lc($usebinary) eq "third" ) { print "third position RY symbol"; }
else { print "binary"; }
print " data to $outfile";

open (my $OUTF, ">$outfile") or die "Could not open file $outfile for output.\n";

if ( lc($usebinary) eq "ry" || lc($usebinary) eq "third" ) {
	print $OUTF "#NEXUS\n\n";
	if ( lc($usebinary) eq "third" ) {
		print $OUTF "[ Matrix has every third nucleotide subjected to RY coding ]\n";
	}
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
	if ( lc($usebinary) eq "third" ) {
		for ($jter=0; $jter<$nchar; $jter++) {
			$base = substr($seq,$jter,1);
			if ( ($jter+1) % 3 == 0 ) {
				if ( $base eq "-" || $base eq "?" ) { print $OUTF "$base"; }
				else {
					# NOTE - program is limited to RY coding for third positions
					#      - program assumes the first position of the alignment is a 
					#        first codon position and that all sites in the alignment
					#        are coding and in frame after that
					$tempch = "?";
					if ( uc($base) eq "A" || uc($base) eq "G" ) { $tempch = "$zero"; }
					elsif ( uc($base) eq "C" || uc($base) eq "T" ) { $tempch = "$one"; }
					elsif ( uc($base) eq "R" ) { $tempch = "$zero"; }
					elsif ( uc($base) eq "Y" ) { $tempch = "$one"; }
					print $OUTF "$tempch";
				}
			}
			else {
				print $OUTF "$base";
			}
		}
	} # end of if ( lc($usebinary) eq "third"...
	else { # recode all nucleotides
		
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
	} # end of "else" for "recode all nucleotides"
	
	print $OUTF "\n";
}

if ( lc($usebinary) eq "ry" || lc($usebinary) eq "third" ) {
	print $OUTF "\t;\nEnd;\n\n";
}

close($OUTF) or die "Could not close file $outfile\n";

print " done\n\n";

