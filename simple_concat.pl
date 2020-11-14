#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;

############################################################################
# Set the global variables
############################################################################
my($progname) = $0;

my($iter);
my($jter);
my($kter);
my($lter);
my($mter);
my($nter);

if ( @ARGV != 3 ) {
	print "Usage:\n  \$ $progname <authority> <filelist> <outfile>\n";
	print "  authority = list of taxa\n";
	print "  filelist  = list of phylip format files\n";
	print "  outfile   = concatenated outfile\n";
	print "exiting...\n";
	exit;
}

my($authfile) = $ARGV[0];
my($listfile) = $ARGV[1];
my($outfile)  = $ARGV[2];

############################################################################
# Read the control file
############################################################################
print "Reading the authority file... ";
open (my $AUTHF, $authfile) or die "Could not open file $authfile for input.\n";
my @authlist = <$AUTHF>; # Read the authority file
close($AUTHF) or die "Could not close file $authfile.\n";
my($authnum) = $#authlist + 1;
print "$authnum taxon names read\n\n";

my($maxlen) = 0;
for ($iter=0; $iter<$authnum; $iter++) {
	chomp($authlist[$iter]);
	if ( length($authlist[$iter]) > $maxlen ) { $maxlen=length($authlist[$iter]); }
	print "   $authlist[$iter]\n";
}
print "Maximum taxon name length is $maxlen characters\n\n";

print "Reading the list of phylip input files... ";
open (my $LISTF, $listfile) or die "Could not open file $listfile for input.\n";
my @listlist = <$LISTF>; # Read the list file
close($LISTF) or die "Could not close file $listfile.\n";
my($listnum) = $#listlist + 1;
print "will concatenate data from $listnum files\n\n";

############################################################################
# Iterate through the files and generate concatenated files
############################################################################
my @datalist;
my($name);
my($seq);
my($namelen);

my($ntaxa);
my($nsites);
my($phynum);
my @phylist;
my($found);
my($totalsites)=0;

# First, calculate the total number of sites and write a sets block
open (my $PARTF, ">$outfile.partitions.txt") or die "Could not open file $outfile.partitions.txt for output.\n";
my($firstsite) = 1;
my($lastsite);

print $PARTF "\nBegin sets;\n";

for ($iter=0; $iter<$listnum; $iter++) {
	chomp($listlist[$iter]);
#	print "  -- file = $listlist[$iter]\n";
	
	# Read the phylip input file, calculate the total number of sites
	open (my $PHYF, $listlist[$iter]) or die "Could not open file $listlist[$iter] for input.\n";
	@phylist = <$PHYF>; # Read the phylip input file
	close($PHYF) or die "Could not close file $listlist[$iter]\n";
	$phynum = $#phylist + 1;
	($ntaxa,$nsites) = split(/\s+/, $phylist[0]);
	chomp($nsites);
	$totalsites = $totalsites + $nsites;
	
	$lastsite = $firstsite + $nsites - 1;
	print $PARTF "\tcharset $listlist[$iter] = $firstsite - $lastsite ; [ $nsites sites ]\n";
	$firstsite = $firstsite + $nsites;
}

print $PARTF "End;\n\n";
close($PARTF) or die "Could not close file $outfile.partitions.txt\n";


print "The final concatenated file will have a total of $totalsites sites\n\n";

# Then iterate through the files and output the nexus file
open (my $OUTF, ">$outfile") or die "Could not open file $outfile for output.\n";

print $OUTF "#NEXUS\n\n";
print $OUTF "Begin data;\n";
print $OUTF "\tDimensions ntax=$authnum nchar=";
print $OUTF "$totalsites";
print $OUTF ";\n";
print $OUTF "\tFormat datatype=dna missing=? gap=- interleave;\n";
print $OUTF "Matrix\n\n";

$totalsites=0; # reset $totalsites

for ($iter=0; $iter<$listnum; $iter++) {
	chomp($listlist[$iter]);
	print "  -- file = $listlist[$iter]\n";
	
	# Read the phylip input file
	open (my $INF, $listlist[$iter]) or die "Could not open file $listlist[$iter] for input.\n";
	@phylist = <$INF>; # Read the phylip input file
	close($INF) or die "Could not close file $listlist[$iter]\n";
	$phynum = $#phylist + 1;
	($ntaxa,$nsites) = split(/\s+/, $phylist[0]);
	chomp($nsites);
	$totalsites = $totalsites + $nsites;
	print $OUTF "\[ Data block $listlist[$iter] -- $nsites sites -- $totalsites total sites \]\n";
	
	for ($jter=0; $jter<$authnum; $jter++) {
		print $OUTF "$authlist[$jter]  ";
		$namelen = $maxlen - length($authlist[$jter]);
		for ($kter=0; $kter<$namelen; $kter++) { print $OUTF " "; }
		$found = 0;
		for ($kter=1; $kter<$phynum; $kter++) {
			($name,$seq) = split(/\s+/, $phylist[$kter]);
			if ( $name eq $authlist[$jter] ) {
				chomp($seq);
				print $OUTF "$seq\n";
				$found = 1;
			}
		}
		if ( $found == 0 ) {
			for ($lter=0; $lter<$nsites; $lter++) { print $OUTF "?"; }
			print $OUTF "\n";
		}
	} # end for ($jter=0; $jter<$authnum ...
	print $OUTF "\n\n";

}

print $OUTF "\t;\nEnd;\n\n";

close($OUTF) or die "Could not close file $outfile.\n";

# Append the sets block to the nexus file
system("cat $outfile.partitions.txt >> $outfile");

print "\nFile concatenation complete...\n";

exit;
