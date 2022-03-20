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

my($tempvar);
my($tempstr);
my @temparray;

if ( @ARGV != 3 ) {
	print "Usage:\n  \$ $progname <filelist> <authority> <outfile>\n";
	print "  filelist  = list of phylip format files\n";
	print "  authority = list of taxa\n";
	print "     -- use --fromfiles to generate the authority file\n";
	print "     -- use --getnamesonly to generate the authority file and exit\n";
	print "     -- ontherwise enter the name of a file with taxon names\n";
	print "  outfile   = outfile prefix\n";
	print "     -- concatenated nexus file: <outfile>.nex\n";
	print "     -- nexus sets block:        <outfile>.partitions.txt\n";
	print "     -- RAxML partition file:    <outfile>.raxparts.txt\n";
	print "     -- Complete taxon list:     <outfile>.taxonlist.txt\n";
	print "NOTE:\n";
	print "  The filelist, authority file, and outfile prefix can be passed to\n";
	print "  program in any order if you run the program using the following flags:\n";
	print "     -l for filelist\n";
	print "     -a for authority file\n";
	print "     -o for outfile prefix\n";
	print "  Like this:\n";
	print "  \$ $progname -a=<authority> -l=<filelist> -o=<outfile>\n";
	print "exiting...\n";
	exit;
}

# Read in default order, then check for flags
my($listfile) = $ARGV[0];
my($authfile) = $ARGV[1];
my($outfile)  = $ARGV[2];

for ($iter=0; $iter<3; $iter++) {
	$tempstr = substr($ARGV[$iter],0,3);
	if ( $tempstr eq "-l=") { $listfile = substr($ARGV[$iter],3); }
	if ( $tempstr eq "-a=") { $authfile = substr($ARGV[$iter],3); }
	if ( $tempstr eq "-o=") { $outfile =  substr($ARGV[$iter],3); }
	if ( lc($ARGV[$iter]) eq "--fromfiles" )  { $authfile = "--fromfiles"; }
	if ( lc($ARGV[$iter]) eq "--getnames" )  { $authfile = "--getnames"; }
	if ( lc($ARGV[$iter]) eq "--getnamesonly" )  { $authfile = "--getnamesonly"; }
}

# Echo information about the input files
print "\nRunning $progname\n\n";
if ( $authfile eq "--getnamesonly" ) { print "Identify taxon names "; }
else { print "Concatenate sequences "; }
print "in the relaxed phylip files listed in $listfile\n";
if ( $authfile eq "--fromfiles" || $authfile eq "--fromfiles" ) {
	print "Generate taxon authority file from the taxon names in the input files\n";
}
elsif ( $authfile eq "--getnamesonly" ) {
	print "Generate taxon authority file from the taxon names in the input files and EXIT...\n";
}
else {
print "Use taxon authority file $authfile to organize taxa in concatenated file\n";

}
print "Write data to output files with prefix $outfile";
if ( $authfile eq "--getnamesonly" ) {
print " (authority file is $outfile.taxonlist.txt)\n\n";
}
else {
print " (main output file is $outfile.nex)\n\n";
}

############################################################################
# Read the list of sequence alignment files
############################################################################
print "Reading the list of relaxed phylip input files... ";
open (my $LISTF, $listfile) or die "Could not open file $listfile for input.\n";
my @listlist = <$LISTF>; # Read the list file
close($LISTF) or die "Could not close file $listfile.\n";
my($listnum) = $#listlist + 1;
print "will concatenate data from $listnum files\n\n";

############################################################################
# Iterate through the files and generate sets block
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

my @completenameset;
my($completenamenum)=0;

# First, calculate the total number of sites and write a sets block
open (my $PARTF, ">$outfile.partitions.txt") or die "Could not open file $outfile.partitions.txt for output.\n";
my($firstsite) = 1;
my($lastsite);

# Also write a RAxML partitions block in parallel with the nexus sets block
open (my $RAXF, ">$outfile.raxparts.txt") or die "Could not open file $outfile.raxparts.txt for output.\n";

print $PARTF "#NEXUS\n\nBegin sets;\n";

for ($iter=0; $iter<$listnum; $iter++) {
	chomp($listlist[$iter]);
	print "  -- file = $listlist[$iter]\n";
	
	# Read the phylip input file, calculate the total number of sites
	open (my $PHYF, $listlist[$iter]) or die "Could not open file $listlist[$iter] for input.\n";
	@phylist = <$PHYF>; # Read the phylip input file
	close($PHYF) or die "Could not close file $listlist[$iter]\n";
	$phynum = $#phylist + 1;
	chomp($phylist[0]);
	(@temparray) = split(/\s+/, $phylist[0]);
	$tempvar = @temparray; # check whether ntax and nsites is preceded by a space
	$ntaxa = $temparray[$tempvar-2];
	$nsites = $temparray[$tempvar-1];
#	print "READING phylip file $listlist[$iter] -- $tempvar -- $ntaxa -- $nsites\n";
	$totalsites = $totalsites + $nsites;
	
	$lastsite = $firstsite + $nsites - 1;
	print $PARTF "\tcharset $listlist[$iter] = $firstsite - $lastsite ; [ $nsites sites ]\n";
	$firstsite = $firstsite + $nsites;
	# Print the RAxML partition file
	print $RAXF "DNA, $listlist[$iter]";
	print $RAXF "=";
	print $RAXF "$firstsite";
	print $RAXF "-";
	print $RAXF "$lastsite\n";
	
	# Check the taxon names
	if ( $iter == 0 ) {
		for ($jter=0; $jter<$ntaxa; $jter++) {
			(@temparray) = split(/\s+/, $phylist[$jter+1]);
			$completenameset[$jter] = $temparray[0];
			$completenamenum++;
		}
	}
	else {
		for ($jter=0; $jter<$ntaxa; $jter++) {
			(@temparray) = split(/\s+/, $phylist[$jter+1]);
			$found = 0;
			for ($kter=0; $kter<$completenamenum; $kter++) {
				if ( $temparray[0] eq $completenameset[$kter] ) { $found = 1; }
			}
			if ( $found == 0 ) {
				push(@completenameset, $temparray[0]);
				$completenamenum++;
			}
		}
	}

}

print $PARTF "End;\n\n";
close($PARTF) or die "Could not close file $outfile.partitions.txt\n";
close($RAXF) or die "Could not close file $outfile.raxparts.txt\n";

print "The final concatenated file will have a total of $totalsites sites\n";
print "A total of $completenamenum unique taxon names were found in the input files\n\n";

open (my $NAMEF, ">$outfile.taxonlist.txt") or die "Could not open file $outfile.taxonlist.txt for output.\n";
for ($iter=0; $iter<$completenamenum; $iter++) {
	print $NAMEF "$completenameset[$iter]\n";
}
close($NAMEF) or die "Could not close file $outfile.taxonlist.txt\n";

if ( $authfile eq "--fromfiles" || $authfile eq "--getnames" ) { 
	$authfile = "$outfile.taxonlist.txt"; 
}
elsif ( $authfile eq "--getnamesonly" ) { 
	print "Run in --getnamesonly mode. Taxon names written to $outfile.taxonlist.txt\n";
	print "exiting...\n";
	exit; 
}

############################################################################
# Read the authority file
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

############################################################################
# Iterate through the files and generate the concatenated file
############################################################################

# Then iterate through the files and output the nexus file
open (my $OUTF, ">$outfile.nex") or die "Could not open file $outfile.nex for output.\n";

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
	chomp($phylist[0]);
	(@temparray) = split(/\s+/, $phylist[0]);
	$tempvar = @temparray;
	$ntaxa = $temparray[$tempvar-2];
	$nsites = $temparray[$tempvar-1];
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
				$seq = uc($seq);
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

close($OUTF) or die "Could not close file $outfile.nex.\n";

# Append the sets block to the nexus file
system("sed '1d' $outfile.partitions.txt >> $outfile.nex");

print "\nFile concatenation complete...\n";

exit;
