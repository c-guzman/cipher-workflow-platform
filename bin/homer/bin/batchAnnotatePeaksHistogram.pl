#!/usr/bin/env perl
use warnings;



if (@ARGV < 3) {
	print STDERR "\n\tFor making histograms across multiple samples using the same parameters.\n";
	print STDERR "\tusage:\n";
	print STDERR "\t\tusage:batchAnnotatePeaks.pl <genome> [annotatePeaks.pl options] -f <peak file 1> ...\n";
	print STDERR "\t\t\tEverything after -f will be treated as a peak file\n";
	print STDERR "\n";
	exit;
}


my $cmdStr = "";
my @files = ();
my $histFlag = 0;
my $fileMode = 0;
for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-f') {
		$fileMode = 1;
		next;
	}
	if ($fileMode) {
		push(@files, $ARGV[$i]);
	} else {
		$cmdStr .= " \"" . $ARGV[$i] . '"';
		if ($ARGV[$i] eq '-hist') {
			$histFlag = 1;
		}
	}
}

unless ($histFlag) {
	print STDERR "!!!! warning - \"-hist\" option not used!!!\n";
	exit;
}

my $rand = rand();
my $tmpFile = $rand . ".tmp";
my $tmpFile2 = $rand . ".2.tmp";
my $tmpFile3 = $rand . ".3.tmp";

for (my $i=0;$i<@files;$i++) {
	
	`annotatePeaks.pl "$files[$i]" $cmdStr > $tmpFile`;
	if ($i==0) {
		`cp $tmpFile $tmpFile2`;
	} else {
		`addDataHeader.pl $tmpFile2 $tmpFile > $tmpFile3`;
		`mv $tmpFile3 $tmpFile2`;
	}
}
open IN, $tmpFile2;
while (<IN>) {
	print $_;
}
close IN;

`rm $tmpFile $tmpFile2`;
