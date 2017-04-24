#!/usr/bin/env perl
use warnings;


if (@ARGV < 1) {
	print STDERR "<fastq file>\n";
	print STDERR "Will output a FASTA file\n";
	exit;
}


open IN, $ARGV[0];
my $good = 0;
my $numLines = 0;
while (<IN>) {
	$numLines++;
	chomp;
	if (/^\@/) {
		next if ($good == 0 && $numLines==1);
		$good = 1;
		s/^\@//;
		s/ /_/g;
		print ">$_\n";
	} elsif (/^\+/) {
		$good = 0;
		$numLines = 0;
	} else {
		if ($good) {
			print "$_\n";
		}
	}
}
close IN;
