#!/usr/bin/env perl
use warnings;


if (@ARGV < 1) {
	print STDERR "\n\tusage: getChrLengths.pl <fasta file> [fasta file2] ...\n";
	print STDERR "\tReturns the lengths of each sequence in the fasta files\n";
	print STDERR "\n";
	exit;
}

print "Chromosome\tLength(bp)\n";
for (my $i=0;$i<@ARGV;$i++) {
	open IN, $ARGV[$i];
	my $currentLength = 0;
	my $currentName = "";
	while (<IN>) {
		chomp;
		s/\r//g;
		if (/^>/) {
			if ($currentName ne '') {
				print "$currentName\t$currentLength\n";
			}
			s/^>//;
			s/\s.*//;
			$currentName = $_;
			$currentLength = 0;
			next;
		}
		$currentLength += length($_);
	}
	close IN;
	if ($currentName ne '') {
		print "$currentName\t$currentLength\n";
	}
}
