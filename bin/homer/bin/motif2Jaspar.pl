#!/usr/bin/env perl
use warnings;

if (@ARGV < 1) {
	print STDERR "\n\tUsage: motif2Jaspar.pl <motif file>\n";
	print STDERR "\tOutputs a Jaspar formatted motif file to stdout\n";
	print STDERR "\n";
	exit;
}
my @m = ();
open IN, $ARGV[0];
while (<IN>) {
	chomp;
	s/\r//g;
	next if (/^>/);
	my @line = split /\t/;
	push(@m, \@line);
}
close IN;
my @alpha = ('A','C','G','T');
for (my $j=0;$j<@alpha;$j++) {
	print "$alpha[$j] [";
	for (my $i=0;$i<@m;$i++){ 
		print " $m[$i]->[$j]";
	}
	print "]\n";
}

