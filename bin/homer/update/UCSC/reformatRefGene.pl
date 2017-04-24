#!/usr/bin/perl -w
if (@ARGV < 2) {
	print STDERR "<refGene.txt> <genome>\n";
	exit;
}
my $genome = $ARGV[1];

open IN, $ARGV[0];
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	if ($genome eq 'hg17') {
		for (my $i=0;$i<@line;$i++) {
			print "\t$line[$i]";
		}
		print "\n";
	}
}
close IN;
