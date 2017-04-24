#!/usr/bin/env perl
use warnings;


if (@ARGV < 1) {
	print STDERR "usage: revoppMotif.pl <motif file>\n";
	print STDERR "Prints new motif file (of reverse opposite) to stdout\n";
	exit;
}
my @profile = ();
my $p = \@profile;
open IN, $ARGV[0];
my $start = 1;
while (<IN>) {
	chomp;
	s/\r//g;
	if (/^>/) {
		if ($start != 1) {
			printProfile($p);
			print "$_\n";
		} else {
			print "$_\n";
			$start = 0;
		}
		my @a = ();
		$p = \@a;
		next;
	}
	my @line = split /\t/;
	push(@$p, \@line);
}
printProfile($p);

sub printProfile {
	my ($m) = @_;
	for (my $i=@$m-1;$i>=0;$i--) {
		for (my $j=3;$j>=0;$j--) {
			print "\t" if ($j!=3);
			print $m->[$i]->[$j];
			
		}
		print "\n";
	}


}
	
