#!/usr/bin/env perl
use warnings;
if (@ARGV < 1) {
	print STDERR "<seq file>\n";
	print STDERR "Removes sequences with characters other than ACGTN\n";
	exit;
}
my $total=0;
my $bad = 0;
open IN, $ARGV[0];
while (<IN>) {	
	chomp;
	s/\r//g;
	my @line = split /\t/;
	my $s = $line[1];
	$total++;
	if ($s =~ /[^ACGTN]/) {
		$bad++;
		print STDERR "$line[0]\t$line[1]\n";
	}
	print "$line[0]\t$line[1]\n";
}
close IN;
print STDERR "removed $bad / $total\n";
