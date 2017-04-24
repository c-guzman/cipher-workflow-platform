#!/usr/bin/env perl
use warnings;

if (@ARGV < 1) {
	print STDERR "\n\tgff2pos.pl <gff formated file>\n";
	print STDERR "\tSends peak/pos file to stdout\n";
	print STDERR "\n";
	exit;
}

my %usedIDs = ();
open IN, $ARGV[0];
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	next if (/^\s*#/);
	my $id = $line[1];
	if (exists($usedIDs{$id})) {
		$usedIDs{$id}++;
		$id .= "-" . $usedIDs{$id};
	} else {
		$usedIDs{$id} = 1;
	}
	print "$id\t$line[0]\t$line[3]\t$line[4]\t$line[6]\t$line[5]\n";
}
close IN;
