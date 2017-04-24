#!/usr/bin/env perl
use warnings;
use lib "/home/carlos/tools/homer/.//bin";
my $homeDir = "/home/carlos/tools/homer/./";

use Statistics;

if (@ARGV < 1) {
	print STDERR "\n\tInternal program called by other scripts\n";
	print STDERR "\tprepForR.pl <data file> [options]\n";
	print STDERR "\n";
	exit;
}
my $mode = 'corrMatrix';
if (@ARGV > 1) {
	$mode = $ARGV[1];
}

if ($mode eq 'corrMatrix') {	
	my @good = ();
	open IN, $ARGV[0];
	my $count=0;
	while (<IN>) {
		$count++;
		if ($count == 1) {
			push(@good,1);
			next;
		}
		chomp;
		s/\r//g;
		my @line = split /\t/;
		shift @line;
		shift @line;
		my ($avg,$var) = Statistics::avevar(\@line);
		if ($var < 1e-20) {
			push(@good, 0);
		} else {
			push(@good, 1);
		}
	}
	close IN;

	open IN, $ARGV[0];
	$count=0;
	while (<IN>) {
		$count++;
		chomp;
		s/\r//g;
		my @line = split /\t/;
		shift @line;
		next if ($good[$count-1]==0);
		print "$line[0]" if ($count > 1);
		for (my $i=1;$i<@good;$i++) {
			if ($good[$i] == 1) {
				print "\t$line[$i]";
			}
		}	
		print "\n";
	}
}
