#!/usr/bin/env perl
use warnings;

sub printCMD {
	print STDERR "\n\tUsage: getFocalPeaks.pl <peaks.centered.txt> [options]\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-min <#> (minimum focus ratio threshold, default: 0.75)\n";
	print STDERR "\t\t-max <#> (maximum focus ratio threshold, default: none)\n";
	print STDERR "\n";
	exit;
}

my $minThresh = 0.75;
my $maxThresh = 1e20;
if (@ARGV < 1) {
	printCMD();
}
for (my $i=1;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-min') {
		$minThresh = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-max') {
		$maxThresh = $ARGV[++$i];
	} else {
		printCMD();
	}
}

my $total=0;
my $good =0;
my $count = 0;
open IN, $ARGV[0];
while (<IN>) {
	$count++;
	next if (/^#/);
	if ($count < 2) {
		print $_;
		next;
	}
	chomp;
	s/\r//g;
	my $og = $_;
	my @line = split /\t/;
	next if (@line < 7);
	$total++;
	next if ($line[6] <= $minThresh);
	next if ($line[6] >= $maxThresh);
	print "$og\n";
	$good++;
}
my $percent = sprintf("%.2f",100*$good/$total);
print STDERR "\n\tGood: $good/$total ($percent%)\n\n";
