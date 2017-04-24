#!/usr/bin/env perl
use warnings;

my $col = 1;
@freq = (0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10,0.12,0.14,0.16);

if (@ARGV < 1) {
	print STDERR "<freq file> [column, zero index, default=1] [freq bin 1] [freq bin2]...\n";
	exit;
}
if (@ARGV > 1) {
	$col= $ARGV[1];
}
if (@ARGV > 2) {
	@freq = ();
	for (my $i=2;$i<@ARGV;$i++) {
		push(@freq, $ARGV[$i]);
	}
}

if ($col == 2) {
	@freq = (0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8);
}

print STDERR "\tFrequency Bins: @freq\n";

@freq = sort {$a <=> $b} @freq;
my $m = @freq;
open IN, $ARGV[0];
my %counts = ();
while (<IN>) {
	chomp;
	my @line = split /\t/;
	next if ($line[1] eq 'NA');
	print "$line[0]\t";
	my $c = 0;
	for (my $i=0;$i<@freq;$i++) {
		if ($line[$col] < $freq[$i]) {
			print "$i\n";
			$counts{$i}++;
			$c=1;
			last;
		}
	}
	if ($c == 0) {
		$counts{$m}++;
		print "$m\n";
	}
}
close IN;
push(@freq, $m);
my @z = sort {$a <=> $b} keys %counts;
print STDERR "\tFreq\tBin\tCount\n";
foreach(@z) {
	print STDERR "\t$freq[$_]\t$_\t$counts{$_}\n";
}
