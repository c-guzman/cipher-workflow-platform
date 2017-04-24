#!/usr/bin/env perl
use warnings;



use POSIX;

if (@ARGV < 2) {
	print STDERR "\n\tUsage: cluster2bed.pl <cluster file> <res> [min percent| 0.05]\n";
	print STDERR "\t  res = resolution used to create the file\n";
	print STDERR "\t  min percent = do not output clusters containing fewer than x percent of the data\n";
	print STDERR "\t  i.e. cluster2bed.pl out.clusters.txt 1000000 0.05 > out.bed\n\n";
	exit;
}
my $pthresh = 0.05;
if (@ARGV > 2) {
	$pthresh = $ARGV[2];
}

print "track name=\"domains, $ARGV[0]\" itemRgb=\"On\"\n";

my $max = 0;
my %counts = ();
open IN, $ARGV[0];
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	$line[0] =~ /(.*?)\-(\d+)$/;
	my $chr = $1;
	my $p = $2;
	if ($p > $max) {
		$max = $p;
	}
	my $c = 0;
	if (@line > 1) {
		$c = $line[1];
		if ($c eq '') {
			$c = 0;
		}
	}
	$counts{$c}++;
	$total++;
}
close IN;

my %bad = ();
my $good = 0;
foreach(keys %counts) {
	if ($counts{$_}/$total < $pthresh) {
		$bad{$_} =1;
	} else {
		$good++;
	}
}

print STDERR "\t$good of $total clusters have more than $pthresh of data\n";

my %colors = ();
open IN, $ARGV[0];
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	$line[0] =~ /(.*?)\-(\d+)$/;
	my $chr = $1;
	my $p = $2;
	next if ($p == $max);
	my $p2 = $p+$ARGV[1];
	my $c = 0;
	if (@line > 1) {
		$c = $line[1];
		if ($c eq '') {
			$c = 0;
		}
	}
	next if (exists($bad{$c}));
	if (!exists($colors{$c})) {
		$colors{$c} = floor(rand()*255) . "," . floor(rand()*255) . "," . floor(rand()*255);
	}
	print "$chr\t$p\t$p2\t$line[0]\t0\t+\t$p\t$p2\t$colors{$c}\n";


}
close IN;
