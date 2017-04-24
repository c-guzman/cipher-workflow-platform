#!/usr/bin/env perl
use warnings;



use POSIX;

if (@ARGV < 2) {
	print STDERR "<group file> <bin file> [option 1=disregard weights]\n";

	print STDERR "Only works with ordinal values\n";
	exit;
}

my $lowerLimit = 1;
my $disregard = 0;
if (@ARGV > 2) {
	$disregard = 1;
	print STDERR "\tDisregarding pre-assigned weights!\n";
}

open IN, $ARGV[0];
my %group = ();
while (<IN> ){
	chomp;
	s/\r//g;
	my @line = split /\t/;
	my $w = 1;
	if (@line > 2 && !$disregard) {
		$w = $line[2];
	}
	if (exists($group{$line[0]})) {
		if ($group{$line[0]}->{'g'} < $line[1]) {
			$group{$line[0]} = {w=>$w,g=>$line[1]};
		}
	} else {
		$group{$line[0]} = {g=>$line[1],w=>$w};
	}
}
close IN;


my $numGenes= 0;
my $numPosGenes = 0;
my %bin = ();
open IN, $ARGV[1];
while (<IN> ){
	chomp;
	s/\r//g;
	my @line = split /\t/;
	next if (!exists($group{$line[0]}));
	if (!exists($bin{$line[1]})) {
		my %a = ();
		$bin{$line[1]} = \%a;
	}
	$numPosGenes+=$group{$line[0]}->{'w'} if ($group{$line[0]}->{'g'} == 1);
	$numGenes+=$group{$line[0]}->{'w'};
	$bin{$line[1]}->{$line[0]} = $group{$line[0]};
}
close IN;

#print STDERR "\tTotal regions=$numGenes\n";
my $rate = $numPosGenes / $numGenes;
my %weight = ();
my %totals = ();
my %ptotals = ();
my $totalP = 0;
my $totalA = 0;
my %tweitght = ();
my $tw = 0;
my $numBin = 0;
print STDERR "\tBin\t# Targets\t# Background\tBackground Weight\n";
#open DIST, ">.dist.tsv";
my @bins = sort {$a <=> $b} keys %bin;
my $newNumGenes = 0;
my $newNumPosGenes = 0;
foreach(@bins) {
	my $b = $_;
	my $a = 0;
	my $p = 0;
	foreach(values %{$bin{$b}}) {
		$p += $_->{'w'} if ($_->{'g'} == 1);
		$a += $_->{'w'};
	}
	next unless ($p >= $lowerLimit && $a-$p >= $lowerLimit);
	$newNumGenes += $a;
	$newNumPosGenes += $p;
}
$numGenes = $newNumGenes;
$numPosGenes = $newNumPosGenes;
#print STDERR "$numGenes\t$numPosGenes\n";
foreach(@bins) {
	my $b = $_;
	my $a = 0;
	my $p = 0;
	foreach(values %{$bin{$b}}) {
		$p += $_->{'w'} if ($_->{'g'} == 1);
		$a += $_->{'w'};
	}
	if ($p >= $lowerLimit && $a-$p >= $lowerLimit) {

		my $allocation = $p/$numPosGenes * ($numGenes-$numPosGenes);

		my $w = $allocation / ($a-$p);

		foreach(keys %{$bin{$b}}) {
			if ($bin{$b}->{$_}->{'g'} == 1) {
				my $W = sprintf("%.5f",1.0 * $bin{$b}->{$_}->{'w'});
				print "$_\t1\t$W\n";
			} else {
				my $W = sprintf("%.5f", $w * $bin{$b}->{$_}->{'w'});
				print "$_\t0\t$W\n";
			}
		}
		my $xxx = $a-$p;
		$w = sprintf("%.3f",$w);
		print STDERR "\t$b\t$p\t$xxx\t$w\n";
#print STDERR "BINDIBN = $b\n";
		#print DIST "B=$b\t$p\t$a\t$w\n";
	}
}
#close DIST;
