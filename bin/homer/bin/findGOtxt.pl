#!/usr/bin/env perl
use warnings;
use lib "/home/carlos/tools/homer/.//bin";
my $homeDir = "/home/carlos/tools/homer/./";
my $promoterSeqOffset = -2000;

if (@ARGV < 2) {
	print STDERR "<go.genes> <base> <list> <1=restrict flag, 0=don't restrict>\n";
	#print STDERR " or <go.genes> <group file>\n";
	print STDERR "\tgo.genes: ontology 'genes' file\n";
	print STDERR "\tList of IDs to use as a base list\n";
	print STDERR "\tList of IDs in target group\n";
	print STDERR "\n\tRestrict business: restricts base to be all genes represented in *.genes file\n";
	exit;
}
use POSIX;
use Statistics;

my $restrictFlag = 1;
my $baseFlag = 1;
my %base = ();
my %list = ();
if (open IN, $ARGV[1]) {
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		$base{$line[0]} = 0;
		if (@ARGV < 3) {
			if ($line[1] eq '1') {
				$list{$line[0]} = 0;
			}
		}
	}
	close IN;
} else {
	$baseFlag = 0;
}

if (@ARGV > 2) {
	open IN, $ARGV[2];
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		if ($baseFlag == 1) {
			next if (!exists($base{$line[0]}));
		}
		$list{$line[0]} = 0;
		$base{$line[0]} = 0;
	}
	close IN;
}
if (@ARGV > 3) {
	$restrictFlag= $ARGV[3];
}


my %nbase = ();
my %nlist = ();
my %ontology = ();
open IN, $ARGV[0] or die "Could not open $ARGV[0]\n";
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	next if (@line < 3);
	my @genes = split /\,/, $line[2];
	my %uniq = ();
	foreach(@genes) {
		$uniq{$_}=1;
	}
	@genes = keys %uniq;
	my %genes = ();
	my $n2 = 0;
	my $n = 0;
	my %common = ();
	foreach(@genes) {
		if ($baseFlag == 1) {
			next if (!exists($base{$_}));
		} else {
			$base{$_} = 0;
		}
		$nbase{$_} = 1;
		$n2++;
		if (exists($list{$_})) {
			$n++;
			$nlist{$_} = 1;
			$common{$_} = $_;
		}
	}
	next if ($n2 == 0);
	$ontology{$line[0]} = {id=>$line[1],n2=>$n2,n=>$n,genes=>\%common};
}
close IN;

my $N = scalar(keys %nbase);
my $n1 = scalar(keys %nlist);

foreach(keys %ontology) {
	my $go = $_;
	my $n = $ontology{$go}->{'n'};
	my $n2 = $ontology{$go}->{'n2'};
	$ontology{$go}->{'logp'} = Statistics::loghypergeo($N, $n1, $n2, $n);
	$ontology{$go}->{'p'} = exp($ontology{$go}->{'logp'});
	my $d = $n1;
	$d = 1 if ($n1 == 0);
	$ontology{$go}->{'fraction'} = $n/$d;
}
my @go = sort {$ontology{$a}->{'logp'} <=> $ontology{$b}->{'logp'}} keys %ontology;
print "GO\tTerm\tP-value\tLogP\tNumber in Term\tNumber in common\tFraction of List\tTotal List\tTotal Genes\tCommon Gene IDs\n";
foreach(@go){ 
	print "$_\t$ontology{$_}->{'id'}\t$ontology{$_}->{'p'}\t$ontology{$_}->{'logp'}\t$ontology{$_}->{'n2'}\t$ontology{$_}->{'n'}\t$ontology{$_}->{'fraction'}\t$n1\t$N\t";

	my $z = 0;
	foreach(keys %{$ontology{$_}->{'genes'}}) {
		print "," if ($z>0);
		$z++;
		print "$_";
	}
	print "\n";
	
}
