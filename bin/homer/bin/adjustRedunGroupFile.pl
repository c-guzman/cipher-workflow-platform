#!/usr/bin/env perl
use warnings;



if (@ARGV < 1) {
	print STDERR "<groupFile> [M|H| redundant file]\n";
	exit;
}


my %extra = ();
my $extraFlag = 0;
my %group = ();
my %order = ();
my $Z = 0;
open IN, $ARGV[0];
while (<IN>) {
	chomp;
	$Z++;
	my @line = split /\t/;
	if (@line < 2) {
		$group{$line[0]} = 0;
		$order{$line[0]} = 0;
	} else {
		my $gene = shift @line;
		my $group = shift @line;
		if (exists($group{$gene})) {
			if ($group{$gene} < $group) {
				$group{$gene} = $group;
				$order{$gene} = $Z;
				if (@line > 0) {
					$extraFlag = 1;
					$extra{$gene} = \@line;
				}
			}
		} else {
			$group{$gene} = $group;
			$order{$gene} = $Z;
			if (@line > 0) {
				$extraFlag = 1;
				$extra{$gene} = \@line;
			}
		}
	}

}
close IN;

my %taken = ();

my $f = "/home/cbenner/Transcription/newcpp/q/redunMouse500-100.tsv";
if (@ARGV > 1) {
	if ($ARGV[1] eq 'M') {

	} elsif ($ARGV[1] eq 'H') {
		$f = "/home/cbenner/Transcription/newcpp/q/redunHuman500-100.tsv";
	} else {
		$f = $ARGV[1];
	}
}
open IN, $f;
my %conv = ();
while (<IN>) {
	chomp;
	my @line = split /\t/;
	my @genes = split /\,/, $line[1];
	$conv{$line[0]} = \@genes;
}
close IN;

my @g = sort {$group{$b} <=> $group{$a} || $order{$a}<=>$order{$b}} keys %group;
my $numGood = 0;
my $totalNum =0;
foreach(@g) {
	$totalNum++;
	next if (exists($taken{$_}));
	my $bad = 0;
	foreach(@{$conv{$_}}) {
		if (exists($taken{$_})) {
			$bad = 1;
			last;
		}
	}
	if ($bad) {
		next;
	}
	$taken{$_} = 1;
	$numGood++;
	print "$_\t$group{$_}";
	if ($extraFlag) {
		foreach(@{$extra{$_}}) {
			print "\t$_";
		}
	}
	print "\n";
}
print STDERR "\t\tKept $numGood of $totalNum\n";


