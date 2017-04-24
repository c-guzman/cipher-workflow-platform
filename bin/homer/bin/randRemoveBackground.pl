#!/usr/bin/env perl
use warnings;


use POSIX;
if (@ARGV < 1) {
	print STDERR "<groupfile> <total # of sequences> <bins | null>\n";
	exit;
}
my $max = $ARGV[1];
my @background = ();
my $n = 0;
my %bins = ();
my %backbins = ();
my %posbins = ();


my %bb = ();
if (@ARGV > 2 && $ARGV[2] ne 'null') {
	open IN, $ARGV[2];
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		$bins{$line[0]} = $line[1];
		if (!exists($bb{$line[1]})) {
			my @a = ();
			$bb{$line[1]} = \@a;
		}
	}
	close IN;
} else {
	open IN, $ARGV[0];
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		$bins{$line[0]} = 1;
		if (!exists($bb{1})) {
			my @a = ();
			$bb{1} = \@a;
		}
	}
	close IN;
}

open IN, $ARGV[0];
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	next if (!exists($bins{$line[0]}));
	if ($line[1] == 1) {
		$n++;
		print "$line[0]\t1\n";
		$posbins{$bins{$line[0]}}++;
	} else {
		push(@background, $line[0]);
		$backbins{$bins{$line[0]}}++;
		push(@{$bb{$bins{$line[0]}}}, $line[0]);
	}
}
close IN;

my $num = $max-$n;

if ($num < $n) {
	print STDERR "Less background than target sequences!!!!!\n";
#	exit;
}
if ($num < 1) {
	print STDERR "Something is wrong - not enough background sequences!!!!\n";
	exit;
}
my %numBack = ();
my $totalPos = 0;
my @bins = keys %posbins;
my @realbins = ();
foreach(@bins) {
	next if (!exists($backbins{$_}));
	next if ($backbins{$_} < 1);
	$totalPos += $posbins{$_};
	push(@realbins, $_);
}
@realbins = sort {$b <=> $a} @realbins;
my $leftover = 0;
foreach(@realbins) {
	my $bin = $_;
	my $num2add = floor($posbins{$bin}/$totalPos*$num+0.49999999);
	$num2add+=$leftover;
	$leftover = 0;

	my @ids = @{$bb{$bin}};
	my $nids = scalar(@ids);
	my $z =0;

	if ($num2add > $nids/2) {

		my %ids = ();
		foreach(@ids) {
			$ids{$_} = rand();
		}
		@ids = sort {$ids{$a} <=> $ids{$b}} @ids;
		for (my $i=0;$i<$num2add;$i++) {
			last if ($i>=@ids);
			print "$ids[$i]\t0\n";
			$z++;
		}
	} else {
		my %mask = ();
		while ($z < $num2add) {
			my $index = floor(rand()*$nids);
			next if ($index >= $nids);
			next if (exists($mask{$index}));
			$mask{$index}=1;
			print "$ids[$index]\t0\n";
			$z++;
		}

	}
	if ($z < $num2add) {
		$leftover += $num2add-$z;
	}
	#print STDERR "\t\t$bin\t$posbins{$bin}\t$num2add\t$z\n";
}
