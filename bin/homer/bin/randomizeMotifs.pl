#!/usr/bin/env perl
use warnings;

if (@ARGV < 1) {
	print STDERR "\n\tUsage: randomizeMotifs.pl <motif file> [0|1] [0|1]\n";
	print STDERR "\t1st [0|1]: enter 1 randomization of positions, 0 to skip\n";
	print STDERR "\t2nd [0|1]: enter 1 invert bases A->T, C->G, G->C, T->A, 0 to skip\n";
	print STDERR "\tDefault: is both set to 1\n";
	print STDERR "\n";
	exit;
}
my $randOrderFlag = 1;
my $invertFlag = 1;

if (@ARGV > 1) {
	$randOrderFlag = $ARGV[1];
}
if (@ARGV > 2) {
	$invertFlag = $ARGV[2];
}

my @motifs = ();
my @m = ();
my $mm = \@m;
open IN, $ARGV[0];
while (<IN>) {
	chomp;
	s/\r//g;

	my @line = split /\t/;
	if ($line[0] =~ /^>/) {
		if (scalar(@$mm) > 1) {
			push(@motifs, $mm);
			my @m = ();
			$mm = \@m;
		}
	}
	push(@$mm, \@line);
}
close IN;
push(@motifs, $mm);

foreach(@motifs) {
	my $m = $_;
	my $len = scalar(@$m)-1;
	my %randOrder = ();
	for (my $i=1;$i<=$len;$i++) {
		$randOrder{$i} = rand();
	}
	print "$m->[0][0]\t$m->[0][1]";
	if ($randOrderFlag) {
		print "-ColumnRandomized";
	}
	if ($invertFlag) {
		print "-Inverted";
	}
	for (my $i=2;$i<@{$m->[0]};$i++) {
		print "\t$m->[0][$i]";
	}
	print "\n";
	my @index = sort {$a <=> $b} keys %randOrder;
	if ($randOrderFlag) {
		@index = sort {$randOrder{$a} <=> $randOrder{$b}} keys %randOrder;
	}
	foreach(@index) {
		if ($invertFlag) {
			print "$m->[$_][3]";
			for (my $i=2;$i>=0;$i--) {
				print "\t$m->[$_][$i]";
			}
			print "\n";
		} else {
			print "$m->[$_][0]";
			for (my $i=1;$i<4;$i++) {
				print "\t$m->[$_][$i]";
			}
			print "\n";
		}
	}
				
}
	
