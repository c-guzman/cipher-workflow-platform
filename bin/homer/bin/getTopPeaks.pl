#!/usr/bin/env perl
use warnings;

if (@ARGV < 2) {
	print STDERR "<peak file> <# of peaks> \n";
	print STDERR "outputs to stdout - uses 6th column to find best peaks\n";
	exit;
}

my %peaks = ();
open IN, $ARGV[0];
while (<IN>) {
	chomp;
	s/\r//g;
	my $og = $_;
	my @line = split /\t/;
	next if ($line[0] =~ /^#/);
	next if (@line < 6);
	next unless ($line[5] =~ /^\-*[\d\.e\+\-]+$/);
	$peaks{$line[0]} = {s=>$line[5], og=>$og};
	
}
close IN;

my @peaks = sort {$peaks{$b}->{'s'} <=> $peaks{$a}->{'s'}} keys %peaks;
for (my $i=0;$i<@peaks;$i++) {
	last if ($i >= $ARGV[1]);
	print "$peaks{$peaks[$i]}->{'og'}\n";
}
