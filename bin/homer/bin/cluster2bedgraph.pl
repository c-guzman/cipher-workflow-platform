#!/usr/bin/env perl
use warnings;



use POSIX;

if (@ARGV < 2) {
	print STDERR "\n\tUsage: cluster2bedgraph.pl <cluster distance file> [res] [name]\n";
	print STDERR "\t  res = resolution used to create the file\n";
	print STDERR "\n";
	exit;
}

my $trackName = "Correlation to Cluster";
my $res = $ARGV[1];
if (@ARGV > 2) {
	$trackName = $ARGV[2];
}
my $max = 0;
my %counts = ();
my %regions = ();
my $count = 0;
my @header = ();
open IN, $ARGV[0];
while (<IN>) {
	$count++;
	chomp;
	s/\r//g;
	s/\"//g;
	my @line = split /\t/;
	if ($count == 1) {
		shift @line;
		@header = @line;
		next;
	}
		
	$line[0] =~ /(.*?)\-(\d+)$/;
	my $chr = $1;
	my $p = $2;
	if (!exists($max{$chr})) {
		$max{$chr}=$p;
	} elsif ($p > $max{$chr}) {
		$max{$chr} = $p;
	}
	shift @line;
	$data{$line[0]} = {c=>$chr,p=>$p,d=>\@line};
}
close IN;

my @regions = sort {$data{$a}->{'c'} cmp $data{$b}->{'c'} || $data{$a}->{'p'} <=> $data{$b}->{'p'}} keys %data;

for (my $i=0;$i<@header;$i++) {
	my $name = $header[$i];
	print "track name=\"$trackName $name\"  yLineMark=\"0.0\" alwaysZero=on maxHeightPixels=100:75:11 visibility=full viewLimits=-1:1 autoScale=on type=bedGraph\n";
	foreach(@regions) {
		my $c = $data{$_}->{'c'};
		my $p = $data{$_}->{'p'};
		my $p2  = $p + $res;
		next if ($p >= $max{$c});
		my $v = $data{$_}->{'d'}->[$i];
		print "$c\t$p\t$p2\t$v\n";
	}
}
