#!/usr/bin/env perl
use warnings;



#
use POSIX;

if (@ARGV < 1) {
	print STDERR "<Sequence file> [% match i.e. 0.5]\n";
	exit;
}

my $score = 0.5;
if (@ARGV > 1) {
	$score = $ARGV[1];
}
my $name = rand();

`tab2fasta.pl "$ARGV[0]" > "$name.fa"`;
`blat "$name.fa" "$name.fa" "$name.psl"`;

my %data = ();

open IN, "$name.psl";
my $count = 0;
while (<IN>) {
	chomp;
	$count++;
	next if ($count < 6);
	my @line = split /\t/;
	
	$minScore = $score*$line[10];

	next if ($line[0] < $minScore);

	my $id1 = $line[9];
	my $id2 = $line[13];
	next if ($id1 eq $id2);
	if (!exists($data{$id1})) {
		my %a = ();
		$data{$id1} = \%a;
	}
	if (!exists($data{$id2})) {
		my %a = ();
		$data{$id2} = \%a;
	}
	$data{$id1}->{$id2}=1;
	$data{$id2}->{$id1}=1;
}
foreach(keys %data) {
	print "$_\t";
	my $c=0;
	foreach(keys %{$data{$_}}) {
		print "," if ($c>0);
		$c++;
		print "$_";
	}
	print "\n";
}

`rm "$name.fa" "$name.psl"`;
