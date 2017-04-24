#!/usr/bin/perl -w
if (@ARGV < 2) {
	print STDERR "<gene interactions> <gene_info>\n";
	exit;
}

my %pi = ();

open IN, $ARGV[0];
while (<IN>) {
	chomp;
	my @line = split /\t/;
	my $id1 = $line[1];
	my $id2 = $line[6];
	next if ($line[7] ne 'GeneID');
	next if ($id1 eq '-' || $id1 eq '');
	next if ($id2 eq '-' || $id2 eq '');
	if (!exists($pi{$id1})){ 
		my %a = ();
		$pi{$id1}=\%a;
	}
	if (!exists($pi{$id2})) {
		my %a = ();
		$pi{$id2} =\%a;
	}
	$pi{$id1}->{$id2}=1;
	$pi{$id2}->{$id1}=1;
}
close IN;

my %terms = ();
open IN, $ARGV[1];
while (<IN>) {
	chomp;
	my @line= split /\t/;
	next if (@line < 9);
	my $gid = $line[1];
	next if (!exists($pi{$gid}));
	my $name = $line[2];
	my $desc = $line[8];
	my $info = $line[2] . " (" . $line[8] . ")";

	print "$gid\t$info\t";
	my $c = 1;
	foreach(keys %{$pi{$gid}}) {
		print "," if ($c>1);
		$c++;
		print "$_";
	}
	print "\n";
}
close IN;
