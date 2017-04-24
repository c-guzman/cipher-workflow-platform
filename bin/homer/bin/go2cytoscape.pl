#!/usr/bin/env perl
use warnings;


if (@ARGV < 2) {
	print STDERR "\n\tUsage: go2cytosape.pl <output prefix> <Homer GO results file> <organism>\n";
	print STDERR "\n\tThis program will take GO results from a specific organism and output a series\n";
	print STDERR "\tof file that can be loaded into cytoscape (using import functions)\n";
	print STDERR "\n\tInput arguments\n";
	print STDERR "\t\t<output prefix> : The start of the output file names\n";
	print STDERR "\t\t<homer GO results> : Subset of homer GO results you want to visulize.  Try not to use\n";
	print STDERR "\t\t\tthe whole file as this will create a very large cytoscape network\n";
	print STDERR "\t\t<organism> : Used to translate gene ids/accession numbers into gene names\n";
	print STDERR "\n\tOutput Files: (for cytoscape)\n";
	print STDERR "\t\t<prefix>.network.sif.txt - load in cytoscape using import->network from text\n";
	print STDERR "\t\t<prefix>.node.sizes.txt - load in cytoscape as a node attribute\n";
	print STDERR "\n";
	exit;
}

my $prefix = $ARGV[0];
my $input = $ARGV[1];
my $organism = $ARGV[2];
my $rand = rand();
my $tmpFile = $rand . ".tmp";
my $tmpFile2 = $rand . ".2.tmp";
open OUT, ">$tmpFile";
my %genes = ();
open IN, $input;
while (<IN>) {
	chomp;
	s/\r//g;
	s/\"//g;
	my @line = split /\t/;
	next if (@line < 10);
	next if ($line[0] eq 'GO');
	my @genes = split /\,/, $line[9];
	foreach(@genes) {
		if (!exists($genes{$_})) {
			print OUT "$_\t$_\n";
		}
		$genes{$_} = "";
	}
	my $name = $line[1];
	$name =~ s/\s/_/g;
	my $go = {name=>$name,id=>$line[0],logp=>$line[3],g=>\@genes};
	push(@GO, $go);
}
close IN;

`convertIDs.pl "$tmpFile" $organism name no no yes > "$tmpFile2"`;
open IN, $tmpFile2;
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	next if (@line < 2);
	next if (!exists($genes{$line[1]}));
	$genes{$line[1]} = $line[0];
}
close IN;
`rm "$tmpFile" "$tmpFile2"`;

open NETWORK, ">$prefix.network.sif.txt";
open NODES, ">$prefix.node.sizes.txt";
print NODES "Enrichment\n";

foreach(@GO) {
	my $go = $_;
	my $name = $go->{'name'};
	foreach(@{$go->{'g'}}) {
		next if (!exists($genes{$_}));
		my $n = $genes{$_};
		print NETWORK "$name\tpp\t$n\n";
	}
	my $logp = $go->{'logp'};
	print NODES "$name = $logp\n";
}
foreach(values %genes) {
	print NODES "$_ = 0\n";
}
close NODES;
close NETWORK;




