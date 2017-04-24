#!/usr/bin/env perl
use warnings;
use lib "/home/carlos/tools/homer/.//bin";
my $homeDir = "/home/carlos/tools/homer/./";


my $homologeneFile = $homeDir . "/data/accession/homologene.data";
my $taxFile = $homeDir . "/data/accession/taxids.tsv";

my %taxid = ();
open IN, $taxFile or die "!! Could not find the taxonomy ID file in $homeDir/data/accession/taxid.tsv\nMight need to update the software\n";
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	$taxid{$line[1]} = $line[0];
}
close IN;

sub printCMD() {
	print STDERR "\n\tUsage: convertOrganismID.pl <file> <current organsism> <new organism> <output ID type> [header: yes/no]\n";
	print STDERR "\t\tPrints new file to stdout\n";
	print STDERR "\n\tAvailable Organisms:\n";
	print STDERR "\t\t";
	my $c =0;
	foreach(keys %taxid) {
		print STDERR ", " if ($c > 0);
		$c++;
		print STDERR "$_";
	}
	print STDERR "\n";
	print STDERR "\n\tAvailable Output ID types:\n";
	print STDERR "\t\tgene, refseq, unigene, ensembl\n";
	print STDERR "\n";
	exit;
}

if (@ARGV < 4) {
	printCMD();
}

my $inputfile = $ARGV[0];
my $org1 = $ARGV[1];
my $org2 = $ARGV[2];
my $idtype = $ARGV[3];
my $header = "no";
if (@ARGV > 4) {
	$header = $ARGV[4];
}

if (!exists($taxid{$org1})) {
	print STDERR "!!!! Cannot recognize current organism ($org1)\n";
	printCMD();
}
if (!exists($taxid{$org2})) {
	print STDERR "!!!! Cannot recognize new organism ($org2)\n";
	printCMD();
}
my $taxid1 = $taxid{$org1};
my $taxid2 = $taxid{$org2};

my $rand = rand();
my $tmpfile = $rand . ".tmp";
my $tmpfile2 = $rand . ".2.tmp";
`convertIDs.pl "$inputfile" $org1 gene $header > "$tmpfile"`;

## create conversion table from homologene
my %conv1 = ();
my %conv2 = ();
open IN, $homologeneFile;
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	if ($line[1] eq $taxid1) {
		$conv1{$line[2]} = $line[0];
	}
	if ($line[1] eq $taxid2) {
		$conv2{$line[0]} = $line[2];
	}
}
close IN;

open OUT, ">$tmpfile2";
open IN, $tmpfile;
my $count = 0;
while (<IN>) {
	$count++;
	chomp;
	s/\r//g;
	if ($count == 1 && $header eq 'yes') {
		print OUT "$_\n";
		next;
	}
	my @line = split /\t/;
	if (exists($conv1{$line[0]})) {
		my $homoID= $conv1{$line[0]};
		if (exists($conv2{$homoID})) {
			my $newID = $conv2{$homoID};
			print OUT "$newID";
			for (my $i=1;$i<@line;$i++) {
				print OUT "\t$line[$i]";
			}
			print OUT "\n";
		}
	}
}
close IN;
close OUT;

`convertIDs.pl "$tmpfile2" $org2 $idtype $header > "$tmpfile"`;
open IN, $tmpfile;
while (<IN>) {
	print $_;
}
close IN;
`rm "$tmpfile" "$tmpfile2"`;

