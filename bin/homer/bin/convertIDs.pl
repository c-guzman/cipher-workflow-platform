#!/usr/bin/env perl
use warnings;
use lib "/home/carlos/tools/homer/.//bin";
my $homeDir = "/home/carlos/tools/homer/./";

my $promoterSeqOffset = -2000;
my $accDir = $homeDir . "/data/accession/";

if (@ARGV < 3) {
	print STDERR "\n\tUsage: convertIDs.pl <input file> <organism> <ID-type> [header] [keep original] [keep all]\n";
	print STDERR "\tInput file ids: gene, unigene, ensembl, refseq, orf, name, uniprot, affy, agilent\n";
	print STDERR "\tPossible Organisms: mouse, human\n";
	print STDERR "\tPossible ID-type: unigene, gene, refseq, ensembl, orf (use for arabidopsis/yeast), name (e.g. symbol)\n";
	print STDERR "\theader: no, yes (default=no)\n\n";
	print STDERR "\tkeep original: no, yes (default=no)\n\n";
	print STDERR "\tkeep all: no, yes (default=no, excludes lines that aren't in the database)\n\n";
	print STDERR "\texample: convertIDs.pl genelist.txt mouse unigene no no > genelist.unigene.txt\n";
	print STDERR "\n";
	exit;
}

my $inputfile = $ARGV[0];
my $org = $ARGV[1];
my $type = $ARGV[2];
my $header = 0;
my $keep = 0;
if (@ARGV > 3) {
	if ($ARGV[3] eq 'yes' || $ARGV[3] eq 'y') {
		$header = 1;
	}
}
if (@ARGV > 4) {
	if ($ARGV[4] eq 'yes' || $ARGV[4] eq 'y') {
		$keep=1;
	}
}
my $keepAll =0;
if (@ARGV > 5) {
	if ($ARGV[5] eq 'yes' || $ARGV[5] eq 'y') {
		$keepAll=1;
	}
}
my $col = 1;
if ($type eq 'gene') {
	$col=1;
} elsif ($type eq 'unigene') {
	$col=2;
} elsif ($type eq 'refseq') {
	$col =3;
} elsif ($type eq 'ensembl') {
	$col=4;
} elsif ($type eq 'orf' || $type eq 'tair') {
	$col=5;
} elsif ($type eq 'name' || $type eq 'symbol') {
	$col=6;
}

my $convertFile = $accDir . '/' . $org . "2gene.tsv";
open IN, $convertFile or die "could not find ID file for organism = $org ($convertFile does not exist)\n";
close IN;

my $tmpFile = rand() . ".tmp";
my $tmpFile2 = rand() . ".tmp";
my $tmpFile3 = rand() . ".tmp";

if ($keep == 1) {
	`duplicateCol.pl "$inputfile" > "$tmpFile2"`;
	$inputfile = $tmpFile2;
}
`removeAccVersion.pl "$inputfile" > "$tmpFile3"`;
`genericConvertIDs.pl "$tmpFile3" "$convertFile" 0 $col $header 0 $keepAll > "$tmpFile"`;
open IN, $tmpFile;
while (<IN>) {
	print $_;
}
close IN;
`rm "$tmpFile" "$tmpFile3"`;
if ($keep == 1) {
	`rm "$tmpFile2"`;
}

