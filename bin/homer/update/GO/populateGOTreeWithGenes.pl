#!/usr/bin/perl -w



sub printCMD {
	print STDERR "\n\tUsage: populateGOTreeWithGenes.pl <tree> <terms file> <gene2go> [options]\n";
	print STDERR "\n\tCreates *.genes file assigning genes to each ontology group up the tree.\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-org <#> (Organism tax id or \"HG\")\n";
	print STDERR "\t\t-format <gene2go|do_rif|db|HG> (format of gene2go file, default: gene2go from NCBI gene db)\n";
	print STDERR "\n";
	exit;
}

if (@ARGV < 3) {
	printCMD();
}

my $format = "gene2go";
my $org = "";
my $treeFile = $ARGV[0];
my $termsFile = $ARGV[1];
my $gene2goFile = $ARGV[2];

for (my $i=3;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-org') {
		$org = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-format') {
		$format = $ARGV[++$i];
	} else {
		printCMD();
	}
}

my %acceptableEvidenceCodes=();
$acceptableEvidenceCodes{''}=1;

my %def = ();
open IN, $termsFile;
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	$def{$line[0]} = $line[1];
}
close IN;
	

my %ontology = ();
my %terms = ();
my %alt = ();

open IN, $treeFile;
while (<IN>) {
	chomp;
	my @line = split /\t/;
	my @nodes = split (/\,/, $line[1]);
	$ontology{$line[0]} = \@nodes;
	if (!exists($terms{$line[0]})) {
		my %a = ();
		$terms{$line[0]} = \%a;
		$alt{$line[0]} = $line[0];
	}
	foreach(@nodes) {
		if (!exists($terms{$_})) {
			my %a = ();
			$terms{$_} = \%a;
			$alt{$_} = $_;
		}
	}
		
}
close IN;

open IN, $termsFile;
while (<IN>) {
	chomp;
	my @line = split /\t/;
	my $id = $line[0];
	next if (!exists($terms{$id}));
	next if (@line < 3);
	next if ($line[2] eq '');
	my @alt = split /\,/, $line[2];
	foreach(@alt) {
		$alt{$_} = $id;
	}
}
close IN;

my $numDefs = 0;
open IN, $gene2goFile;
while (<IN>) {
	chomp;
	my @line = split /\t/;

	#for ebi
	#my $geneID = $line[0];
	#my $goID = $line[1];

	#for format == gene2go
	next if (@line < 3);
	my $taxID = $line[0];
	my $geneID = $line[1];
	my $goID = $line[2];
	my $evidenceCode = '';
	
	if ($format eq 'gene2go') {
		$evidenceCode = $line[3];
	} elsif ($format eq 'HG') {
		next if (@line < 2);
		$geneID = $line[0];
		$goID = $line[1];
	} elsif ($format eq 'db') {
		next if (@line < 2);
		$geneID = $line[1];
		$goID = $line[0];
	} elsif ($format eq 'do_rif') {
		$geneID = $line[0];
		$goID = $line[4];
	} elsif ($taxID ne '') {
		next if ($taxID ne $org);
	}

	if (!exists($acceptableEvidenceCodes{$evidenceCode})) {
		#next;
	}

	if (!exists($alt{$goID})) {
		if (!exists($def{$goID})) {
			#print STDERR "@line\n";
		}
		next;
	}
$numDefs++;
	$terms{$alt{$goID}}->{$geneID} = 1;
}
close IN;
#print STDERR "$numDefs\n";

my %GENES = ();
printGenes('root');


sub printGenes {
	my ($node) = @_;
	if (exists($GENES{$node})) {
		return $GENES{$node};
	}
	my %genes = ();
	if (exists($ontology{$node})) {
		foreach(@{$ontology{$node}}) {
			my $g = printGenes($_);
			foreach(keys %$g) {
				$genes{$_} =1;
			}
		}
	}
	foreach(keys %{$terms{$node}}) {
		$genes{$_} =1;
	}
	my @genes = keys %genes;
	if (@genes > 0) {
		my $term = "NA";
		if (exists($def{$node})) {
			$term = $def{$node};
		}
		print "$node\t$term\t";
		my $c = 0;
		foreach(keys %genes) {
			print "," if ($c > 0);
			$c++;
			print "$_";
		}
		print "\n";
	}
	$GENES{$node} = \%genes;
	return \%genes;
}
