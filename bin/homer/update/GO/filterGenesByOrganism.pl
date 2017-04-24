#!/usr/bin/perl -w
#

if (@ARGV < 3) {
	print STDERR "\n\tUsage: filterGenesByOrganism.pl <organism name> <organ.description file> <ontolgoy.genes file> [ontology.genes file2] ...\n";
	print STDERR "\n\tWill create new files with the organism name as the prefix for each file.\n";
	print STDERR "\n";
	exit;
}

my $orgname= $ARGV[0];
my $orgfile = $ARGV[1];

my %ids = ();
open IN, $orgfile;
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	$ids{$line[0]} = 1;
}
close IN;

for (my $i=2;$i<@ARGV;$i++) {
	my $newfile = $orgname . "." . $ARGV[$i];
	open OUT, ">$newfile";
	open IN, $ARGV[$i] or print STDERR "!!! Could not open $ARGV[$i] (update/GO/filterGenesByOrganism.pl)\n";
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		next if (@line < 3);
		my @genes = split /\,/,$line[2];
		my $c = 0;
		my $newgenes = "";
		foreach(@genes) {
			if (exists($ids{$_})) {
				if ($c>0) {
					$newgenes .= ",";
				}
				$c++;
				$newgenes .= $_;
			}	
		}
		if ($c>0) {
			print OUT "$line[0]\t$line[1]\t$newgenes\n";
		}
	}
	close IN;
	close OUT;
}
				
