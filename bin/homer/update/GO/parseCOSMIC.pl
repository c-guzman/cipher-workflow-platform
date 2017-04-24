#!/usr/bin/perl -w
#


my $minNumGenesPerGroup = 3;
if (@ARGV < 1) {
	print STDERR "\n\tUsage: parseCOSMIC.pl CosmicCompleteExport.tsv human2gene.tsv\n";
	print STDERR "\n\t\tMinimum number of genes per group: $minNumGenesPerGroup\n";
	print STDERR "\n";
	exit;
}
my %genes = ();

my $cosmicDatabase = $ARGV[0];
my $geneFile = $ARGV[1];

#`wget -O CosmicCompleteExport_v67_241013.tsv.gz ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/data_export/CosmicCompleteExport_v67_241013.tsv.gz`;
#`gunzip CosmicCompleteExport_v67_241013.tsv.gz`;


my %conv = ();
open IN, $geneFile or die "!! Could not open $geneFile !!\n";
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	$conv{$line[0]} = $line[1];
}
close IN;

open IN, $cosmicDatabase;
my $c = 0;
while (<IN>) {
	$c++;
	next if ($c < 2);
	chomp;
	s/\r//g;
	my @line = split /\t/;
	next if (@line < 10);
	my @ids = ();
	push(@ids, $line[7], $line[8], $line[9], $line[10]);
	my $metaid = $line[7] . "-" . $line[8] . "-" . $line[9] . "-" . $line[10];
	$metaid =~ s/\-NS//g;
	push(@ids, $metaid);
	#push(@ids, "$line[4]_" . $metaid);
	
	my $g = $line[0];

	if (!$conv{$g}) {
		next;
	}
	$g = $conv{$g};

	foreach(@ids) {
		my $id = $_;
		next if ($id eq 'NS');
		next if ($id eq '');
		if (!exists($genes{$id})) {
			my %a = ();
			$genes{$id} = \%a;
		}
		$genes{$id}->{$g}=1;
	}

}
close IN;

foreach(keys %genes) {
	my $id = $_;
	my $c = 0;
	my @genes = keys %{$genes{$id}};
	my $N = scalar(@genes);
	next if ($N < 3);
	print "$id\t$id\t";
	foreach(@genes) {
		print "," if ($c > 0);
		$c++;
		print "$_";
	}
	print "\n";
	#print STDERR "$id\t$c\n";
}

#`rm CosmicCompleteExport_v67_241013.tsv`;
