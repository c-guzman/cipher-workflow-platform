#!/usr/bin/perl -w
#

if (@ARGV < 1) {
	print STDERR "Usage: human2gene.tsv\n";
	exit;
}
my %genes = ();
`wget -O smpdb_proteins.csv.zip http://www.smpdb.ca/downloads/smpdb_proteins.csv.zip`;
`unzip -u smpdb_proteins.csv.zip`;

my %conv = ();
open IN, $ARGV[0] or die "!! Could not open $ARGV[0] !!\n";
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	$conv{$line[0]} = $line[1];
}
close IN;

open IN, "proteins.csv" || die "!!! Could not open proteins.csv\n";
my $c = 0;
while (<IN>) {
	$c++;
	next if ($c < 2);
	chomp;
	s/\r//g;
	s/\"(.*?)\,(.*?)\"/$1\-$2/g;
	my @line = split /\,/;
	next if (@line < 5);
	my $id = $line[0] . "\t" . $line[1];
	#my $g = $line[4];
	next if (@line < 9);
	my $g = $line[8];

	if (!$conv{$g}) {
		next;
	}
	$g = $conv{$g};

	if (!exists($genes{$id})) {
		my %a = ();
		$genes{$id} = \%a;
	}
	$genes{$id}->{$g}=1;

}
close IN;
#print STDERR "$c\n";

foreach(keys %genes) {
	my $id = $_;
	print "$id\t";
	my $c = 0;
	foreach(keys %{$genes{$id}}) {
		print "," if ($c > 0);
		$c++;
		print "$_";
	}
	print "\n";
}

`rm proteins.csv smpdb_proteins.csv.zip`;
