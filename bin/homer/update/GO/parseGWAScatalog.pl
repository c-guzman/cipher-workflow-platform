#!/usr/bin/perl -w
#

my %genes = ();
`wget -O GWASCatalogDump.txt ftp://ftp.ncbi.nlm.nih.gov/dbgap/GWASCatalog/GWASCatalogDump.txt`;
open IN, "GWASCatalogDump.txt";
my $c = 0;
while (<IN>) {
	$c++;
	next if ($c < 2);
	chomp;
	s/\r//g;
	my @line = split /\t/;
	next if (@line < 23);
	my $g = $line[22];
	my $id = $line[11] . "(PMID:$line[5])";
	if (!exists($genes{$id})) {
		my %a = ();
		$genes{$id} = \%a;
	}
	$genes{$id}->{$g}=1;

}
close IN;

foreach(keys %genes) {
	my $id = $_;
	print "$id\t$id\t";
	my $c = 0;
	foreach(keys %{$genes{$id}}) {
		print "," if ($c > 0);
		$c++;
		print "$_";
	}
	print "\n";
}

`rm GWASCatalogDump.txt`;
