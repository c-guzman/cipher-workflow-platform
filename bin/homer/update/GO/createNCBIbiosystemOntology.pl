#!/usr/bin/perl -w
#
print STDERR "\tNo command line arguments needed...\n";

my $keepFiles = 0;
if ($keepFiles == 0) {
	`wget -O bsid2info.gz ftp://ftp.ncbi.nih.gov/pub/biosystems/CURRENT/bsid2info.gz`;
	`wget -O biosystems_gene.gz ftp://ftp.ncbi.nih.gov/pub/biosystems/CURRENT/biosystems_gene.gz`;
	
	`gunzip -f bsid2info.gz`;
	`gunzip -f biosystems_gene.gz`;
}

my %genes = ();
open IN, "biosystems_gene";
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	if (!exists($genes{$line[0]})) {
		my @a = ();
		$genes{$line[0]} = \@a;
	}
	push(@{$genes{$line[0]}},$line[1]);
}
close IN;
my %data = ();
open IN, "bsid2info";
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	my $id =$line[0];
	my $source = $line[1];
	my $sid = $line[2];
	my $name = $line[3];
	next if (!exists($genes{$id}));
	if (!exists($data{$source})) {
		$data{$source} = '';
		print STDERR "\tFound $source\n";
	}
	my $str = "$sid\t$name\t";
	my $c = 0;
	foreach(@{$genes{$id}}) {
		$str .= "," if ($c > 0);
		$c++;
		$str .= $_;
	}
	$str .= "\n";
	$data{$source} .= $str;
}
close IN;

foreach(keys %data) {
	my $source =$_;
	my $s = $source;
	$s =~ s/\s/_/g;
	open OUT, ">$s.biosystems.genes";
	print OUT $data{$source};
	close OUT;
}
if ($keepFiles == 0) {
	`rm bsid2info biosystems_gene`;	
}
