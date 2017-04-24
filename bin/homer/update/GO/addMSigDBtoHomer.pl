#!/usr/bin/perl -w

if (@ARGV < 3) {
	print STDERR " <msigdb.v2.5.symbols.gmt file> <human2gene.tsv> <homologene.dat>\n";
	print STDERR " Uses human gene symbols\n";
	exit;
}

my %conv = ();
my %groups = ();
my %lookup= ();

open IN, $ARGV[0];
while (<IN>) {
	chomp;	
	s/\r//g;
	my @line = split /\t/;
	my $name = shift @line;
	my $na = shift @line;
	$groups{$name} = \@line;
	foreach(@line) {
		$lookup{$_}=1;
	}
}

open IN, $ARGV[1];
while (<IN>) {
	chomp;
	my @line = split /\t/;	
	if (exists($lookup{$line[0]})) {
		next if ($line[1] eq '-');
		$gid{$line[0]} = $line[1];
		my @a = ();
		$conv{$line[1]} = \@a;
	}
}
close IN;

open IN, $ARGV[2];
my @current = ();
my $id = "";
my $mapid = "";
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	if ($line[0] ne $id) {
		if ($mapid ne '') {
			push(@{$conv{$mapid}},@current);
		}	
		$mapid = '';
		$id = $line[0];
		@current = ();
	}
	push(@current, $line[2]);
	if (exists($conv{$line[2]})) {
		$mapid = $line[2];
	}
}
close IN;

foreach(keys %groups) {
	print "$_\t$_\t";
	my $g = $_;
	my $c = 0;
	foreach(@{$groups{$g}}) {
		next if (!exists($gid{$_}));
		my $gid = $gid{$_};
		next if (!exists($conv{$gid}));
		foreach(@{$conv{$gid}}) {
			print "," if ($c > 0);
			$c++;
			print $_;
		}
	}
	print "\n";
}


