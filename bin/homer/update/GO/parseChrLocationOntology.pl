#!/usr/bin/perl -w

print STDERR "This file will automatically generate chromosomes.genes from human.description etc.\n";

my $tmpfile = rand() . ".tmp";
`ls -1 *.description > $tmpfile`;
my @files = ();
open IN, $tmpfile;
while (<IN>) {
	chomp;
	s/\r//g;
	push(@files, $_);
}
close IN;
`rm $tmpfile`;

my %chr = ();

foreach(@files) {

	my $file = $_;
	/^(.*?)\.description/;
	my $org = $1;
#print STDERR "$file $org\n";

	open IN, $file;
	while (<IN>) {
		chomp;
		my @line = split /\t/;
		my $gid = $line[0];
		my $loc = $line[7];
		my @details = split / /, $loc;
		my $single = $details[0];
		my $chr1 = "$org chr$loc";
		my $chr2 = "$org chr$single";
		if (!exists($chr{$chr1})) {
			my %a = ();
			$chr{$chr1}=\%a;
		}
		if (!exists($chr{$chr2})) {
			my %a = ();
			$chr{$chr2}=\%a;
		}
		$chr{$chr1}->{$gid}=1;
		$chr{$chr2}->{$gid}=1;
		
	}
	close IN;
}

open OUT, ">chromosome.genes";
foreach(keys %chr) {
	my $chr= $_;
	print OUT "$chr\t$chr\t";
	my $c = 0;
	foreach(keys %{$chr{$chr}}) {
		print OUT "," if ($c >0);
		$c++;
		print OUT "$_";
	}
	print OUT "\n";
}
close OUT;

