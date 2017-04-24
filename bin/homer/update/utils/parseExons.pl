#!/usr/bin/perl -w
if (@ARGV < 1) {
	print STDERR "<refseq alignment> [options]\n";
	exit;
}
open OUT3, ">3p.exonJunc.tsv";
open OUT5, ">5p.exonJunc.tsv";
open CDS, ">ATG.tsv";

my %lens = ();
my %done = ();
open IN, $ARGV[0];
while (<IN>) {
	chomp;
	my @line = split /\t/;
	my $id = $line[10];
	next if (exists($done{$id}));
	my $chr = $line[14];

	my $tss = $line[16]+1;
	my $tes = $line[17];
	my @len = split /\,/, $line[19];
	foreach(@len) {
		$len{$_}++;
	}
	my @starts = split /\,/, $line[21];
	if ($line[9] eq '+') {
		for (my $i=0;$i<@starts;$i++) {
			my $p5 = $starts[$i]+1;
			my $p3 = $starts[$i]+$len[$i];
			my $s = $p5-100;
			my $e = $p5+100;
			my $c = $i+1;
			if ($i!=0) {
				print OUT5 "$id-$c-5p\t$chr\t$s\t$e\t0\n";
			}
			$s = $p3-100;
			$e = $p3+100;
			if ($i!=@line-1) {
				print OUT3 "$id-$c-3p\t$chr\t$s\t$e\t0\n";
			}
		}
	}
	if ($line[9] eq '-') {
		for (my $i=0;$i<@starts;$i++) {
			my $p5 = $starts[$i]+1;
			my $p3 = $starts[$i]+$len[$i];
			my $s = $p5-100;
			my $e = $p5+100;
			my $c = $i+1;
			if ($i!=0) {
				print OUT3 "$id-$c-3p\t$chr\t$s\t$e\t1\n";
			}
			$s = $p3-100;
			$e = $p3+100;
			if ($i!=@line-1) {
				print OUT5 "$id-$c-5p\t$chr\t$s\t$e\t1\n";
			}
		}
	}
	$done{$id} = 1;
}
close IN;

my @lens = sort {$a <=> $b} keys %len;
for (my $i=0;$i<=$lens[@lens-1];$i++) {
	if (exists($len{$i})) {
		print "$i\t$len{$i}\n";
	} else {
		print "$i\t0\n";
	}
}
