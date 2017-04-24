#!/usr/bin/perl -w

my %info = ();

open IN, $ARGV[0];
my $c= 0;
while (<IN>) {
	$c++;
	chomp;
	s/\r//g;
	my @line = split /\t/;
	next if ($c ==1);
	$info{$line[0]} = \@line;
}
close IN;

foreach(keys %info) {
	my $id = $_;
	my $file = "top10align_motifs/" . $_ . "_top10align_pfm.txt";
	open IN, $file or die "Couldn't open $file\n";
	my $c = 0;
	while (<IN>) {
		$c++;
		chomp;
		s/\r//g;
		my @line = split /\t/;
		if ($c==1) {
			my $name= $info{$id}->[2] . "(" . $info{$id}->[5] . ")/" . $info{$id}->[3] . "-" . $id . "-PBM/HughesRNA";
			print ">$id\t$name\t0\n";
		} else {
			my $v = sprintf("%.3f",$line[1]);
			print "$v";
			for (my $i=2;$i<@line;$i++) {
				$v = sprintf("%.3f",$line[$i]);
				print "\t$v";
			}
			print "\n";
		}
	}
}
		
