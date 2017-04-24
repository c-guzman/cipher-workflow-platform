#!/usr/bin/perl -w

foreach(@ARGV) {
	my $file = $_;
	my $person = $file;
	$person =~ s/ready-//;
	$person =~ s/\.txt//;

	open IN, $file;
	my @lines = ();
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line = split /\s/;
		push(@lines, \@line);
	}
	close IN;
	for (my $i=0;$i<@lines;$i+=6) {
		my $name = $lines[$i][1];
		my $tfname = "$name/dmmpmm($person)/fly";
		print ">$tfname\t$tfname\t0\n";
		for (my $j=0;$j<@{$lines[$i+1]};$j++) {
			print "$lines[$i+1][$j]";
			for (my $k=2;$k<5;$k++) {
				print "\t$lines[$i+$k][$j]";
			}
			print "\n";
		}
	}
}
		
