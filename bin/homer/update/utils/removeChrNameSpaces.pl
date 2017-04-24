#!/usr/bin/perl -w

for (my $i=0;$i<@ARGV;$i++) {
	my $file =$ARGV[$i];
	print STDERR "\tProcessing $file: ";
	open OUT, ">.tmp";
	open IN, $ARGV[$i];
	while (<IN>) {
		chomp;	
		s/\r//g;
		if (/^>/) {
			print STDERR "$_ -> ";
			s/\s.*$//;
			s/\..*$//;
			s/\-.*$//;
			print STDERR "$_\n";
			print OUT "$_\n";
		} else {
			print OUT "$_\n";
		}
	}
	close IN;
	close OUT;
	`mv .tmp $file`;
}

