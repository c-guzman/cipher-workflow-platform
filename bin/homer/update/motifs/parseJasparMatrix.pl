#!/usr/bin/perl -w
print STDERR "\n\tParsing JASPAR matrix file $ARGV[0] into a HOMER style motif file\n\n";

open IN, $ARGV[0] or die "Could not open file $ARGV[0]\n";

my @data = ();
while (<IN>) {
	chomp;
	if (/^>/) {
		s/^>//;
		s/ /_/g;
		my @line = split /\t/;
		my $n = $line[1] . "/" . $line[0];
		push(@data, $n);
	} elsif ($_ eq '') {
		next;
	} else {
		s/^.+\[\s*//;
		s/\s*\].*$//;
		my @d = split /\s+/;
		push(@data, \@d);
	}
}
close IN;


for (my $i=0;$i<@data;$i+=5) {
	my $name= $data[$i+0];
	next if ($name =~ /^CN/);
	next if ($name =~ /^PF/);
	print ">$name\t$name/Jaspar\t0\n";
	my $len = scalar(@{$data[$i+1]});
	for (my $j=0;$j<$len;$j++) {
		my $v = $data[$i+1][$j];
		print "$v";
		for (my $k=1;$k<4;$k++) {
			my $vv = $data[$i+1+$k][$j];
			print "\t$vv";
		}
		print "\n";
	}
}
