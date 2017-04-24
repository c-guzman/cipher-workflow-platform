#!/usr/bin/perl -w

if (@ARGV < 1) {
	print STDERR "<atha table>\n";
	print STDERR "copy and paste table from http://www.athamap.de/documentation_matrix_based.php into a file\n";
	exit;
}


open IN, $ARGV[0];
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	foreach(@line) {
		s/\^\s+//;
		s/\s+$//;
	}
	my $name = $line[0];
	next if ($line[0] eq 'Factor');
print STDERR "Name = \"$name\"\n";
	my $family = $line[1];
	my $species = $line[2];
	my $urlName = $name;
	$urlName =~ s/\(/\%28/g;
	$urlName =~ s/\)/\%29/g;
	`wget -O tmp.txt http://www.athamap.de/detail.php?name=$urlName`;
	open IN2, "tmp.txt";
	my $in = "";
	while (<IN2>) {
		$in .= $_;
	}
	close IN2;
		
#print STDERR "\nin=$in\n\n";
	my @lines = split /\n/,$in;
	my @matrix = ();
	foreach(@lines) {
		chomp;
		my $og = $_;
		if (s/^.*[ACGT] \|\s+//) {
			s/\<.+$//g;
			my @counts = split /\s+/;
			push(@matrix, \@counts);
		} else {
		}
	}
	print ">$name\t$name($family)/$species/AthaMap\t0\n";
	for (my $i=0;$i<@{$matrix[0]};$i++) {
		print "$matrix[0][$i]";
		for (my $j=1;$j<4;$j++) { 
			print "\t$matrix[$j][$i]";
		}
		print "\n";
	}
	`rm tmp.txt`;
}
close IN;
