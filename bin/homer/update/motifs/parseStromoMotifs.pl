#!/usr/bin/perl -w

my @matrix = ();

my $person = "Harbison";
my $consensus = "";
open IN, $ARGV[0];
while (<IN>) {
	chomp;
	s/\r//g;
	if (/Source:\s*(.*?)$/) {
		my $name = $1;
		$name =~ s/\s+\d+\s+/ /g;
		$name =~ s/\s+/\//g;
		print ">$consensus\t$name($person)/Yeast\t0\n";

		my $m = fixMatrix(\@matrix);
		foreach(@$m) {
			my $c = 0;
			foreach(@$_) {
				$c++;
				print "\t" if ($c > 1);
				my $v = sprintf("%.3f",$_);
				print "$v";
			}
			print "\n";
		}


		@matrix = ();
	} elsif (/^Log.*\s\d+\s(.*?)\s\(/) {
		$consensus = $1;
		@matrix = ();
	} elsif (/^#A/ || /^#C/ || /^#T/ || /^#G/) {
		my @p = split /\s+/;
		shift @p;
		push(@matrix, \@p);
	}
}
close IN;


sub fixMatrix {
	my ($matrix) = @_;

	my @new = ();
	my $len = scalar(@{$matrix->[0]});
print STDERR "Len = $len\n";
	for (my $i=0;$i<$len;$i++) {
		my $total = 0;
		my $min = 0;
		my @d = ();
		for (my $j=0;$j<4;$j++) {
			$matrix[$j][$i] = exp($matrix[$j][$i]);
			if ($matrix[$j][$i] < $min) {
				$min = $matrix[$j][$i];
			}
		}
		$total = 0;
		for (my $j=0;$j<4;$j++) {
			$matrix[$j][$i] = $matrix[$j][$i]-$min;
			$total += $matrix[$j][$i];
		}
		for (my $j=0;$j<4;$j++) {
			$matrix[$j][$i] = $matrix[$j][$i]/$total;
		}
		my @n = ($matrix[0][$i],$matrix[1][$i],$matrix[3][$i],$matrix[2][$i]);
		push(@new,\@n);
		next;
	}
	return \@new;
}
