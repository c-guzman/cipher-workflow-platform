#!/usr/bin/perl -w
if (@ARGV < 1) {
	print STDERR "\n\tUsage: transfac2homer.pl <transfac matrix.dat file> > output.homer.motifs\n\n";
	exit;
}
my $id = "";
my $acc = "";
my $name = "";
my @m = ();
my $matrix = \@m;

open IN, $ARGV[0];
while (<IN>) {
	chomp;
	s/\r//g;
	if (/^\/\//) {
		if ($id ne '') {
			printMatrix($id,$acc,$name,$matrix);
		}
		$id = "";
		$acc = "";
		$name = "";
		my @m = ();
		$matrix = \@m;
	} elsif (/^AC\s+(.+)/) {
		$acc = $1;
	} elsif (/^NA\s+(.+)/) {
		$name = $1;
	} elsif (/^ID\s+(.+)/) {
		$id = $1;
	} elsif (s/^\d\d\s+//) {
		my @freq = split /\s+/;
		push(@$matrix, \@freq);
	}
}
close IN;

if ($id ne "")  {
	printMatrix($id,$acc,$name,$matrix);
}

sub printMatrix {
	my ($id, $acc,$name,$matrix) = @_;

	my $con = '';
	foreach(@$matrix) {
		$con .= $_->[4];
	}
	print ">$con\t$name($id/$acc)\t0\n";
	foreach(@$matrix) {
		print "$_->[0]";
		for (my $i=1;$i<4;$i++){ 
			print "\t$_->[$i]";
		}
		print "\n";
	}
	
}
