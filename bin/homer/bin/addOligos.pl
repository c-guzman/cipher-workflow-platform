#!/usr/bin/env perl
use warnings;



if (@ARGV < 2) {
	print STDERR "<oligo file1> <oligo file 2>\n";
	print STDERR "joins files that contain oligos, looking for reverse opposites\n";
	exit;
}

open IN, $ARGV[1];
while (<IN>) {
	chomp;
	s/\r//g;
	my @line= split /\t/;
	my $id = shift @line;
	my $rev = revopp($id);
	$data{$id} = $_;
	$data{$rev} = $_;
}

my $maxCol = 0;
open IN, $ARGV[0];
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	if ($maxCol < @line) {
		$maxCol = @line;
	}
}
close IN;
print STDERR "Max Col = $maxCol\n";

open IN, $ARGV[0];
while (<IN>) {
	chomp;
	s/\r//g;
	print $_;
	my @line = split /\t/;
	if (@line < $maxCol) {
		for (my $i=@line;$i<$maxCol;$i++) {
			print "\t";
		}
	}
	
	if (exists($data{$line[0]})) {
		print "\t$data{$line[0]}";
	}
	print "\n";
	

}
close IN;


sub revopp {
	my ($seq)  = @_;
	$seq = reverse($seq);
	$seq =~ s/A/X/g;
	$seq =~ s/T/A/g;
	$seq =~ s/X/T/g;
	$seq =~ s/C/X/g;
	$seq =~ s/G/C/g;
	$seq =~ s/X/G/g;

	return $seq;
}
