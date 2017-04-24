#!/usr/bin/env perl
use warnings;



if (@ARGV < 3) {
	print STDERR "\n\tUsage: batchFindMotifs.pl [promoter set] [options...] -f list1.txt list2.txt ...\n";
	print STDERR "\n";
	exit;
}

my $promoterSet = $ARGV[0];
my $options = "";
my @files = ();
for (my $i=1;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-f') {
		for (my $j=$i+1;$j<@ARGV;$j++) {
			push(@files, $ARGV[$j]);
		}
		last;
	} else {
		$options .= " " . $ARGV[$i];
	}
}

foreach(@files) {
	my $dir = $_;
	$dir =~ s/\.txt$//;
	$dir =~ s/\.bed$//;
	$dir = "Motifs-$dir";
	print STDERR "findMotifs.pl \"$_\" $promoterSet \"$dir\" $options\n";
	`findMotifs.pl "$_" $promoterSet "$dir" $options`;
}
