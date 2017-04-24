#!/usr/bin/env perl
use warnings;

if (@ARGV < 1) {
	print STDERR "<result directory> [output zip file name, default=directoryName.zip]\n";
	exit;
}

my $dir = $ARGV[0];
my $output = $dir;
if (@ARGV > 1) {
	$output = $ARGV[0];
} else {
	$output =~ s/\///g;
	$output =~ s/\.//g;
	$output =~ s/\~//g;
	$output .= ".zip";
}

`zip -r $output $dir/ -x $dir/*.tags.tsv $dir/putativePeaks.txt $dir/*target* $dir/*seq.tsv $dir/*adj`; 


