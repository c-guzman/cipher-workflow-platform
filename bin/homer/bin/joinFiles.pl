#!/usr/bin/env perl
use warnings;

if (@ARGV < 2) {
	print STDERR "[-h] <file1> <file2> <file3> ...\n";
	exit;
}
my $program = "addData.pl";
if ($ARGV[0] eq '-h') {
	$program  = "addDataHeader.pl";
	shift @ARGV;
} 

my $tmp = rand().  ".tmp";
my $tmp2 = rand().  ".tmp";
`cp $ARGV[0] $tmp2`;
for (my $i=1;$i<@ARGV;$i++) {
	`$program $tmp2 $ARGV[$i] > $tmp`;
	`mv $tmp $tmp2`;
}
open IN, $tmp2;
while (<IN> ){
	print $_;
}
close IN;
`rm $tmp $tmp2`;
