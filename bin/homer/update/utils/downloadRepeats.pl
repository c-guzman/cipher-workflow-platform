#!/usr/bin/perl -w

my $maxChr=22;
if (@ARGV < 1) {
	print STDERR "\n\tUsage: downloadRepeats.pl <url base>\n";
	print STDERR "\t\tWill download chr1 through chr$maxChr (and chrX, chrY, chrM)\n";
	print STDERR "\tExample of URL base: http://hgdownload.cse.ucsc.edu/goldenPath/mm8/database/\n";
	print STDERR "\n";
	exit;
}
my $base = $ARGV[0];
my @chr = ();
for (my $i=1;$i<=$maxChr;$i++) {
	push(@chr, "chr" . $i);
}
push(@chr, "chrX", "chrY", "chrM");
foreach(@chr) {
	my $file = $base . "/" . $_ . "_rmsk.txt.gz";
	`wget $file`;
}
