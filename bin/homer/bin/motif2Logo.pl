#!/usr/bin/env perl
use warnings;

if (@ARGV < 1) {
	printCMD();
}

my $motifFile = $ARGV[0];


my $styleFlag = "-S";
my $height = 2.0;
my $widthFactor = 0.75;
my $widthOffset = 3.0;

my $format = "PNG";
my $numSeq = 100;
my $output = $motifFile;

for (my $i=1;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-pdf') {
		$format = "PDF";
	} elsif ($ARGV[$i] eq '-ns') {
		$numSeq = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-bit') {
		$styleFlag = "";
	} elsif ($ARGV[$i] eq '-o') {
		$output = $ARGV[++$i];
	} else {
		print STDERR "!!!Cannot recognize command line option \"$ARGV[$i]\"\n";
		printCMD();
	}
}



my $randnum = rand();
my $tmpfile = $randnum . ".tmp";
my $tmpfile2 = $randnum . ".2.tmp";

open OUT, ">$tmpfile";
open IN, $motifFile;
my $start = 0;
my $length = 0;
while (<IN>) {
	my $og = $_;
	if (/^>/) {
		if ($start == 0) {
			$start = 1;
			print OUT $_;
			next;
		} else {
			#stop
			last;
		}
	}
	if ($start == 1) {
		print OUT $_;
		$length++;
	}
}
close IN;
close OUT;

my $width = $widthFactor * ($widthOffset+$length);
	
`profile2seq.pl "$tmpfile" $numSeq > "$tmpfile2"`;	
`seqlogo -f "$tmpfile2" -F $format $styleFlag -c -o $output -h $height -w $width 2> /dev/null`;
`rm "$tmpfile" "$tmpfile2"`;


sub printCMD {
	print STDERR "\n\tUsage: motif2Logo.pl <motif file> [options]\n";
	print STDERR "\t\tBy default produces \"motif file\".png\n";
	print STDERR "\t\t!! Only processes the first motif in the motif file !!\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-pdf (create a pdf: \"motif file\".pdf: default creates a PNG image)\n";
	print STDERR "\t\t-ns <#> (Number of sequences to feed seqlogo: default 100)\n";
	print STDERR "\t\t-bit (scale logo by information content: default scales nucleotides to probability)\n";
	print STDERR "\t\t-o <OutputPrefix> (prefix of output file, i.e. outputprefix.png : default, motif file)\n";
	print STDERR "\n";
	exit;
}
