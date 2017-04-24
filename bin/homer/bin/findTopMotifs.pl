#!/usr/bin/env perl
use warnings;


if (@ARGV < 3) {
	printCMD();
}
my $peakFile = $ARGV[0];
my $genome = $ARGV[1];
my $searchGenome = $ARGV[2];
my $size = 200;
my $searchSize = '';

my $motifFindingOptions = '';
my $prefix = $peakFile;
my $numTopMotifs = 5;

for (my $i=3;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-top') {
		$numTopMotifs = $ARGV[++$i];
		print STDERR "Number of top motifs to find set to $numTopMotifs\n";;
		next;
	} elsif ($ARGV[$i] eq '-prefix') {
		$prefix = $ARGV[++$i];
		print STDERR "Top Motifs will be stored with prefix $prefix\n";
		next;
	} elsif ($ARGV[$i] eq '-searchSize') {
		$searchSize = $ARGV[++$i];
		print STDERR "Fragment Search size = $searchSize\n";
		next;
	} elsif ($ARGV[$i] eq '-size') {
		$size = $ARGV[++$i];
		print STDERR "Fragment size = $size\n";
		next;
	} else {
		$motifFindingOptions .= " " . $ARGV[$i];
		next;
	}
}

if ($searchSize eq '') {
	$searchSize = $size;
}

my $rand = rand();
my $tmpfile = $rand . ".tmp";
my $tmpfile2 = $rand . ".2.tmp";
my $tmpfile3 = $rand . ".3.tmp";
my $tmpDir = $rand . ".tmpDir";
`cp $peakFile $tmpfile`;
my $currentPeakFile = $tmpfile;
for (my $i=0;$i<$numTopMotifs;$i++) {
	`findMotifsGenome.pl $tmpfile $genome $tmpDir -noknown -nocheck -S 5 -size $size $motifFindingOptions`;
	`cat $tmpDir/homerMotifs.motifs* > $tmpfile2`;

	my ($pvalue, $motifStr) = openMotifFile($tmpfile2);
	my $iter = $i+1;
	my $bestmotif = $prefix . ".Top" . $iter . ".motif";
	open OUT, ">$bestmotif";
	print OUT $motifStr;
	close OUT;

	`annotatePeaks.pl $tmpfile $searchGenome -size $searchSize -center $bestmotif > $tmpfile2`;
	`filterListBy.pl $tmpfile $tmpfile2 > $tmpfile3`;
	`rm $tmpfile2`;
	`rm -r $tmpDir`;
	`mv $tmpfile3 $tmpfile`;
}
`rm $tmpfile`;

sub printCMD {
	print STDERR "\n\tUsage: ./getTopMotifs.txt <peak file> <genome> <search genome> [options]\n";
	print STDERR "\n\t\tThis program finds the top motif, removes the peaks that contain it, \n";
	print STDERR "\t\tthen repeats this for the -top <#> times - motifs are named prefix.top1.motif...\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-top <#> (number of times to find top motif default:5)\n";
	print STDERR "\t\t-size <#> (size of fragment | 200)\n";
	print STDERR "\t\t-searchSize <#> (size of fragment to search to remove for future rounds | 200)\n";
	print STDERR "\t\t-prefix <file name prefix> (Name of motif files: prefix.top1.motif... | peakfile.top1.motif)\n";
	print STDERR "\t\tALL other options will be passed to findMotifsGenome for motif finding!!\n";
	print STDERR "\n";
	exit;
}


sub openMotifFile {
	my ($file) = @_;
	my $motifString = $_;
	my $pvalue = 10000;
	open IN, $file;
	while (<IN>) {
		my $og = $_;
		if (/^>/) {
			chomp;
			$good = 0;
			my @line = split /\t/;
			my $p = $line[3];
			if ($p < $pvalue) {
				$motifString = $og;
				$good = 1;
				$pvalue = $p;
			}
			next;
		}
		$motifString .= $og if ($good ==1);
	}
	close IN;
	return ($pvalue, $motifString);
}
