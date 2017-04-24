#!/usr/bin/env perl
use warnings;



if (@ARGV < 3) {
	print STDERR "\n\tUsage: batchFindMotifsGenome.pl [genome] [options...] -d <TagDirectory1> [TagDirectory 2] ...\n";
	print STDERR "\tUsage: batchFindMotifsGenome.pl [genome] [options...] -f <peak/BED file1> [peak/BED file2] ...\n";
	print STDERR "\n\t\t-dist <#> (Will only analyze promoter-distal regions ># away from TSS)\n";
	print STDERR "\n";
	exit;
}

my $genome = $ARGV[0];
my $options = "";
my $dist = '';
my @files = ();
my @ffiles = ();
for (my $i=1;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-d') {
		for (my $j=$i+1;$j<@ARGV;$j++) {
			push(@files, $ARGV[$j]);
		}
		last;
	} elsif ($ARGV[$i] eq '-f') {
		for (my $j=$i+1;$j<@ARGV;$j++) {
			push(@ffiles, $ARGV[$j]);
		}
		last;
	} elsif ($ARGV[$i] eq '-dist') {
		$dist = $ARGV[++$i];
	} else {
		$options .= " " . $ARGV[$i];
	}
}

foreach(@files) {
	my $dir = $_;
	my $outdir = $_ . "/" . "Motifs-BatchOutput/";
	if ($dist ne '') {
		$outdir = $_ . "/" . "Motifs-BatchOutput.distal/";
	}
		

	my $inputFile = "$dir/peaks.txt";
	$inputFile = "$dir/regions.txt" unless (-e $inputFile);
	$inputFile = "$dir/tss.txt" unless (-e $inputFile);
	unless (-e $inputFile) {
		print STDERR "!!! Error, could not find peak file in directory $dir!!!\n";
		next;
	}
	my $inputFile2 = $inputFile . ".dist.txt";

	if ($dist ne '') {
		`getDistalPeaks.pl $inputFile $genome -d $dist > $inputFile2`;
		$inputFile = $inputFile2;
	}

	print STDERR "findMotifsGenome.pl \"$inputFile\" $genome \"$outdir\" $options\n";
	`findMotifsGenome.pl "$inputFile" $genome "$outdir" $options`;
	
}
foreach(@ffiles) {
	my $inputFile = $_;
	my $outdir = "Motifs-" . $inputFile;
	my $inputFile2 = $inputFile . ".dist.txt";
	if ($dist ne '') {
		$outdir = "Motifs-" . $_ . ".distal/";
		`getDistalPeaks.pl $inputFile $genome -d $dist > $inputFile2`;
		$inputFile = $inputFile2;
	}
	print STDERR "findMotifsGenome.pl \"$_\" $genome \"$outdir\" $options\n";
	`findMotifsGenome.pl "$_" $genome "$outdir" $options`;
}
