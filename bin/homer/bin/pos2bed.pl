#!/usr/bin/env perl
use warnings;
#
use POSIX;

if (@ARGV < 1) {
	print STDERR "\n\tUsage: pos2bed.pl [options] <peak/pos file>\n";
	print STDERR "\n\tThis will output a BED-format file to stdout\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-o <filename> (Output to file)\n";
	print STDERR "\t\t-bed (Output to file with same name as input with *.bed extension)\n";
	print STDERR "\t\t-track <name> (Include track line with name for uploading to UCSC Genome Browser)\n";
	print STDERR "\t\t-5 (Set 5th column to the value 1 instead of value in 6th column of pos file)\n";
	print STDERR "\t\t-float (Allow the 5th column to be a floating point number, default: integer)\n";
	print STDERR "\t\t-color strand (color strands red and blue, will also add a track line to file)\n";
	print STDERR "\n";
	exit;
}
my $inputFile = "";
my $outputFile = "";
my $oneFlag = 0;
my $trackname = "";
my $bedFlag = 0;
my $integerFlag = 1;
my $colorFlag = 0;
my $colorType = '';

for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-o') {
		$outputFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-bed') {
		$bedFlag = 1;
	} elsif ($ARGV[$i] eq '-track') {
		$trackname = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-float') {
		$integerFlag = 0;
	} elsif ($ARGV[$i] eq '-color') {
		$colorFlag = 1;
		$colorType = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-5') {
		$oneFlag = 1;
	} else {
		$inputFile = $ARGV[$i];
	}
}
if ($inputFile eq '') {
	print STDERR "!!! Need to specify an input file!!!\n";
	exit;
}
if ($bedFlag) {
	$outputFile = $inputFile;
	$outputFile =~ s/(\..*?)$//;
	$outputFile .= ".bed";
}

my $filePtr = *STDOUT;
if ($outputFile ne '') {
	print STDERR "\tOutput File: $outputFile\n";
	open OUT, ">$outputFile";
	$filePtr = *OUT;
}

if ($trackname ne '' || $colorFlag==1) {
	if ($trackname eq '') {
		$trackname = $inputFile;
	}
	print $filePtr "track name=\"$trackname\"";
	if ($colorFlag) {
		print $filePtr " itemRgb=\"On\"";
	}
	print $filePtr "\n";
}

	

my $total = 0;

open IN, $inputFile or die "!!! Could not open input file: $inputFile!!!!\n";
while (<IN>) {
	chomp;
	s/\r//g;
	if (/^#/) {
		print $filePtr "$_\n";
		next;
	}
	my @line = split /\t/;
	my $dir = '+';
	my $v = 0;
	next if (@line < 5);
	next if ($line[2] =~ /^[^\d\.\-]/);

	my $chr = $line[1];
	my $start = $line[2]-1;
	my $end = $line[3];

	if (@line > 5) {
		$v = $line[5];
	}
	if ($line[4] eq '1' || $line[4] eq '-') {
		$dir = '-';
	}
	if ($oneFlag ==1) {
		$v = 1;
	}
	if ($v =~ /^[\de\+\-]+$/) {
		if ($integerFlag) {
			$v = floor($v+0.5);
		}
	} else {
		$v = 1;
	}
	$start = 0 if ($start < 0);
	$end = 0 if ($end < 0);
	$total++;
	print $filePtr "$chr\t$start\t$end\t$line[0]\t$v\t$dir";

	if ($colorFlag && $colorType eq 'strand') {
		my $color = "0,0,255";
		if ($dir eq '-') {
			$color = "255,0,0";
		}
		print $filePtr "\t$start\t$end\t$color";

	}
	print $filePtr "\n";
}
close IN;
print STDERR "\n\tConverted $total peaks total\n\n";

