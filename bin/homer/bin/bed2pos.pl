#!/usr/bin/env perl
use warnings;




if (@ARGV < 1) {
	print STDERR "\n\tUsage: bed2pos.pl [options] <BED file>\n";
	print STDERR "\n\tThis outputs a position/peak file to stdout\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-check (Checks if the file is already peak/pos formatted)\n";
	print STDERR "\t\t-unique (Make peaks names unique by adding numbers to replicate names)\n";
	print STDERR "\t\t-o <filename> (Send output to this file, default: stdout)\n";
	print STDERR "\t\t-pos (Send output to file with same name as input file with *.pos extension)\n";
	print STDERR "\n";
	exit;
}
my $check = 0;
my $outputFile = '';
my $extensionFlag = 0;
my $inputFile = "";
my $makeUnique = 0;

for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-check') {
		$check = 1;
	} elsif ($ARGV[$i] eq '-unique') {
		$makeUnique = 1;
	} elsif ($ARGV[$i] eq '-o') {
		$outputFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-pos') {
		$extensionFlag = 1;
	} else {
		$inputFile = $ARGV[$i];	
	}
}
if ($inputFile eq '') {
	print STDERR "!!! Missing input BED file!!!\n";
	exit;
}
if ($extensionFlag) {
	$outputFile = $inputFile;
	$outputFile =~ s/(\..*?)$//;
	$outputFile .= ".pos";
}

my $numDuplicates = 0;
my %names = ();
my %dupCount = ();

my $filePtr = *STDOUT;
if ($outputFile ne '') {
	print STDERR "\tOutput File: $outputFile\n";
	open OUT, ">$outputFile";
	$filePtr = *OUT;
}
	
my $bedFormat = 0;
my $posFormat = 0;

	
open IN, $inputFile or die "!!! Could not open file $inputFile !!!\n";
my $ID = 1;
while (<IN>) {
	chomp;
	s/\r//g;
	my $og = $_;
	my @line = split /\t/;
	next if (@line < 3);
	next if ($line[0] =~ /^\s*browser/);
	next if ($line[0] =~ /^\s*track/);
	next if ($line[0] =~ /^\s*#/);
	my $pflag = 0;
	my $name = $ID;
	my $start = 0;
	my $end = 0;
	my $dir = 0;
	my $v = 0;
	my $chr = '';

	if ($check) {
		if (@line > 4) {
			if ($line[2] =~ /^\d+$/ && $line[3] =~ /^\d+$/ && $line[4] =~ /^[01+-\.]$/) {
				$pflag = 1;
				$name = $line[0];
			} elsif ($line[1] =~ /^\d+$/ && $line[2] =~ /^\d+$/) {
				# probably a BED line
			} else {
				#probably a header line...
				next;
			}
		}
	}



	if ($pflag == 0) {
		#parse bed file

		$dir = 0;
		$v = 0;
		$chr= $line[0];
		$start= $line[1]+1;
		$end= $line[2];
	
		if (@line > 3) {
			if ($line[3] eq '-') {
				$dir = 1;
			} elsif ($line[3] eq '+') {
				$dir = 0;
			} elsif ($line[3] eq '.') {
				$dir = 0;
			} else {
				$name = $line[3];
			}
		} else {
			$ID++;
		}

		if (@line > 4) {
			if ($line[4] =~ /^[\d\.\-\+\e]+$/) {
				$v = $line[4];
			}
		}
		if (@line > 5) {
			if ($line[5] eq '+') {
				$dir = 0;
			} elsif ($line[5] eq '.') {
				$dir = 0;
			} elsif ($line[5] eq '-') {
				$dir = 1;
			}
		}
	} 

	if (exists($names{$name})) {
		$numDuplicates++;
		if (exists($dupCount{$name})) {
			$dupCount{$name}++;
		} else {
			$dupCount{$name} = 2;
		}
 		if ($makeUnique) {
			$name = $name . "-" . $dupCount{$name};
			$names{$name}=1;
		}
	} else {
		$names{$name}=1;
	}

	if ($pflag == 1) {
		print $filePtr "$name";
		for (my $i=1;$i<@line;$i++) {
			print $filePtr "\t$line[$i]";
		}
		print $filePtr "\n";
		$posFormat++;
	} else {
		print $filePtr "$name\t$chr\t$start\t$end\t$dir\t$v\n";
		$bedFormat++;
	}
}
close IN;
if ($check) {
	print STDERR "\tPeak/BED file conversion summary:\n";
	print STDERR "\t\tBED/Header formatted lines: $bedFormat\n";
	print STDERR "\t\tpeakfile formatted lines: $posFormat\n";
}
if ($makeUnique) {
	print STDERR "\t\tDuplicated Peak IDs: $numDuplicates\n";
}
