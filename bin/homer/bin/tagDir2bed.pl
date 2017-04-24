#!/usr/bin/env perl
use warnings;
use POSIX;

sub printCMD() {
	print STDERR "\n\tUsage: tagDir2bed.pl <tag directory> [options]\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-len <#> (read length to report, default: given sizes in tags.tsv file)\n";
	print STDERR "\t\t-separate (report separate BED reads if there are multiple reads per position)\n";
	print STDERR "\n";
	exit;
}

if (@ARGV < 1) {
	printCMD();
}
my $setLen = -1;
my $sepFlag = 0;
for (my $i=1;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-len') {
		$setLen = $ARGV[++$i];
		print STDERR "\tSetting read lengths to $setLen\n";
	} elsif ($ARGV[$i] eq '-separate') {
		$sepFlag = 1;
		print STDERR "\tWill output separate reads if more than one tag is found per position\n";
	} else {
		print STDERR "What is $ARGV[$i]??\n";
		printCMD();
	}
}
	
my $dname = "Tag";
my $count = 1;

my $file = rand() . ".tmp";

`ls "$ARGV[0]/"*.tags.tsv > "$file"`;	
my @files = ();
open IN, "$file";
while (<IN>) {
	chomp;
	s/\r//g;
	push(@files, $_);
}
close IN;
`rm "$file"`;



for (my $i=0;$i<@files;$i++) {
	my $tagFile = $files[$i];
	open IN, "$tagFile";
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		my $name = $line[0];
		my $chr = $line[1];
		my $pos = $line[2];
		my $dir = $line[3];
		my $v = $line[4];
		my $L = 1;
		if (@line > 5 && $line[5] > -1) {
			$L = $line[5];
		}
		if ($setLen > -1) {
			$L = $setLen;
		}

		if ($v == floor($v)) {
			$v = floor($v);
		}
	
		my $start = $pos;	
		my $end = $start + $L;
		my $d = "+";
		if ($dir == 1) {
			$d = "-";
			$start = $pos - $L;
			$end = $start +  $L;
		} else {
			$start = $start -1;
			$end = $start +  $L;
		}
		
		if ($name eq '') {
			$name = $dname . $count++;
		}
		if ($sepFlag) {
			for (my $j=0;$j<$v;$j++) {
				print "$chr\t$start\t$end\t$name-$j\t1\t$d\n";
			}
		} else {
			print "$chr\t$start\t$end\t$name\t$v\t$d\n";
		}
	
	}
	close IN;
}

