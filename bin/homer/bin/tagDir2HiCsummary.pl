#!/usr/bin/env perl
use warnings;
use POSIX;

sub printCMD() {
	print STDERR "\n\tUsage: tagDir2HiCsummary.pl <tag directory>\n";
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
	
my $dname = "ID";
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
		next if (@line < 10);
		my $name = $line[0];
		my $chr = $line[1];
		my $pos = $line[2];
		my $dir = $line[3];
		my $v = $line[4];
		my $len = $line[5];
		my $chr2 = $line[6];
		my $pos2 = $line[7];
		my $dir2 = $line[8];
		my $len2 = $line[9];


		my $cmp1 = $chr cmp $chr2;
		next if ($cmp1 > 0);
		if ($cmp1 == 0) {
			my $cmp2 = $pos <=> $pos2;
			next if ($cmp2 > 0);
		}

		if ($v == floor($v)) {
			$v = floor($v);
		}
	
		if ($dir == 1) {
			$dir = "-";
		} else {
			$dir = '+';
		}
		if ($dir2 == 1) {
			$dir2 = "-";
		} else {
			$dir2 = '+';
		}
		print "$dname$count\t$chr\t$pos\t$dir\t$chr2\t$pos2\t$dir2\n";
		$count++;

	
	}
	close IN;
}

