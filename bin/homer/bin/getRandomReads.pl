#!/usr/bin/env perl
use warnings;


if (@ARGV < 2) {
	print STDERR "\n\tUsage: getRandomReads.pl <tag directory> <# of reads>\n";
	print STDERR "\tThis program will send new tag file to stdout.\n";
	print STDERR "\tMake new directory with -t option\n";
	print STDERR "\n";
	exit;
}
my $tagDir = $ARGV[0];
my $count = $ARGV[1];

$peFlag = 0;
open IN, "$ARGV[0]/tagInfo.txt";
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	if ($line[0] =~ /^genome/) {
		$totalReads = $line[2];
	}
	if ($line[0] eq 'pairedEnd=true') {
		$peFlag = 1;
	}
}
close IN;

print STDERR "\tTotal Reads in $tagDir: $totalReads\n";
if ($count > $totalReads) {
	print STDERR "\tHold on bro: $count > $totalReads\n";
	print STDERR "\tTry requesting less reads...\n";
	exit;
}
if ($peFlag==1) {
	print STDERR "\tSampling from a Paired-end directory...\n";
}
if ($totalReads < 1) {
	print STDERR "\ttotalReads = $totalReads!! Something is probably wrong...\n";
	exit;
}
my $ratio = $count/$totalReads;
my $randID = rand();
my $tmpFile = "$randID.tmp";
`ls -1 "$tagDir"/*.tags.tsv > "$tmpFile"`;
open IN, $tmpFile;
my @files = ();
while (<IN>) {
	chomp;
	push(@files, $_);
}
close IN;
`rm "$tmpFile"`;

for (my $i=0;$i<@files;$i++) {
	open IN, $files[$i];
	print STDERR "\t\t$files[$i]\n";
	while (<IN>) {
		chomp;
		my @line = split /\t/;
		my $v = $line[4];
		if ($peFlag) {
			my $cmp = $line[1] cmp $line[6];
			next if ($cmp < 0);
			next if ($cmp == 0 && $line[2] < $line[7]);
		}
		for (my $j=0;$j<$v;$j++) {
			my $r = rand();
			if ($r < $ratio) {
				if ($peFlag==1) {
					print "\t$line[1]\t$line[2]\t$line[3]\t1\t$line[5]\t$line[6]\t$line[7]\t$line[8]\t$line[9]\n";
					print "\t$line[6]\t$line[7]\t$line[8]\t1\t$line[9]\t$line[1]\t$line[2]\t$line[3]\t$line[5]\n";
				} else {
					print "\t$line[1]\t$line[2]\t$line[3]\t1";
					for (my $i=5;$i<@line;$i++) {
						print "\t$line[$i]";
					}
					print "\n";
				}
			}
		}
	}
	close IN;
}
	
