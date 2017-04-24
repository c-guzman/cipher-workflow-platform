#!/usr/bin/env perl
use warnings;

my $homeDir = "/bioinformatics/homer";


use POSIX;


my $check = `which bedToBigBed`;
chomp $check;
unless ($check =~ /\//) {
	print STDERR "!!! $check !!!\n";
	print STDERR "Could not find bedToBigBed (download from UCSC Genome browser and place in executable path)\n";
}

	

if (@ARGV < 3) {
	print STDERR "\n\tMake known motif bigBed track:\n";
	print STDERR "\n\tUsage: makeBigBedMotifTrack.pl <track name> <motif file> <genome>\n";
	print STDERR "\n\tWill create files <track name>.<genome>.track.txt and <track name>.<genome>.bigBed";
	print STDERR "\n\tbedToBigBed found: $check\n";
	print STDERR "\n";
	exit;
}
exit;
my $name = $ARGV[0];
my $knownMotifs = $ARGV[1];
my $genome = $ARGV[2];
my $maskFlag = "";
if ($genome =~ s/r$//) {
	$maskFlag = " -mask ";
}
my $rand = rand();
my $tmpfile = $rand . ".tmp";
my $tmpfile2 = $rand . ".2.tmp";
my $tmpfile3 = $rand . ".3.tmp";
`scanMotifGenomeWide.pl $knownMotifs $genome -bed $maskFlag > $tmpfile`;

open IN, $tmpfile;
open OUT, ">$tmpfile2";
while (<IN>) {
	chomp;
	my @line = split /\t/;
	$line[4]=floor($line[4]);
	print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\n";
}
close IN;
close OUT;

`sort -k1,1 -k2,2n $tmpfile2 > $tmpfile`;

my $gdir = "$homeDir" . "/data/genomes/$genome/";
`homerTools extract stats "$gdir" > $tmpfile2`;

open IN, $tmpfile2;
open OUT, ">$tmpfile3";
my $count = 0;
while (<IN>) {
	$count++;
	next if ($count < 2);
	chomp;
	s/\r//g;
	my @line = split /\t/;
	next if ($line[0] eq 'genome');
	print OUT "$line[0]\t$line[1]\n";
}
close IN;
close OUT;

`bedToBigBed $tmpfile $tmpfile3 \"$name.$genome.bigBed\"`;
open OUT, ">$name.$genome.track.txt";
print OUT "track type=bigBed name=\"$name ($genome, $knownMotifs)\" description=\"$name ($genome, $knownMotifs)\" bigDataUrl=http://YOUR_SERVER/PATH_TO_FILE/$name.$genome.bigBed visibility=3\n";
close OUT;
`rm -f $tmpfile $tmpfile2 $tmpfile3`;
print STDERR "\tAll finished - need copy $name.$genome.bigBed to web server, and then change\n";
print STDERR "\tthe URL in the $name.$genome.track.txt file to reflect the location of the bigBed file\n";
