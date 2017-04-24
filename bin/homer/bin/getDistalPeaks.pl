#!/usr/bin/env perl
use warnings;
use lib "/home/carlos/tools/homer/.//bin";
my $homeDir = "/home/carlos/tools/homer/./";


my $distance = 3000;
my $minTTSdistance = 10000;
if (@ARGV < 2) {
	printCMD();
}
my $peakFile = $ARGV[0];
my $genome = $ARGV[1];

my $proxFlag = 0;
my $intergenicFlag = 0;
my $intragenicFlag = 0;
my $annFlag = " -noann";
my $ttsFlag = "";
my $ttsColumn = 19;
my $gtfFile = "";
my $gtfFlag = '';
my $gidFlag = "";

for (my $i=2;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-prox') {
		$proxFlag = 1;
	} elsif ($ARGV[$i] eq '-gff') {
		$gtfFile = '"' . $ARGV[++$i] . '"';
		$gtfFlag = " -gff";
	} elsif ($ARGV[$i] eq '-gid') {
		$gidFlag = " -gid";
	} elsif ($ARGV[$i] eq '-gff3') {
		$gtfFile = '"' . $ARGV[++$i] . '"';
		$gtfFlag = " -gff3";
	} elsif ($ARGV[$i] eq '-gtf') {
		$gtfFile = '"' . $ARGV[++$i] . '"';
		$gtfFlag = " -gtf";
	} elsif ($ARGV[$i] eq '-d') {
		$distance = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-intergenic') {
		$intergenicFlag = 1;
		$annFlag = "";
	} elsif ($ARGV[$i] eq '-intragenic') {
		$intragenicFlag = 1;
		$annFlag = "";
	} elsif ($ARGV[$i] eq '-d') {
		$distance = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-noTTS') {
		$ttsFlag  = " -pdist -p \"$homeDir/data/genomes/$genome/$genome.tts\" ";
	} elsif ($ARGV[$i] eq '-targets') {
		$targetFlag = 1;
	} else {
		printCMD();
	}
}

print STDERR "\tDistal peaks defined as >$distance from TSS\n";



my $total = 0;
my $good = 0;
my $tmp = rand() . ".tmp";
my $tmp2 = rand() . ".2.tmp";

if ($ttsFlag ne '' && $gtfFile ne '') {
	my $x = $gtfFlag;
	if ($gtfFlag eq '-gtf') {
		$x = "";
	}
	`parseGTF.pl "$gtfFile" tts $gtfFlag $gidFlag > "$tmp2"`;
	$ttsFlag = " -pdist -p \"$tmp2\" ";
}

`annotatePeaks.pl $peakFile $genome $gtfFlag $gtfFile $gidFlag $annFlag $ttsFlag > $tmp`;
open IN, $tmp;
my $c = 0;
while (<IN>){ 
	$c++;
	chomp;
	my @line = split /\t/;
	if ($c==1) {
		if ($targetFlag) {
		} else {
			print "$line[0]";
			for (my $i=1;$i<6;$i++) {
				print "\t$line[$i]";
			}
			print "\n";
		}
		next;
	}
	$total++;
	next if ($line[9] eq 'NA');
	if ($proxFlag == 1) {
		next if (abs($line[9]) >= $distance);
	} else {
		next if (abs($line[9]) <= $distance);
	}
	if ($intergenicFlag == 1) {
		next if ($line[7] ne 'Intergenic');
	}
	if ($intragenicFlag == 1) {
		next if ($line[7] eq 'Intergenic' || $line[7] =~ /^promoter/);
	}
	if ($ttsFlag ne "") {
		if (@line > $ttsColumn) {
			next if ($line[$ttsColumn] < $minTTSdistance);
		}
	}

	if ($targetFlag) {
		print "$line[10]\n";
	} else {
		print "$line[0]";
		for (my $i=1;$i<6;$i++) {
			print "\t$line[$i]";
		}
		print "\n";
	}
	$good++;

}
close IN;
`rm "$tmp"`;
my $rate = sprintf("%.4f", $good/$total);
my $irate = 1-$rate;

$rate = $rate*100 . '%';
$irate = $irate*100 . '%';
print STDERR "\tKept $good of $total ($rate / $irate)\n";

sub printCMD {
	print STDERR "\n\tUsage: ./getDistalPeaks.txt <peakfile> <genome> [options]\n";
	print STDERR "\tOptions:\n";
	print STDERR "\t\t-d <#> (Absolute Distance from TSS to keep, default: $distance)\n";
	print STDERR "\t\t-prox (keep proximal peaks intead of distal peaks)\n";
	print STDERR "\t\t-intergenic (keep only intergenic, distal peaks)\n";
	print STDERR "\t\t-intragenic (keep only intragenic, distal peaks)\n";
	print STDERR "\t\t-noTTS (Exclude peaks within $minTTSdistance bp of the Transcription termination site)\n";
	print STDERR "\t\t-gtf <GTF file> (custom gene annotation file, -gff/-gff3 can work, but GTF is better)\n";
	print STDERR "\t\t\t-gid (parse GTF file by gene_id instead of transcript_id)\n";
	print STDERR "\t\t-targets (output target genes instead of peaks)\n";
	print STDERR "\n";
	exit;
}

