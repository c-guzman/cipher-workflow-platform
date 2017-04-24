#!/usr/bin/perl -w
use POSIX;

my $minNumGOGenes = 10;
my $generichfile = "geneRichRegions.ann.txt";
my $genepoorfile = "geneDeserts.ann.txt";

if (@ARGV < 1) {
	print STDERR "\n\t<refGene.txt> <max distance| 100000> <gaps.ann.txt> <ID mapping file> [GO.genes files...]\n";
	print STDERR "\tCreates annotation file for each GO term with >$minNumGOGenes genes\n";
	print STDERR "\tAlso creates geneRichRegions.ann.txt and geneDeserts.ann.txt \n";
	print STDERR "\tand individual regions greater than max distance in size\n";
	print STDERR "\n";
	exit;
}

my $refGeneFile = $ARGV[0];
my $maxDistance = $ARGV[1];
my $gapFile = $ARGV[2];
my $convFile = $ARGV[3];
my @gofiles = ();
for (my $i=4;$i<@ARGV;$i++) {
	push(@gofiles, $ARGV[$i]);
}

my %genes = ();
my $refGeneFlag = 1;

open IN, $refGeneFile;
while (<IN>) {
	chomp;
	my @line= split /\t/;

	my $id = "";
	my $chr = "";
	my $dir = 0;
	my $gstart = 0;
	my $gend = 0;
	my $tss = 0;

	if ($refGeneFlag == 0) {
		$id = $line[0];
		$chr = $line[1];
		$dir = $line[2];
		$gstart = $line[3];
		$gend = $line[4];
	} else {
		$id = $line[1];
		$chr = $line[2];
		$dir = $line[3];
		$gstart = $line[4];
		$gend = $line[5];
	}

	if ($dir eq '+') {
		$tss= $gstart;
		$dir = 0;
	} else {
		$tss= $gend;
		$dir = 1;
	}
	
	if (!exists($genes{$chr})) {
		my %a = ();
		$genes{$chr} = \%a;
	}
	$genes{$chr}->{$id} = {tss=>$tss,d=>$dir,id=>$id,chr=>$chr,t=>"g"};

}
close IN;

open IN, $gapFile;
while (<IN>) {
	chomp;
	my @line =split /\t/;
	my $chr = $line[1];
	my $id = $line[0];
	my $start = $line[2];
	my $end = $line[3];
	my $dir = $line[4];
	if (!exists($genes{$chr})) {
		my %a = ();
		$genes{$chr} = \%a;
	}
	$genes{$chr}->{$id} = {tss=>$start,e=>$end,d=>$dir,id=>$id,chr=>$chr,t=>"gap"};
}
close IN;


my %ann = ();
my %deserts = ();
my %generich = ();

my $genepoorID = 1;
my $generichID = 1;

foreach(keys %genes) {
	my $chr = $_;

	my @ids = sort {$genes{$chr}->{$a}->{'tss'} <=> $genes{$chr}->{$b}->{'tss'}} keys %{$genes{$chr}};

	my $genepoorStart = 0;
	my $lastGeneStart = -1;

	my $curStart = 0;

my $NN = @ids;
#print STDERR "$chr\t$NN\n";
	for (my $i=0;$i<@ids;$i++) {

		my $printRegion = 0;
		my $gstart = $genes{$chr}->{$ids[$i]}->{'tss'}-$maxDistance;
		my $gend = $genes{$chr}->{$ids[$i]}->{'tss'}+$maxDistance;

		if ($genes{$chr}->{$ids[$i]}->{'t'} eq 'gap') {
			if ($i<@ids-1 && $genes{$chr}->{$ids[$i+1]}->{'t'} eq 'gap') {
				my $s = $genes{$chr}->{$ids[$i]}->{'e'};
				my $e = $genes{$chr}->{$ids[$i+1]}->{'tss'};
				my $poorStr = "GeneDesert$genepoorID\t$chr\t$s\t$e\t0\tGeneDesert\n";
				my $ssize = $e-$s;
				my $iid = "GeneDesert_$chr" . "_$s" . "_$e";
				$deserts{$iid} = {o=>$poorStr,s=>$ssize};
				$genepoorID++;
			}
			next;
		}
		if ($i==0) {
			if ($genes{$chr}->{$ids[$i]}->{'tss'} > $maxDistance) {
				$curStart = $genes{$chr}->{$ids[$i]}->{'tss'} - $maxDistance;
				my $s = 0;
				my $e = $curStart;
				my $poorStr = "GeneDesert$genepoorID\t$chr\t$s\t$e\t0\tGeneDesert\n";
				my $ssize = $e-$s;
				my $iid = "GeneDesert_$chr" . "_$s" . "_$e";
				$deserts{$iid} = {o=>$poorStr,s=>$ssize};
				$genepoorID++;
			} else {
				$curStart = 1;
			}
			$gstart = $curStart;
		} else {
			if ($genes{$chr}->{$ids[$i-1]}->{'t'} eq 'gap') {
				if ($genes{$chr}->{$ids[$i]}->{'tss'} - $genes{$chr}->{$ids[$i-1]}->{'e'} > $maxDistance) {
					$curStart = $genes{$chr}->{$ids[$i]}->{'tss'} - $maxDistance;
					my $s = $genes{$chr}->{$ids[$i-1]}->{'e'};
					my $e = $curStart;
					my $poorStr = "GeneDesert$genepoorID\t$chr\t$s\t$e\t0\tGeneDesert\n";
					my $ssize = $e-$s;
					my $iid = "GeneDesert_$chr" . "_$s" . "_$e";
					$deserts{$iid} = {o=>$poorStr,s=>$ssize};
					$genepoorID++;
				} else {
					$curStart = $genes{$chr}->{$ids[$i-1]}->{'e'};
				}
				$gstart = $curStart;
			} else {
				if ($genes{$chr}->{$ids[$i]}->{'tss'} - $genes{$chr}->{$ids[$i-1]}->{'tss'} > $maxDistance*2) {
					$curStart = $genes{$chr}->{$ids[$i]}->{'tss'} - $maxDistance;
					my $s = $genes{$chr}->{$ids[$i-1]}->{'tss'}+$maxDistance;
					my $e = $curStart;
					my $poorStr = "GeneDesert$genepoorID\t$chr\t$s\t$e\t0\tGeneDesert\n";
					my $ssize = $e-$s;
					my $iid = "GeneDesert_$chr" . "_$s" . "_$e";
					$deserts{$iid} = {o=>$poorStr,s=>$ssize};
					$genepoorID++;
					$gstart = $curStart;
				} else {
					$gstart = floor(($genes{$chr}->{$ids[$i]}->{'tss'} + $genes{$chr}->{$ids[$i-1]}->{'tss'})/2);
				}
			}
		}
		if ($i<@ids-1) {
			if ($genes{$chr}->{$ids[$i+1]}->{'t'} eq 'gap') {
				if ($genes{$chr}->{$ids[$i+1]}->{'tss'} - $genes{$chr}->{$ids[$i]}->{'tss'} > $maxDistance) {
					$curEnd = $genes{$chr}->{$ids[$i]}->{'tss'} + $maxDistance;
					my $s = $curEnd;
					my $e = $genes{$chr}->{$ids[$i+1]}->{'tss'};
					my $poorStr = "GeneDesert$genepoorID\t$chr\t$s\t$e\t0\tGeneDesert\n";
					my $ssize = $e-$s;
					my $iid = "GeneDesert_$chr" . "_$s" . "_$e";
					$deserts{$iid} = {o=>$poorStr,s=>$ssize};
					$genepoorID++;
					$printRegion = 1;
				} else {
					$curEnd = $genes{$chr}->{$ids[$i+1]}->{'tss'};
					$printRegion = 1;
				}
				$gend = $curEnd;
			} else {
				if ($genes{$chr}->{$ids[$i+1]}->{'tss'} - $genes{$chr}->{$ids[$i]}->{'tss'} > $maxDistance*2) {
					$curEnd = $genes{$chr}->{$ids[$i]}->{'tss'} + $maxDistance;
					$printRegion = 1;
				} else {
					$curEnd = floor(($genes{$chr}->{$ids[$i]}->{'tss'} + $genes{$chr}->{$ids[$i+1]}->{'tss'})/2);
				}
			}
			$gend = $curEnd;
		} else {
			$curEnd = $genes{$chr}->{$ids[$i]}->{'tss'} + $maxDistance;
			$printRegion = 1;
		}

		if ($printRegion == 1) {
			my $s = $curStart;
			my $e = $curEnd;
			my $richStr = "GeneRich$generichID\t$chr\t$s\t$e\t0\tGeneRich\n";
			my $ssize = $e-$s;
			my $iid = "GeneRich_$chr" . "_$s" . "_$e";
			$generich{$iid} = {o=>$richStr,s=>$ssize};
			$generichID++;
		}

		my $dir = $genes{$chr}->{$ids[$i]}->{'d'};
		my $annStr = "$ids[$i]\t$chr\t$gstart\t$gend\t$dir\tGeneTSS\n";

		$ann{$ids[$i]} = $annStr;
	}

}

open OUT, ">$generichfile";
foreach (keys %generich) {
	my $id = $_;
	my $size = $generich{$id}->{'s'};
	my $output = $generich{$id}->{'o'};
	print OUT $output;
	if ($size > $maxDistance) {
		open OUT2, ">$id.ann.txt";
		print OUT2 $output;
		close OUT2;
	}
}
close OUT;

open OUT, ">$genepoorfile";
foreach (keys %deserts) {
	my $id = $_;
	my $size = $deserts{$id}->{'s'};
	my $output = $deserts{$id}->{'o'};
	print OUT $output;
	if ($size > $maxDistance) {
		open OUT2, ">$id.ann.txt";
		print OUT2 $output;
		close OUT2;
	}
}
close OUT;

my %conv = ();
open IN, $convFile;
while (<IN>) {
	chomp;
	my @line = split /\t/;
	next if (@line < 2);
	$conv{$line[0]}=$line[1];
	$conv{$line[1]}=$line[1];
}
close IN;
my %mapping = ();
foreach(keys %ann) {
	my $id = $_;
	if (exists($conv{$id})) {
		my $gid = $conv{$id};
		$mapping{$gid} = $id;
	}
}


$ZZ=0;
for (my $i=0;$i<@gofiles;$i++) {
	my $file = $gofiles[$i];
	open IN, $file;
	while (<IN>) {
		chomp;
		my @line = split /\t/;
		my $goid = $line[0];
		my $term = $line[1];
		my @gids = split /\,/, $line[2];
		my %good = ();
		foreach(@gids) {
			if (exists($mapping{$_})) {
				$good{$mapping{$_}} = 1;
$ZZ++;
			}
		}
		my @refs = keys %good;
		if (scalar(@refs) < $minNumGOGenes) {
			next;
		}
		my $filename = $goid . "_" . $term . ".ann.txt";
		$filename =~ s/\ /\_/g;
		$filename =~ s/\:/\_/g;
		$filename =~ s/\(/\_/g;
		$filename =~ s/\)/\_/g;
		$filename =~ s/\,//g;
		$filename =~ s/\'//g;
		$filename =~ s/\//\_/g;
		open OUT, ">$filename";
		foreach(@refs) {
			print OUT $ann{$_};
		}
		close OUT;
	}
	close IN;
}
print STDERR "\tTotal GO genes used in genome ontology groups: $ZZ\n";
