#!/usr/bin/env perl
use warnings;
use lib "/home/carlos/tools/homer/.//bin";
my $homeDir = "/home/carlos/tools/homer/./";
# 
# This file is part of HOMER
#
# HOMER is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HOMER is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

use HomerConfig;


my $baseUrl = "http://hgdownload.cse.ucsc.edu/goldenPath/";


sub printCMD {
	print STDERR "\n\tUsage: updateUCSCGenomeAnnotations.pl <ucsc.txt file> [options]\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-tssStart <#> (default: -1000)\n";
	print STDERR "\t\t-tssEnd <#> (default: 100)\n";
	print STDERR "\t\t-ttsStart <#> (default: -100)\n";
	print STDERR "\t\t-ttsEnd <#> (default: 1000)\n";
	print STDERR "\t\t-baseURL <http...> (default: $baseUrl)\n";
	print STDERR "\n";
	exit;
}

if (@ARGV < 1) {
	printCMD();
}
my $ucscFile = $ARGV[0];
my $promoterStart = -1000;
my $promoterEnd = 100;
my $ttsStart = -100;
my $ttsEnd = 1000;
for (my $i=1;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-tssStart') {
		$promoterStart = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-tssEnd') {
		$promoterEnd = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-ttsStart') {
		$ttsStart = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-ttsEnd') {
		$ttsEnd = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-baseURL') {
		$baseUrl = $ARGV[++$i];
	} else {
		printCMD();
	}
}

my $config = HomerConfig::loadConfigFile();


print STDERR "\tGenome\tOrganism\tVersion\tExisting Version\n";
my @genomes = ();
open IN, $ucscFile;
while (<IN>) {
	chomp;
	s/\r//g;
	next if (/^#/);
	my @line = split /\t/;
	my $genome = $line[0];
	my $org = $line[1];
	my $sFlag = 0;
	my $rFlag = 1;
	my $version = 'custom';
	if (@line > 2) {
		$version = $line[2];
	}
	my $g = {genome=>$genome,org=>$org,single=>$sFlag,repeats=>$rFlag,version=>$version};

	my $existingVersion = "Not Installed";
	if (exists($config->{'GENOMES'}->{$genome})) {
		$existingVersion = $config->{'GENOMES'}->{$genome}->{'version'};
	}
	print STDERR "\t$genome\t$org\t$version\t$existingVersion\n";
	push(@genomes, $g);
}
close IN;

print STDERR "\tWating 10 seconds in case you want to review the changes (hit ctrl+C to cancel)\n\t\t";
for (my $i=10;$i>0;$i--) {
	print STDERR " $i";
	`sleep 1`;
}
print STDERR "\n";
	
my $cpgIslandFilename = "cpgIslandExt.txt";
my $cpgOutput = "cpgIsland.ann.txt";

foreach(@genomes) {
	my $gInfo = $_;
	my $genome = $gInfo->{'genome'};
	print STDERR "\tUpdating $genome...\n";
	my $orgName= $gInfo->{'org'};
	my $repeatsFlag = 1;
	my $singleFlag = 0;
	my $orgUrl = $baseUrl . $genome . "/database/";

	my $refGeneFilename = "refGene.txt"; #mir RNA definitions in the refGene file
	if ($orgName eq 'yeast') {
		$refGeneFilename = "sgdGene.txt";
	}

	#create needed output directories
	my $outputdir = "../data/genomes/$genome/";
	`mkdir -p "$outputdir"`;
	$outputdir = "../data/genomes/$genome/preparsed/";
	`mkdir -p "$outputdir"`;
	`chmod 775 "$outputdir"`;
	$outputdir = "../data/genomes/$genome/annotations/";
	`rm -r "$outputdir"`;
	`mkdir -p "$outputdir"`;
	$outputdir = "../data/genomes/$genome/annotations/repeats/";
	`mkdir -p "$outputdir"`;
	$outputdir = "../data/genomes/$genome/annotations/basic/";
	`mkdir -p "$outputdir"`;
	$outputdir = "../data/genomes/$genome/annotations/custom/";
	`mkdir -p "$outputdir"`;


	my $genomeURL = $baseUrl . $genome . "/bigZips/";
	my $files = getFTPfiles($genomeURL);

	my $genomeFileName = '';
	if (exists($files->{"chromFa.tar.gz"})) {
		$singleFlag = 0;
		$genomeFileName = "chromFa.tar.gz";
	} elsif (exists($files->{"chromFa.zip"})) {
		$singleFlag = 0;
		$genomeFileName = "chromFa.zip";
	} elsif (exists($files->{"ScaffoldFa.zip"})) {
		$singleFlag = 0;
		$genomeFileName = "ScaffoldFa.zip";
	} elsif (exists($files->{"GroupFa.zip"})) {
		$singleFlag = 0;
		$genomeFileName = "GroupFa.zip";
	} elsif (exists($files->{"$genome.fa.gz"})) {
		$singleFlag = 1;
		$genomeFileName = "$genome.fa.gz";
	} else {
		print STDERR "Error! - Could not find one of the expected genome FASTA files (chromFa.tar.gz, chromFa.zip, $genome.fa.gz)\n";
		print STDERR "Available files:\n";
		foreach(keys %{$files}) {
			print STDERR "\t\t$_\n";
		}
		print STDERR "Error! - Could not find one of the expected genome FASTA files (chromFa.tar.gz, chromFa.zip, $genome.fa.gz)\n";
		exit;
	}

	$files = getFTPfiles($orgUrl);
	$repeatsFlag = 1;
	if (exists($files->{"rmsk.txt.gz"}) || exists($files->{"gap.txt.gz"})) {
		$repeatsFlag = 2;
	}
	print "\t$genome\tsingleFlag=$singleFlag\t$genomeFileName\trepeatsFlag=$repeatsFlag\n";

	#download genome FASTA files
	my $url = "";
	my @chr = ();

	if ($singleFlag) {
		$url = $baseUrl . $genome . "/bigZips/" . $genomeFileName;
		`wget -O genome.fa.gz "$url"`;
		`gunzip genome.fa.gz`;
		`mv genome.fa "../data/genomes/$genome/"`;
		@chr = ("");
	} else {
		$url = $baseUrl . $genome . "/bigZips/" . $genomeFileName;
		if ($genomeFileName =~ /tar\.gz/) {
			`wget -O genome.tar.gz "$url"`;
			`gunzip -f genome.tar.gz`;
			`tar xvf genome.tar > "$genome.chrListOut.txt" 2> "$genome.chrListErr.txt"`;
		} elsif ($genomeFileName =~ /\.zip/) {
			`wget -O genome.zip "$url"`;
			`unzip -o genome.zip > "$genome.chrListOut.txt" 2> "$genome.chrListErr.txt"`;
		}
		`cat "$genome.chrListOut.txt" "$genome.chrListErr.txt" > "$genome.chrList.txt"`;
			
		open IN, "$genome.chrList.txt";
		while (<IN>) {
			chomp;
			my $file = $_;
			next if ($file =~ /^Archive/);
			$file =~ s/\s*Inflating\:\s*//;
			$file =~ s/\s*inflating\:\s*//;
			$file =~ s/^\s*x\s+//;
			$file =~ s/^\s*X\s//;
			$file =~ s/^\s*//;
			$file =~ s/\s*$//;
			my $chrFileName = $file;
			next unless ($file =~ /\.fa/);
			if ($chrFileName =~ /\//) {
				my $chrDir = $chrFileName;
				$chrDir =~ s/\/.*?$/\//;
				$chrFileName =~ s/^.*\///;
				`mv "$file" "$chrFileName"`;
				`rmdir "$chrDir"`;
			}
			my $chrName = $chrFileName;
			$chrName =~ s/\.fa//i;
			push(@chr, $chrName);
			`mv "$chrFileName" "../data/genomes/$genome/"`;
		}
		#print STDERR "\tChromosomes: @chr\n";
		`rm -f genome.tar genome.zip "$genome.chrList.txt" "$genome.chrListOut.txt" "$genome.chrListErr.txt"`;
	}

	#old: my @chr = getChr($genome);

	print STDERR "\tProcessing genome $genome\n";
	if ($singleFlag) {
	} 
	foreach(@chr) {
		print STDERR "\t\t$genome\t$_\n";
	}

	my $repeatFileStr = ();
	my @gapFiles = ();
	open OUT, ">gaps.ann.txt";
	open OUT2, ">centromeres.ann.txt";
	my $gapID = 1;
	if ($repeatsFlag == 0 || $repeatsFlag == 1) {
		foreach(@chr) {
			my $chr = $_;
			my $repeatFile = $chr . "_rmsk.txt";
			if ($singleFlag) {
				$repeatFile = 'rmsk.txt';
			}
			my $url = $orgUrl . $repeatFile . ".gz";
			if (!exists($files->{"$repeatFile.gz"})) {
				print STDERR "\tWarning: Could not find repeats for $chr ($genome)\n";
			} else {
				`wget -O $repeatFile.gz $url`;
				`gunzip -f $repeatFile.gz`;
				$repeatFileStr .= " " . $repeatFile;
			}
	
			my $gapFile = $chr . "_gap.txt";
			if ($singleFlag) {
				$gapFile = 'gap.txt';
			}
			$url = $orgUrl . $gapFile . ".gz";
			`wget -O $gapFile.gz $url`;
			`gunzip -f $gapFile.gz`;
			open IN, $gapFile;
			while (<IN>) {
				chomp;
				my @line = split /\t/;
				my $str = "$line[7]-$gapID\t$line[1]\t$line[2]\t$line[3]\t0\t$line[7]\n";
				print OUT $str;
				if ($line[7] eq 'centromere') {
					print OUT2 $str;
				}
				$gapID++;
			}
			close IN;
			`rm $gapFile`;
		}
		close OUT;
		close OUT2;
	} else {
		my $repeatFile = "rmsk.txt";
		my $url = $orgUrl . $repeatFile . ".gz";
		`wget -O $repeatFile.gz $url`;
		`gunzip -f $repeatFile.gz`;
		$repeatFileStr .= " " . $repeatFile;
		
		my $gapFile = "gap.txt";
		$url = $orgUrl . $gapFile . ".gz";
		`wget -O $gapFile.gz $url`;
		`gunzip -f $gapFile.gz`;
		open IN, $gapFile;
		while (<IN>) {
			chomp;
			my @line = split /\t/;
			my $str = "$line[7]-$gapID\t$line[1]\t$line[2]\t$line[3]\t0\t$line[7]\n";
			print OUT $str;
			if ($line[7] eq 'centromere') {
				print OUT2 $str;
			}
			$gapID++;
		}
		close IN;
		`rm $gapFile`;
	}

	$outputdir = "../data/genomes/$genome/";


	makeChromSizeFile($outputdir, $outputdir . "chrom.sizes");

	$url = $orgUrl . $refGeneFilename . ".gz";
	`wget -O $refGeneFilename.gz $url`;
	`gunzip -f $refGeneFilename.gz`;
	`./UCSC/makeGenomeAnnotationFromRefGenes.pl $refGeneFilename -prefix $genome -org $orgName -tssStart $promoterStart -tssEnd $promoterEnd -ttsStart $ttsStart -ttsEnd $ttsEnd`;
	`mv *.ann.txt "$outputdir"/annotations/basic/`;

	$url = $orgUrl . $cpgIslandFilename . ".gz";
	`wget -O $cpgIslandFilename.gz $url`;
	`gunzip -f $cpgIslandFilename.gz`;
	processCpG($cpgIslandFilename, $cpgOutput);
	`rm $cpgIslandFilename`;

	#my $goFiles = "../data/GO/biological_process.genes";
	#my $dsize = $desertSize{$orgName};
	#my $accessionFile = "../data/accession/" . $orgName . "2gene.tsv";
	#`./createGOtoGenomeOntology.pl $refGeneFilename $dsize gaps.ann.txt $accessionFile $goFiles`;

	#repeats
	`./UCSC/makeRepeatGroups.pl ./ ../data/genomes/$genome/annotations/repeats/`;

	`cat $genome.raw.annotation $genome.intron.annotation > $genome.basic.raw.annotation `;
	`cat $genome.raw.annotation cpgIsland.ann.txt repeats.ann.txt $genome.intron.annotation > $genome.full.raw.annotation `;

	`assignGenomeAnnotation $genome.basic.raw.annotation $genome.basic.raw.annotation -prioritize .tmp`;
	`mv .tmp $genome.basic.annotation`;

	`assignGenomeAnnotation $genome.basic.annotation $genome.full.raw.annotation -prioritize .tmp`;
	`mv .tmp $genome.full.annotation`;

	`UCSC/getIntergenicAnn.pl $genome.basic.annotation ../data/genomes/$genome/chrom.sizes > intergenic.ann.txt`;
	`mv cpgIsland.ann.txt ../data/genomes/$genome/annotations/basic/`; 
	`mv intergenic.ann.txt ../data/genomes/$genome/annotations/basic/`; 
	
	`rm $genome.raw.annotation`;
	`rm $genome.intron.annotation`;
	`rm $genome.basic.raw.annotation`;
	`rm $genome.full.raw.annotation`;

	`mv $genome.basic.annotation "$outputdir"`; 
	`mv $genome.full.annotation "$outputdir"`; 
	`mv $genome.tss "$outputdir"`; 
	`mv $genome.tts "$outputdir"`; 
	`mv $genome.splice5p "$outputdir"`; 
	`mv $genome.splice3p "$outputdir"`; 
	`mv $genome.aug "$outputdir"`; 
	`mv $genome.stop "$outputdir"`; 
	`mv $genome.rna "$outputdir"`; 
	`mv $genome.miRNA "$outputdir"`; 
	`mv repeats.rna "$outputdir"/$genome.repeats`; 



	`rm $refGeneFilename`;
	`rm *rmsk.txt`;
	`rm *.ann.txt`;

	#update config settings:
	if (!exists($config->{'GENOMES'}->{$genome})) {
		print STDERR "\tMaking new entry for $genome in the config.txt file\n";
	}
	my @params = ($orgName,"default");
	my $g = {org=>$gInfo->{'org'},version=>$gInfo->{'version'},location=>"data/genomes/$genome/",
						url=>"LocalUpdate",desc=>"$orgName genome and annotation for UCSC $genome",
						parameters=>\@params};
	$config->{'GENOMES'}->{$genome} = $g;

	HomerConfig::printConfigFile($config);

}

sub processMIRNA {
	my ($inFile, $outFile) = @_;
	open OUT, ">$outFile";
	open IN, $inFile;
	while (<IN>) {
		chomp;
		my @line = split /\t/;
		$line[2]++;
		print OUT "$line[4]\t$line[1]\t$line[2]\t$line[3]\t$line[6]\tmiRNA\n";
	}
	close OUT;
	close IN;
}
sub processCpG {
	my ($file,$output,$col) = @_;
	open IN, $file;
	my $id = 1;
	my $prefix = "CpG";
	open OUT, ">$output";
	while (<IN>) {
		chomp;
		my @line = split /\t/;
		if ($line[0] =~ /^\d+$/) {
			$col = 1;
		} else {
			$col = 0;
		}
		#print STDERR "$_\n";
		print OUT $prefix . "-" . $id++ . "\t$line[$col]\t$line[$col+1]\t$line[$col+2]\t0\tCpG-Island\n";
	}
	close IN;
	close OUT;
}

sub getChr {
	my ($genome) = @_;
	`ls ../data/genomes/$genome/*.fa* > .ls`;
	open IN, ".ls";
	my @chr = ();
	while (<IN>) {
		chomp;
		s/^.*chr/chr/;
		s/\.fa.*$//;
		push(@chr, $_);
	}
	close IN;
	`rm .ls`;
	return @chr;
}

sub getFTPfiles {
	my ($url) = @_;

	`wget -O .ls $url`;
	my %files = ();
	open LS, ".ls";
	while (<LS>) {
		chomp;
		if (/a href=\"(.*?)\"/i) {
			my $file = $1;
			$files{$file} = 1;
		}
	}
	close LS;
	`rm .ls`;

	return \%files;
}
sub makeChromSizeFile {
	my ($dir, $chromSizeFile) = @_;
	`homerTools extract stats "$dir" > .tmp`;
	open IN, ".tmp";
	open OUT, ">$chromSizeFile";
	while (<IN>) {
		chomp;
		my @line = split /\t/;
		next if ($line[0] eq 'ChrName');
		next if ($line[0] eq 'genome');
		print OUT "$line[0]\t$line[1]\n";
	}
	close IN;
	close OUT;
	`rm ".tmp"`;
}

