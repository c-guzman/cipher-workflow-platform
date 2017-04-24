#!/usr/bin/env perl
use warnings;
use lib "/home/carlos/tools/homer/.//bin";
my $homeDir = "/home/carlos/tools/homer/./";


# Copyright 2009 - 2014 Christopher Benner <cbenner@salk.edu>
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


use POSIX;
use HomerConfig;

my $defaultdirectory = $homeDir . "/data/genomes/";
my $config = HomerConfig::loadConfigFile();

sub printCMD {

	print STDERR "\n\tProgram will prepare a custom genome for use with HOMER\n\n";
	print STDERR "\tUsage: loadGenome.pl <Required Parameters ...> [options]\n";
	print STDERR "\n\tNOTE: If your genome is available at UCSC, consider using the update scripts\n";
	print STDERR "\t\tlocated in the $homeDir/update directory\n";
	print STDERR "\n\tRequired Parameters:\n";
	print STDERR "\t\t-name <genome name> (i.e. hg19, tair10, etc.)\n";
	print STDERR "\t\t-fasta <genome fasta file> (Single genome sequence, preferrabley soft masked, unzipped)\n";
	print STDERR "\t\t-gtf <gene annotation file> (Transcript annotation in gtf format, -gff/-gff3 to use them)\n";
	print STDERR "\t\t\t(Always best to use a gtf file whenever possible, gffs can be organized differently...)\n";
	print STDERR "\t\t-org <organism name, ok to use 'null'>\n";
	print STDERR "\n\tOther options:\n";
	print STDERR "\t\t-force (Overwrite any existing genome with the same name)\n";
	print STDERR "\t\t-promoters <promoter set name> (Create promoter set with genome and gtf files)\n";
	print STDERR "\t\t\t-id <idtype> (options: gene, refseq, unigene, ensembl, or custom, default: custom)\n";
	print STDERR "\t\t-version <version id> (Assign version, versions starting with 'v' are managed\n";
	print STDERR "\t\t\tby the configureHomer.pl script and my be overwritten, default: custom)\n";
	print STDERR "\t\t-tid (Use transcript_id instead of gene_id to identify the transcripts from GTF files)\n";
	#print STDERR "\t\t-repeats <repeat peak file> (tab-delimited file of repeats, column name & example:\n";
	#print STDERR "\t\t\t\t1:id 2:chr 3:start 4:end 5:strand 6:name 7:class/family\n";
	#print STDERR "\t\t\t\tR1375 Scaffold000426 29577 29851 + MITE-Tourist-like MITE/rice\n";
	print STDERR "\n";
	exit;
}

my $genome = "";
my $fastaFile = "";
my $gtfFile = '';
my $gtfOption = "";
my $org = '';
my $repeatsFile = '';
my $version = "custom";
my $forceFlag = 0;
my $promoterSet = "";
my $idtype = "custom";

if (@ARGV < 1) {
	printCMD();
}

for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-name') {
		$genome = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-fasta' || $ARGV[$i] eq '-fa') {
		$fastaFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-gid') {
		$gtfOption .= " -gid";
	} elsif ($ARGV[$i] eq '-id') {
		$idtype = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-tid') {
		$gtfOption .= " -tid";
	} elsif ($ARGV[$i] eq '-gtf') {
		$gtfFile = $ARGV[++$i];
		$gtfOption .= " ";
	} elsif ($ARGV[$i] eq '-gff') {
		$gtfFile = $ARGV[++$i];
		$gtfOption .= " -gff";
	} elsif ($ARGV[$i] eq '-gff3') {
		$gtfFile = $ARGV[++$i];
		$gtfOption .= " -gff3";
	} elsif ($ARGV[$i] eq '-org') {
		$org = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-version') {
		$version = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-repeats') {
		$repeatsFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-promoters') {
		$promoterSet = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-force') {
		$forceFlag = 1;
	} 
}

my $genomeDir = $defaultdirectory . $genome . "/";


print STDERR "\tCurrent Settings:\n";
print STDERR "\t\tGenome name = $genome\n";
print STDERR "\t\tFASTA file = $fastaFile\n";
print STDERR "\t\tGTF file = $gtfFile $gtfOption\n";
print STDERR "\t\tOrganism = $org\n";
print STDERR "\t\tVersion = $version\n";
if ($repeatsFile ne '') {
	print STDERR "\t\tRepeats = $repeatsFile\n";
}
print STDERR "\n\tGenome will be stored in directory:\n";
print STDERR "\t\t$genomeDir\n";
if ($promoterSet ne '') {
	print STDERR "\tWill create promoter set for genome (promoter set name: $promoterSet)\n";
}

if ($genome eq '' || $gtfFile eq '' || $fastaFile eq '' || $org eq '') {
	print STDERR "!! Options -name, -fasta, -gtf, and -org are REQUIRED !!\n";
	printCMD();
	exit;
}

if (exists($config->{'GENOMES'}->{$genome})) {
	print STDERR "!! Warning: Genome \"$genome\" already exists in config.txt !!\n";
	if ($forceFlag == 0) {
		print STDERR "!! Add the option \"-force\" to the command to overwrite !!\n";
		exit;
	} else {
		print STDERR "\tOverwriting...\n";
	}
}

if (exists($config->{'PROMOTERS'}->{$promoterSet})) {
	print STDERR "!! Warning: Promoter Set \"$promoterSet\" already exists in config.txt !!\n";
	if ($forceFlag == 0) {
		print STDERR "!! Add the option \"-force\" to the command to overwrite !!\n";
		exit;
	} else {
		print STDERR "\tOverwriting...\n";
	}
}

print STDERR "\n\tWating 10 seconds in case you want to review the changes (hit ctrl+C to cancel)\n\t\t";
for (my $i=10;$i>0;$i--) {
	print STDERR " $i";
	`sleep 1`;
}
print STDERR "\n";


# create genome directory, copy FASTA, and update config file
`mkdir -p "$genomeDir"`;
`cp "$fastaFile" "$genomeDir/genome.fa"`;
`mkdir -p "$genomeDir/preparsed"`;
`chmod 775 "$genomeDir/preparsed"`;
`mkdir -p "$genomeDir/annotations/"`;
`mkdir -p "$genomeDir/annotations/basic/"`;
`mkdir -p "$genomeDir/annotations/repeats/"`;
`mkdir -p "$genomeDir/annotations/custom/"`;

my @params = ($org,"default");
my $g = {org=>$org,version=>$version,location=>"data/genomes/$genome/",
			url=>"LocalUpdate",desc=>"$org genome and annotation ($genome)",
			parameters=>\@params};
$config->{'GENOMES'}->{$genome} = $g;
HomerConfig::printConfigFile($config);

# parse the GTF file a bunch of times...
`parseGTF.pl "$gtfFile" rna $gtfOption > "$genomeDir/$genome.rna"`;
`parseGTF.pl "$gtfFile" tss $gtfOption > "$genomeDir/$genome.tss"`;
`parseGTF.pl "$gtfFile" tts $gtfOption > "$genomeDir/$genome.tts"`;
`parseGTF.pl "$gtfFile" ann $gtfOption > "$genomeDir/$genome.ann"`;

`assignGenomeAnnotation "$genomeDir/$genome.ann" "$genomeDir/$genome.ann" -prioritize "$genomeDir/$genome.basic.annotation"`;
if ($repeatsFile ne '') {
	open OUT, ">>$genomeDir/$genome.ann";
	open IN, $repeatsFile;
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		my @r = ();
		next if (@line < 6);
		my $name = "$line[5]|";
		my $family = $line[5];
		if (@line > 6) {
			@r = split /\//,$line[6];
		}
		if (@r < 1) {
			$name .= "$line[5]|$line[5]";
		} elsif (@r < 2) {
			$name .= "$r[0]|$r[0]";
			$family = $r[0];
		} else {
			$name .= "$r[0]|$r[1]";
			$family = $r[0];
		}
			
		print OUT "$name\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$family\t0\n";
	}
	close IN;
	close OUT;
	`assignGenomeAnnotation "$genomeDir/$genome.ann" "$genomeDir/$genome.ann" -prioritize "$genomeDir/$genome.full.annotation"`;
}

createAnnotationFiles("$genomeDir/$genome.ann");
`rm "$genomeDir/$genome.ann"`;


if ($promoterSet ne '') {
	`loadPromoters.pl -name $promoterSet -org $org -id $idtype -genome $genome -tss "$genomeDir/$genome.tss" -version $version -force`

}

sub createAnnotationFiles {
	my ($file) = @_;
	my %ann = ();	
	my %map = ();
	$map{'P'} = "promoters";
	$map{'E'} = "exons";
	$map{'I'} = "introns";
	$map{'3UTR'} = "utr3";
	$map{'5UTR'} = "utr5";
	$map{'TTS'} = "tts";
	my %rflags = ();
	
	open IN, $file;
	while (<IN>) {
		chomp;
		s/\r//g;
		my $og = $_;
		my @line = split /\t/;
		next if (@line < 6);
		my $name = $line[0];
		$name =~ s/\-HOMER\d+$//;
		my @r = split /\|/,$name;
		my $flag = 0;
		my $ann = $line[5];
		if (exists($map{$ann})) {
			$ann = $map{$ann};
		} elsif (@r > 2) {
			$ann = $name;
			$flag = 1;
		} else {
		}
		
		if (!exists($ann{$ann})) {
			my @a = ();
			$ann{$ann} = \@a;
			$rflags{$ann} = $flag;
		}
		push(@{$ann{$ann}},$og);
	}
	close IN;
	
	#now print out all annotations
	foreach(keys %ann) {
		my $ann = $_;
		if ($rflags{$ann}) {	
			open OUT, ">$genomeDir/annotations/repeats/$ann.ann.txt";
		} else {
			open OUT, ">$genomeDir/annotations/basic/$ann.ann.txt";
		}
		foreach(@{$ann{$ann}}) {
			print OUT "$_\n";
		}
		close OUT;
	}
}
