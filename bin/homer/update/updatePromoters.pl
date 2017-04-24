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
my $config = HomerConfig::loadConfigFile();

sub printCMD {
	print STDERR "\n\tUsage: updatePromoters.pl <promoters.txt file> [options]\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-rna (Generate -RNA, -3UTR, and -5UTR sets from UCSC refseq mRNA)\n";
	print STDERR "\t\t-dist <#> (max distance between promoters to be considered redundant, def:500)\n";
	print STDERR "\t\t-size <#> (max size of promoters for analysis with findMotifs.pl,def:4000)\n";
	print STDERR "\n";
	exit;
}


if (@ARGV < 1) {
	printCMD();
}

my $baseURL = "http://hgdownload.cse.ucsc.edu/goldenPath/";
my $redunDist = 500;
my $size = 4000;
my $rnaFlag = 0;

my $promoterFile = $ARGV[0];
for (my $i=1;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-rna') {
		$rnaFlag = 1;
	} elsif ($ARGV[$i] eq '-size') {
		$size = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-dist') {
		$redunDist = $ARGV[++$i];
	} else {
		printCMD();
	}
}

print STDERR "\tPromoterSet\tGenome\tOrganism\tIDtype\tUCSC\tVerison\tInstalled Version\n";
my @sets = ();
open IN, $promoterFile;
while (<IN>) {
	chomp;
	s/\r//g;
	next if (/^#/);
	my @line = split /\t/;
	my $name = $line[0];
	my $genome = $line[1];
	my $organism = $line[2];
	my $idtype = $line[3];
	my $UCSC = $line[4];
	my $version = $line[5];

	my $g = {name=>$name,genome=>$genome,org=>$organism,version=>$version,source=>$UCSC,idtype=>$idtype};
	my $oldVersion = 'Not Installed';
	if (exists($config->{'PROMOTERS'}->{$name})) {
		$oldVersion = $config->{'PROMOTERS'}->{$name}->{'version'};
	}
	my $warning='';
	if (!exists($config->{'GENOMES'}->{$genome})) {
		$warning = "!! Could NOT find genome \"$genome\"!  Need to install that first!\n";
	}
	print STDERR "\t$name\t$genome\t$organism\t$idtype\t$UCSC\t$version\t$oldVersion\n";
	print STDERR $warning;

	if ($warning eq '') {
		push(@sets, $g);
	}
}
close IN;

print STDERR "\n";
if ($rnaFlag) {
	print STDERR "\tRNA set will also be generated from UCSC genomes\n\n";
}

print STDERR "\tWating 10 seconds in case you want to review the changes (hit ctrl+C to cancel)\n\t\t";
for (my $i=10;$i>0;$i--) {
	print STDERR " $i";
	`sleep 1`;
}
print STDERR "\n";

foreach(@sets) {
	my $g = $_;
	my $name = $g->{'name'};
	my $org = $g->{'org'};
	my $genome = $g->{'genome'};
	my $version = $g->{'version'};
	my $idtype = $g->{'idtype'};
	my $mgenome = $genome . "r";

	`loadPromoters.pl -name $org -org $org -id $idtype -force -version $version -genome $genome -tss ../data/genomes/$genome/$genome.tss -size $size -dist $redunDist`;

	next if ($rnaFlag == 0 || $g->{'source'} ne 'UCSC');

	my $url = $baseURL . "$genome" . "/bigZips/refMrna.fa.gz";	
	if ($org eq 'yeast') {
		$url = $baseURL . "$genome" . "/bigZips/mrna.fa.gz";	
	}
	`wget -O refMrna.fa.gz $url`;
	`gunzip -f refMrna.fa.gz`;

	$url = $baseURL . "$genome" . "/database/refGene.txt.gz";
	if ($org eq 'yeast') {
		$url = $baseURL . "$genome" . "/database/sgdGene.txt.gz";	
	}
	`wget -O refGene.txt.gz $url`;
	`gunzip -f refGene.txt.gz`;

	`./promoters/fasta2tabMrna.pl refMrna.fa > $org-mRNA.seq`;
	`cut -f1 $org-mRNA.seq > $org-mRNA.base`;
	`addGeneAnnotation.pl $org-mRNA.base $org no > tmp.ann`;
	`./promoters/getRedunMrna.pl tmp.ann > $org-mRNA.redun`;
	`cut -f1 tmp.ann | sort | uniq > $org-mRNA.base.gene`;

	#UTRs
	`./promoters/parseUTRs.pl $org $org-mRNA.seq refGene.txt`;

	`cp $org-mRNA-3UTR.seq ../data/promoters/$org-mRNA-3UTR.seq`;
	`mv $org-mRNA-3UTR.seq ../data/promoters/$org-mRNA-3UTR.mask`;
	`cp $org-mRNA-5UTR.seq ../data/promoters/$org-mRNA-5UTR.seq`;
	`mv $org-mRNA-5UTR.seq ../data/promoters/$org-mRNA-5UTR.mask`;
	`cp $org-mRNA.seq ../data/promoters/$org-mRNA.seq`;
	`mv $org-mRNA.seq ../data/promoters/$org-mRNA.mask`;

	`cp $org-mRNA.base ../data/promoters/$org-mRNA-3UTR.base`;
	`cp $org-mRNA.base ../data/promoters/$org-mRNA-5UTR.base`;
	`mv $org-mRNA.base ../data/promoters/$org-mRNA.base`;

	`cp $org-mRNA.base.gene ../data/promoters/$org-mRNA-3UTR.base.gene`;
	`cp $org-mRNA.base.gene ../data/promoters/$org-mRNA-5UTR.base.gene`;
	`mv $org-mRNA.base.gene ../data/promoters/$org-mRNA.base.gene`;

	`cp $org-mRNA.redun ../data/promoters/$org-mRNA-3UTR.redun`;
	`cp $org-mRNA.redun ../data/promoters/$org-mRNA-5UTR.redun`;
	`mv $org-mRNA.redun ../data/promoters/$org-mRNA.redun`;


	my @params = ($org,$genome,$idtype,0,0);
	my $gg = {org=>$org,version=>$version,url=>"LocalUpdate",location=>"data/promoters/",
					desc=>"$org refseq RNA",parameters=>\@params};
	$config->{'PROMOTERS'}->{"$org-mRNA"} = $gg;

	$gg = {org=>$org,version=>$version,url=>"LocalUpdate",location=>"data/promoters/",
					desc=>"$org refseq RNA 3'UTR",parameters=>\@params};
	$config->{'PROMOTERS'}->{"$org-mRNA-3UTR"} = $gg;

	$gg = {org=>$org,version=>$version,url=>"LocalUpdate",location=>"data/promoters/",
					desc=>"$org refseq RNA 5'UTR",parameters=>\@params};
	$config->{'PROMOTERS'}->{"$org-mRNA-5UTR"} = $gg;

	HomerConfig::printConfigFile($config);

	`rm tmp.ann refGene.txt refMrna.fa`;

}
