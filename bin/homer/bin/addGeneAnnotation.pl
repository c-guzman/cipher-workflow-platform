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

use HomerConfig;

if (@ARGV < 2) {
	print STDERR "\n\tusage: addGeneAnnotation.pl <data file> <organism> [header: yes|no] [keep unknownIDs: yes|no]\n";
	print STDERR "\n\tOutputs new data file with annotation informtion at the end of each row.\n\n";
	exit;
}


my $config = HomerConfig::loadConfigFile();

my $file = $ARGV[0];
my $organism = $ARGV[1];
my $header = "no";
my $keep = "no";
if (@ARGV > 2) {
	if ($ARGV[2] eq 'yes') {
		$header = "yes";
	}
}
if (@ARGV > 3) {
	if ($ARGV[3] eq 'yes') {
		$keep = "yes";
	}
}

if (!exists($config->{'ORGANISMS'}->{$organism})) {
	print STDERR "!!! Warning - can't seem to find your organism in the config file !!!\n";	
}

my $rand = rand();
my $tmpfile = $rand . ".tmp";
my $tmpfile2 = $rand . ".2.tmp";

`convertIDs.pl "$file" $organism gene $header yes $keep > $tmpfile`;
if ($header eq 'yes') {
	`addDataHeader.pl $tmpfile "$homeDir/data/accession/$organism.description" > $tmpfile2`;
} else {
	`addData.pl $tmpfile "$homeDir/data/accession/$organism.description" > $tmpfile2`;
}

open IN, $tmpfile2;
while (<IN>) {
	chomp;
	s/\r//g;
	print "$_\n";
}
close IN;

`rm $tmpfile $tmpfile2`;
