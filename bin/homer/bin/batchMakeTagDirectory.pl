#!/usr/bin/env perl
use warnings;



if (@ARGV < 1) {
	print STDERR "\n\tbatchMakeTagDirectory.pl <key file> [makeTagDirectory options]...\n";
	print STDERR "\n\tKey File: Tag delimited: <directory>TAB<alignment file>\n";
	print STDERR "\t\tDuplicate directory name for multiple alignment files for same experiment\n";
	print STDERR "\t\tOption: -cpu <#> to run more than one at a time...\n";
	print STDERR "\n";
	exit;
}
my $maxCPUs = 1;
my %key = ();
open IN, $ARGV[0];
while (<IN>) {
	chomp;
	s/\r//g;
	next if (/^#/);
	my @line = split /\t/;
	if (!exists($key{$line[0]})) {
		my @a = ();
		$key{$line[0]}= \@a;
	}
	push(@{$key{$line[0]}}, $line[1]);
}
close IN; 

my $otherOptions = "";
for (my $i=1;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-cpu') {
		$maxCPUs = $ARGV[++$i];
		next;
	}
	$otherOptions .= " $ARGV[$i]";
}

my @pids = ();
my $cpus = 0;
foreach(keys %key) {
	my $dirName = $_;
	my $fileList = "";
	foreach(@{$key{$dirName}}) {
		$fileList .= " \"" . $_ . "\"";
	}
	print STDERR "`makeTagDirectory $dirName $otherOptions $fileList`\n";

	my $pid = fork();
	$cpus++;

	if ($pid == 0) {
		`makeTagDirectory "$dirName" $otherOptions $fileList `;
		exit(0);
	}

	push(@pids, $pid);
	if ($cpus >= $maxCPUs) {
		my $id = wait();
		$cpus--;
	}
}
my $id = 0;
while ($id >= 0) {
	$id= wait();
}
	

