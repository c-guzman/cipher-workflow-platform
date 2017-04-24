#!/usr/bin/perl -w
#

my $hubDir = "/sanger/nfs_data/chris/hubs/";
my $set = "Data-";

if (@ARGV < 3) {
	print STDERR "<genome> [-merge <prefixName>] [options] -f <dir1> [dir2] ...\n";
	print STDERR "<genome> [-merge <prefixName>] [options] -k <keyfile.txt>\n";
	exit;
}

my $merge = "";
my $genome = $ARGV[0];
my $options = '';
my $combine = 0;
my $keyFile = '';
my @dirs = ();
for (my $i=1;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-k') {
		$keyFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-f') {
		$i++;
		for (my $j=$i;$j<@ARGV;$j++) {
			push(@dirs, $ARGV[$j]);
		}
		last;
	} elsif ($ARGV[$i] eq '-merge') {
		$merge = $ARGV[++$i];
		$set = "$merge-";
		print STDERR "\tMerging tracks into a single hub: \"$merge\"\n";
	
	} else {
		$options .= " " . $ARGV[$i];
	}
}


my $hubs = "";
if ($keyFile ne '') {
	open IN, $keyFile;
	my %sets = ();
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		if (!exists($sets{$line[0]})) {
			$sets{$line[0]} = "";
		}
		$sets{$line[0]} .= ' "' . $line[1] . '" ';
	}
	close IN;

	foreach(keys %sets) {
		my $name = $_;
		my $dirs = $sets{$name};
		$name = $set . $name;
		print STDERR "\n`makeMultiWigHub.pl $name $genome $options -d $dirs`;\n";
		`makeMultiWigHub.pl $name $genome $options -d $dirs`;
		$hubs .= " $hubDir/$name/";
	}

}

foreach(@dirs) {
	my $dir = $_;
	$name = $dir;
	$name =~ s/\.\.//g;
	$name =~ s/\///g;
	$name =~ s/ //g;
	$name = $set . $name;
	print STDERR "\n`makeMultiWigHub.pl $name $genome $options -d $dir`;\n";
	`makeMultiWigHub.pl $name $genome $options -d "$dir"`;
	$hubs .= " $hubDir/$name/";
}
	
`combineHubs.pl $hubDir/$merge $hubs`;



