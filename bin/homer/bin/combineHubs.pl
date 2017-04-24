#!/usr/bin/perl -w


sub printCMD {
	print STDERR "\n\tusage: combineHubs.pl <new hub dir> <hub dir1> [hub dir2] ...\n";
	print STDERR "\n";
	exit;
}

if (@ARGV < 2) {
	printCMD();
}

my $newDir = $ARGV[0];
my $firstHub = $ARGV[1];
`cp -r "$firstHub" "$newDir"`;

my $genome = "";
open IN, "$firstHub/genomes.txt";
while (<IN>) {
	chomp;
	if (/genome\s+(.*)/) {
		$genome = $1;
	}
}
close IN;

my $newName = $newDir;
$newName =~ s/\/$//;
$newName =~ s/^.*\///;
print STDERR "$newName $genome\n";

open OUT, ">" . $newDir . "/hub.txt";

print OUT "hub $newName\n";
print OUT "shortLabel $newName\n";
print OUT "longLabel $newName\n";
print OUT "genomesFile genomes.txt\n";
print OUT "email cbenner\@ucsd.edu\n";

close OUT;

my $files = '';
for (my $i=1;$i<@ARGV;$i++) {
	my $dir = $ARGV[$i];
	`cp "$dir/$genome/"* "$newDir/$genome"`;
	$files .= " \"$dir/$genome/trackDb.txt\"";
}
`cat $files > "$newDir/$genome/trackDb.txt"`;


