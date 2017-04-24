#!/usr/bin/env perl
use warnings;



if (@ARGV < 1) {
	print STDERR "\n\tusage: changeNewLine.pl <text file>\n";
	print STDERR "\ti.e. changeNewLine.pl data.txt\n";
	print STDERR "\n\tThis program changes text files such that they use newline characters of the form \"\\n\"\n";
	print STDERR "\t(i.e. creates a UNIX style text file)\n"; 
	print STDERR "\t!! The new file will overwrite the old file !!\n";
	print STDERR "\n";
	exit;
}
my $rand = rand() . ".tmp";
`perl -pe 'if (s/\r\n/\n/g) {} else { s/\r/\n/g; }' "$ARGV[0]" > "$rand"`;
`mv "$rand" "$ARGV[0]"`;
