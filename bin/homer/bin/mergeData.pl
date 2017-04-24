#!/usr/bin/env perl
use warnings;


# example program that adds data from one file to the other file, printing out a new file

#first check if there are any command line arguments - perl stores these in a special variable
# named @ARGV
if (@ARGV < 2) { # check if there are less than 2 arguments
	#if so tell the user they are stupid and must give 2 file names
	print STDERR "\n\tUsage: mergeData.pl <file1> <file2> [# lines of HEADER] [-accVer]\n";

	print STDERR "\n\tMust supply to file names on the command line: i.e. 'perl mergeData.pl file1 file2'\n";
	#the print statement will send the string(text) to a file stream.
	#In this case the file stream is called "STDERR", which prints directly to the screen
	print STDERR "\tThis program adds the contents the the 2nd file to that of the 1st\n";
	print STDERR "\tAppend -accVer to ignore \".#\" at the end of identifiers\n";
	#descriptions are always nice.
	print STDERR "\n";
	exit;  #exit the program
}

# - if we get here we know we have 2 command line arguments.
# lets read in the 2nd file to start - that way we can add it to the 1st file while we read the 1st
# file to make things easier.
my $header = 0;
if (@ARGV > 2) {
	$header = $ARGV[2];
}
my $accFlag = 0;
if (@ARGV > 3) {
	if ($ARGV[3] eq "-accVer") {
		$accFlag = 1;
	}
}
my $filename1 = $ARGV[0];
my $filename2 = $ARGV[1];
# get the names of the files from the @ARGV variable that hold command line information

%file2data = ();
%file2dataAcc = ();
# this initializes a the hash "file2data" to be empty.  This isn't necessary but good practice.

open FILE2, $filename2 or die "Could not open $filename2\n";
# this opens $filename2 and assigns it to the filestream FILE2
# the "or die..." part is optional and will kill the program if it can't open the file

#this next part is a loop that will read in each line of the file
# <FILE2> will actually read a line and store it in a special variable named $_
# This part is confusing at first, but you get used to it in perl
my $hstr= '';
my @header = ();
while (<FILE2>) {
	$count++;
	# $_ now contains a line from the file
	chomp; #this removes the "\n" character from the end of the line
	s/\r//g; # this is needed if the file was made in windows - they end lines with "\r\n";

	# now we have a line of text without a newline character(s) at the end
	# Assuming its a tab separated file, lets split it into it's columns based on the tab "\t"
	my @line = split /\t/;

	if ($count <= $header) {
		my $h = '';
		for (my $i=1;$i<@line;$i++) {
			$h .= "\t$line[$i]";	
		}
		push(@header, $h);
		next;
	}

	# all of these commands will act on the default variable $_ - you could act on a different 
	# variable using: my @line = split /\t/, $var2;
	
	# we assume that the first element in the array is an identifier
	# we can remove that element using "shift", which removes the first element of an array 
	# and returns that value.  This function was made for this purpose
	my $id = shift @line;
	
	# finally lets store the data into a hash that associates the data with the id
	$file2data{$id} = \@line;
	if ($accFlag) {
		$id =~ s/\.\d+?$//;
		$file2dataAcc{$id} = \@line;
	}
}
# at this point there are no more lines in the file
close FILE2;
# close the FILE2 stream

#basically do the same thing with file1
open FILE1, $filename1 or die "Could not open $filename1\n";
$count = 0;
while (<FILE1>) {
	$count++;
	chomp;
	s/\r//g;
	if ($count <= $header) {
		print "$_" . $header[$count-1] . "\n";
		next;
	}
	my @line = split /\t/;
	my $id = shift @line;


	# now we want to try to merge them
	# lets assume we only want data that is common between each file
	# First step is to check if the data exists in the other file

	my @file2 = ();
	my $found = 0;
	next if (!defined($id));
	if (exists($file2data{$id})) {
		@file2 = @{$file2data{$id}};
		$found = 1;
	} elsif ($accFlag) {
		my $idAcc = s/\.\d+?$//;
		if (exists($file2dataAcc{$id})) {
			@file2 = @{$file2dataAcc{$id}};
			$found = 1;
		} elsif (exists($file2data{$idAcc})) {
			@file2 = @{$file2data{$idAcc}};
			$found = 1;
		} elsif (exists($file2dataAcc{$idAcc})) {
			@file2 = @{$file2dataAcc{$idAcc}};
			$found = 1;
		}
	}

	if ($found) {
		#good the id exists!
		# all we need to do now is print the data out
		print "$id"; #this prints the id to stdout i.e. the screen
		for (my $i=0;$i<@line;$i++) {
			print "\t$line[$i]";  #print each data item from the first file out
		}
		# at this point we've just repeated what we found in the 1st file
		# now time to add the data from the 2nd file
		for (my $i=0;$i<@file2;$i++) {
			print "\t$file2[$i]"
		}
		print "\n"; # need to add a new file character
	} else {
		#id isn't in the other dataset - just skip this one
		next;
	}
}
close FILE1;

#all done - could add "exit;" if you want, but not needed
