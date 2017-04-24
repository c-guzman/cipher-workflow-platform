#!/usr/bin/perl -w
if (@ARGV < 1) {
	print STDERR "<fasta file>\n";
	exit;
}
open IN, $ARGV[0];
my $prefixx = "P";

my $good = 1;

$count = 0;
my $ZZ = 0;
while (<IN>) {
	$ZZ++;
	chomp;
	s/\r//g;
	next if (/^#/);
	next if ($_ eq '');
	if (/>/) {
		if ($count>0 && $good==1) {
			print "\n";
		}
		$good = 1;


		my $n = '';


		if (0) {
#			/^>gi\|.*?\|ref\|(.*?)\./;
			if (/^>(.*?)\s/) {
				$n = $1;
			} elsif(/^>(.*?)$/) {
				$n = $1;
			} else {
				$n= $_;
			}
		}
		if (0) {
			/^>(.*?) .+DIRECTION: (.+?)$/;
			if ($2 eq 'rev') {
				$good = 0;
				next;
			} else {
				$good = 1;
			}
		}
		#/^>.+\|ti\|(\d*?) /;
		#/^>.+? (.+?) /;
		#/Plasmodium_falciparum\|(.+?)\|/;
		# JGI format
		if (0) {
			/^>.*? (.*?)$/;
			$n = $1;
			if ($n =~  /fwd/) {
				$good = 0;
				next;
			} else {
				$good = 1;
			}
		}

		# Unigene Format
		#/\/gb=(.*?) /;
		#$n = $1;

		#ncbi protein
		#/^>.*?ref\|(.*?)\|/;
		#$n = $1;
		#$n =~ s/\..*$//;

		#JGI/EMBL protein
		#/^>(.*?) /;
		#/^>(.*?)$/;
		#/^>(.*?)\|Anno/;
		#$n = $1;
	
		#genbank dbEST format
		if (0) {
			#/^>.*? (.*?) /;
			/^>.*?gb\|(.*?)\|/;
			$n = $1;
			if (/ 3\'/) {
				$good = 0;
				next;
			}
			if (/ 5\'/) {
				$good = 1;
			} else {
				$good = 0;
				next;
			}
		}
		
		#/^>.*\|(.*?) /;
		#$n = $1;
		
	
		# BASIC
		#/^>.*? (.*?)$/;
		/^>(.*?)\s/;
		#/^>jgi\|Monbr1\|(.*?)\|/;
		$n = $1;

		$n =~ s/\s*$//;
		if ($n eq '') {
			$n = $prefixx . "$count";
		}
		
		print "$n\t";
		$count++;
	}
	else {
		next if ($good == 0);
		s/a/A/g;
		s/c/C/g;
		s/g/G/g;
		s/t/T/g;
		s/[^ACGT]/N/g;
		print $_;
	}
}
close IN;
print "\n";
print STDERR "\tFound $count sequences\n";
