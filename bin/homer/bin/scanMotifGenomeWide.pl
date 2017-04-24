#!/usr/bin/env perl
use warnings;
use lib "/home/carlos/tools/homer/.//bin";
my $homeDir = "/home/carlos/tools/homer/./";


use POSIX;
use HomerConfig;


my $suffix = ".masked";
my $suffix2 = ".fa";

my $config = HomerConfig::loadConfigFile();

sub printCMD {
	print STDERR "\n\tUsage: scanMotifGenomeWide.pl <motif> <genome> [options]\n";
	print STDERR "\t\tPossible Genomes:\n";
	foreach(keys %{$config->{'GENOMES'}}) {
		print STDERR "\t\t$_\t$config->{'GENOMES'}->{$_}->{'org'}\t$config->{'GENOMES'}->{$_}->{'directory'}\n";
	}
	print STDERR "\t\t\t-- or --\n";
    print STDERR "\t\tCustom: provide the path to genome FASTA files (directory or single file)\n";
	print STDERR "\n\tOutput will be sent to stdout\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-5p (report positions centered on the 5' start of the motif)\n";
	print STDERR "\t\t-bed (format as a BED file, i.e. for UCSC upload)\n";
	print STDERR "\t\t\t-int (round motif scores to nearest integer, use if making bigBed file)\n";
	print STDERR "\t\t-homer1 (use the original homer)\n";
	print STDERR "\t\t-homer2 (use homer2 instead of the original homer, default)\n";
	print STDERR "\t\t-keepAll (keep ALL sites, even ones that overlap, default is to keep one)\n";
	print STDERR "\t\t-mask (search for motifs in repeat masked sequence)\n";
	print STDERR "\t\t-p <#> (Number of CPUs to use)\n";
	print STDERR "\n";
	exit;
}

if (@ARGV < 2) {
	printCMD();
}
$p5flag = 0;
$maskFlag = 0;
$homer2Flag = 1;
$bedFlag = 0;
$numCPUs = 1;
$intFlag = 0;
$keepAllFlag = 0;
for (my $i=2;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-5p') {
		$p5flag = 1;
		print STDERR "Outputing file centered on the 5' start of the motif\n";
	} elsif ($ARGV[$i] eq '-bed') {
		$bedFlag = 1;
	} elsif ($ARGV[$i] eq '-mask') {
		$maskFlag = 1;
	} elsif ($ARGV[$i] eq '-int') {
		$intFlag = 1;
	} elsif ($ARGV[$i] eq '-keepAll') {
		$keepAllFlag = 1;
	} elsif ($ARGV[$i] eq '-p') {
		$numCPUs = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-homer1') {
		$homer2Flag = 0;
		print STDERR "\tUsing original homer to scan for motifs\n";
	} elsif ($ARGV[$i] eq '-homer2') {
		$homer2Flag = 1;
		print STDERR "\tUsing homer2 to scan for motifs\n";
	} else {
		printCMD();
	}
}


$mfile = $ARGV[0];
$genome = $ARGV[1];
$genomeDir = "";
$customGenome = 0;
if (!exists($config->{'GENOMES'}->{$genome})) {
	$customGenome = 1;
	my $asdf = "";
	($genome,$genomeDir,$asdf) = HomerConfig::parseCustomGenome($genome);
} else {
	$genomeDir = $config->{'GENOMES'}->{$genome}->{'directory'};
}


$size = 10000;

$tmpFile = rand() . ".tmp";
$tmpFile2 = $tmpFile . ".2";

if ($customGenome==1 && -f $genomeDir) {
	`ls -1 "$genomeDir" > "$tmpFile"`;
} else {
	`ls -1 "$genomeDir"/*fa* > "$tmpFile"`;
}
open IN, $tmpFile;
@files = ();
while (<IN>) {
	chomp;
	s/\r//g;
	push(@files, $_);
}
close IN;

%idCounts = ();
$chr = '';
#print STDERR "Seqfiles: @files\n";


foreach(@files) {
	my $file = $_;
	#$file =~ s/(\/.*)?$//;
	#$file .= "/bioinformatics/homer/data/genomes/mm9/chr11.fa";
	open IN, $file or die "Couldn't open $file\n";
	my $position = 0;
	my $totalLength = 0;
	$chr = '';
	my $curSeq = '';
	my $startPos = 0;
	my $justPrinted = 0;

	#open SEQFILE, ">$tmpFile2";

	while (<IN>) {
		chomp;
		s/\r//g;
		if (/^>/) {
			/^>(.*?)$/;
			my $nextchr = $1;

			if ($chr ne '' && $curSeq ne '') {
				print SEQFILE "$chr-$startPos\t$curSeq\n";		
				close SEQFILE;
				processChr();
			}

			$chr = $nextchr;
			print STDERR "\n\tProcessing $chr\n";
			$position = 0;
			$totalLength = 0;
			$curSeq = '';
			$startPos = 0;
			open SEQFILE, ">$tmpFile2";
			next;
		}
		if ($maskFlag==1) {
			s/[acgt]/N/g;
		} else {
			s/a/A/g;
			s/c/C/g;
			s/g/G/g;
			s/t/T/g;
		}

		my $len = length($_);
		$totalLength += $len;
		$curSeq .= $_;
		if ($totalLength > $size) {
			print SEQFILE "$chr-$startPos\t$curSeq\n";		
			$totalLength = $len;
			$startPos = $position;
			$curSeq = $_;
			$justPrinted = 1;
		} else {
			$justPrinted = 0;
		}

		$position += $len;

	}
	if ($justPrinted == 0) {
		print SEQFILE "$chr-$startPos\t$curSeq\n";		
	}
	close IN;
	close SEQFILE;
	processChr();
}

`rm -f "$tmpFile" "$tmpFile2"`;
exit;

sub processChr {

	if ($homer2Flag) {
		`homer2 find -s "$tmpFile2" -m "$mfile" -offset 0 -p $numCPUs > "$tmpFile"`;
	} else {
		`homer -s "$tmpFile2" -a FIND -m "$mfile" > "$tmpFile"`;
	}

	open INN, "$tmpFile";
	my %pos = ();
	my @sites = ();
	while (<INN>) {
		chomp;
		my @line = split /\t/;

		$line[0]=~ /^(.*?)\-(\d+)$/;
		my $chr = $1;
		my $gpos = $2;
		my $pos = $line[1];
		my $seq = $line[2];

		my $d = 0;
		my $name = "";
		my $score = "";
		if ($homer2Flag) {
			$name = $line[3];
			$d = $line[4];
			$score = $line[5];
		} else {
			$d = $line[4];
			$name = $line[5];
			$score = $line[6];
		}

		my $start = $gpos+$pos+1;
		my $end = $start + length($seq)-1;
		if ($homer2Flag && ($d eq '-' || $d eq '1') ) {
			$start -= length($seq)-1;
			$end -= length($seq)-1;
		}
		my $mid = floor(($start+$end)/2);

		my $pd = $name . "-" . $chr . "-" . $start . "-" . $d;
		next if (exists($pos{$pd}));
		$pos{$pd}= 1;

		if (!exists($idCounts{$name})) {
			$idCounts{$name} = 1;
		}

		if ($p5flag==1) {
			my $center = $start;
			if ($d eq '-' || $d eq '1') {
				$center = $end;
			}
			$start = $center -100;
			$end = $center +100;
		}
		
		if ($start > $end) {
			next;
		}
		my $ss = {s=>$start,e=>$end,d=>$d,seq=>$seq,m=>$mid,n=>$name,ss=>$score,c=>$chr};
		push(@sites, $ss);
	}
	close INN;
	@sites = sort {$a->{'c'} cmp $b->{'c'} || 
						$a->{'s'} <=> $b->{'s'}} @sites;

	my $Nsites = scalar(@sites);
	my $removed= 0;

	for (my $i=0;$i<@sites;$i++) {
		my $m= $sites[$i]->{'m'};
		my $bad = 0;
		my $mlen = length($sites[$i]->{'seq'});
		if ($keepAllFlag==0) {
			for (my $j=$i-1;$j>=0;$j--) {
				last if ($m - $sites[$j]->{'m'} > $mlen/2);
				if ($sites[$i]->{'n'} eq $sites[$j]->{'n'} && 
						$sites[$i]->{'c'} eq $sites[$j]->{'c'} && 
							$sites[$i]->{'ss'} < $sites[$j]->{'ss'}) {
					$bad = 1;
					last;
				}
			}
			for (my $j=$i+1;$j<@sites;$j++) {
				last if ($sites[$j]->{'m'}-$m > $mlen/2);
				if ($sites[$i]->{'n'} eq $sites[$j]->{'n'} &&
						$sites[$i]->{'c'} eq $sites[$j]->{'c'} && 
							$sites[$i]->{'ss'} <= $sites[$j]->{'ss'}) {
					$bad = 1;
					last;
				}
			}
		}
		if ($bad == 0) {
			

			my $chr = $sites[$i]->{'c'};
			my $start = $sites[$i]->{'s'};
			my $end = $sites[$i]->{'e'};
			my $dir = $sites[$i]->{'d'};
			my $score = $sites[$i]->{'ss'};
			my $name = $sites[$i]->{'n'};
			my $seq = $sites[$i]->{'seq'};
			if ($start < 0) {
				#print STDERR "\t!!! Error, start=$start < 0\n";
				#print STDERR "\t\t$chr\t$start\t$end\t$name\t$score\t$dir\n";
				next;
			}
			if ($start >= $end) {
				#print STDERR "\t!!! Error, start>=end: $start $end\n";
				#print STDERR "\t\t$chr\t$start\t$end\t$name\t$score\t$dir\n";
				next;
			}

			if ($bedFlag) {
				$name =~ s/\/.*$//g;
				$name =~ s/\s/_/g;
				if ($intFlag) {
					$score = floor($score+0.5);
				}
				print "$chr\t$start\t$end\t$name\t$score\t$dir\n";
			} else {
				my $id = $idCounts{$name};
				print "$name-$id\t$chr\t$start\t$end\t$dir\t$score\t$seq\n";
				$idCounts{$name}++;
			}
		} else {
			$removed++;
		}
	}
	print STDERR "\t$chr - total: $Nsites (removed: $removed)\n";
}

