#!/usr/bin/env perl
use warnings;



if (@ARGV < 1) {
	print STDERR "\n\t[options] <s_X_eland_result.txt or BED file> [file 2] ...\n";
	print STDERR "\tprogram will automatically detect eland or BED or bowtie format...\n";
	print STDERR "\tEland Options: (optional)\n";
	print STDERR "\t\t-mis <0,1,2> (Max mismatches, default=2)\n";
	print STDERR "\t\t-len # (only consider mismatches within the first # bp of sequence [for eland_extended]-def: 32)\n";
	print STDERR "\t\t-uniq <#> (Max matches to genome for multiple eland format, default=1)\n";
	print STDERR "\t\t-seq (include sequence with tag)\n";
	print STDERR "\t\t-nonuniq (include a single version of nonuniq tags, bowtie only)\n";
	print STDERR "\n\tEland format expects the following columns:\n";
	print STDERR "\t\tColumn 1: SpotID\n";
	print STDERR "\t\tColumn 2: Sequence\n";
	print STDERR "\t\tColumn 3: MatchInfo (i.e. U0,NM,R2,etc.)\n";
	print STDERR "\t\tColumn 4: # of locations with perfect match\n";
	print STDERR "\t\tColumn 5: # of locations with 1 mis-match\n";
	print STDERR "\t\tColumn 6: # of locations with 2 mis-match\n";
	print STDERR "\t\tColumn 7: chromosome (i.e. chr1 or chr1.fa)\n";
	print STDERR "\t\tColumn 8: chromosome position\n";
	print STDERR "\t\tColumn 9: F or R for stand\n";
	print STDERR "\t\t...\n";
	print STDERR "\tMultiple Eland Format\n";
	print STDERR "\t\tColumn 1: SpotID\n";
	print STDERR "\t\tColumn 2: Sequence\n";
	print STDERR "\t\tColumn 3: MatchInfo (i.e. NM, or 2:12:255)\n";
	print STDERR "\t\tColumn 4: positions (i.e. chr:positionF0,chr:positionR2,... etc.)\n";
	print STDERR "\tBED format expects the following columns:\n";
	print STDERR "\t\tColumn 1: chromosome (looks for the line to start with \"chr\")\n";
	print STDERR "\t\tColumn 2: chromosome position start\n";
	print STDERR "\t\tColumn 3: chromosome position end\n";
	print STDERR "\t\tColumn 4: + or - for stand\n";
	print STDERR "\t\t -or-\n";
	print STDERR "\t\tColumn 6: + or - for stand\n\n";
	exit;
}

my %stats = ();
my $maxMatches = 1;
my $maxMisMatches = 2;
my $maxLength = 32;
my $seqFlag = 0;
my $nonuniqFlag = 0;
my @files = ();
for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-mis') {
		$maxMisMatches = $ARGV[++$i];
		print STDERR "\tAllowing up to $maxMisMatches mismatches in tags\n";
	} elsif ($ARGV[$i] eq '-nonuniq') {
		$nonuniqFlag = 1;
		print STDERR "\tkeeping non-unique tags too\n";
	} elsif ($ARGV[$i] eq '-seq') {
		$seqFlag = 1;
		print STDERR "\tKeeping sequence\n";
	} elsif ($ARGV[$i] eq '-len') {
		$maxLength = $ARGV[++$i];
		print STDERR "\tWill only consider mismatches through base $maxLength\n";
	} elsif ($ARGV[$i] eq '-uniq') {
		$maxMatches = $ARGV[++$i];
		print STDERR "\tAllowing tags to match genome up to  $maxMatches times\n";
	} else {
		push(@files, $ARGV[$i]);
		print STDERR "\tWill parse file: $ARGV[$i]\n";
	}
}





my $bowtieLastID = '';
my $bowtieLastMis = 1e100;
my $bowtieDone = 0;
my $bowtieOutput = '';



my $id = 0;
foreach(@files) {
	my $file = $_;
	open IN, $file;
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		my $dir = 0;
		my $pos = 0;
		my $v = 0;
		my $seq = '';

		# check for bowtie format
		if ($line[1] eq '-' || $line[1] eq '+') {


			if ($line[0] eq $bowtieLastID) {
				next if ($bowtieDone);
				
				my $curMis = 0;
				if (@line > 7) {
					my @a = split /\,/,$line[7];
					$curMis = scalar(@a);
				}
				if ($curMis == $bowtieLastMis) {
					$bowtieOutput = '' if ($nonuniqFlag == 0);
				}
				$bowtieDone = 1;
			} else {

				print $bowtieOutput if ($bowtieOutput ne '');
				$bowtieDone = 0;

				$bowtieLastMis = 0;
				if (@line > 7) {
					my @a = split /\,/,$line[7];
					$bowtieLastMis = scalar(@a);
				}

				$dir=0;
				$dir = 1 if ($line[1] eq '-');
				$chr = $line[2];
				$seq = $line[4];
				$pos = $line[3]+1;
				$pos = $pos += length($seq)-1 if ($dir == 1);

				$bowtieOutput = "$id\t$chr\t$pos\t$dir\t1";
				$bowtieOutput .= "\t$seq" if ($seqFlag==1);
				$bowtieOutput .= "\n";
				$bowtieLastID = $line[0];

			}
			$id++;
			next;
		}

		#check for BED format
		if ($line[0] =~ /^chr/) {
			my $pos = $line[1];
			$v = 1;
			if (@line > 3) {
				if ($line[3] ne '-' && $line[3] ne '+') {
					if (@line > 5 && $line[5] eq '-') {
						$dir = 1;
						#$pos = $line[2]+24;
						$pos = $line[2];
					}
				} elsif ($line[3] eq '-') {
					$dir = 1;
					$pos = $line[2];
				}
				if (@line > 4 && $line[4] =~ /^\d+$/) {
					$v = $line[4];
					if ($v ==0) {
						$v =1;
					}
					if (0) {
						if ($v > 100) {
							print STDERR "Warning: 5th column has value of $v\n";
						}
						if ($v == 0) {
							print STDERR "5th columan has value of zero\n";
						}
					}
				}
			}
    		print "$id\t$line[0]\t$pos\t$dir\t$v\n";
    		$id++;
			next;
		} else {
		#eland format

			$id = $line[0];

			if ($line[2] =~ /\:/) {
				#multiple eland format or eland_extended
				my @stat = split /\:/, $line[2];
				my $bestMis = -1;
				my $mode = 'NM';
				if ($stat[0] > 0) {
					if ($stat[0] <= $maxMatches) {
						$bestMis = 0;
						$mode = 'U0';
					} else {
						$mode = 'R0';
					}
				} elsif ($stat[1] > 0 && $maxMisMatches > 0) {
					if ($stat[1] <= $maxMatches) {
						$bestMis = 1;
						$mode = 'U1';
					} else {
						$mode = 'R1';
					}
				} elsif ($stat[2] > 0 && $maxMisMatches > 1) {
					if ($stat[2] <= $maxMatches) {
						$bestMis = 2;
						$mode = 'U2';
					} else {
						$mode = 'R2';
					}
				}
				$stats{$mode}++;
				next if ($bestMis == -1);
				my @pos = split /\,/, $line[3];
				my $chr = '';
				my $seq = $line[1];
				my $seqLen = length($line[1]);
				my $curMaxLength = $maxLength;
				if ($curMaxLength > $seqLen) {
					$curMaxLength = $seqLen;
				}
				my $pcount = 0;
#print "===============================\n";
				foreach(@pos) {
					$pcount++;
					my $p = $_;
					my $pos = 0;
					my $dir = 0;
					my $good = 0;
					my $goodMis = 0;
					if (/\:/) {
						my @p = split /\:/;
						$chr = $p[0];
						$p = $p[1];
					}
					if ($p =~ /^(\d+)([FR])(\d)$/) {
						#eland multiple format
						$pos = $1;
						$dir = 0;
						my $curMis = $3;
						if ($2 eq 'R') {
							$dir = 1;
							$pos = $pos + $seqLen-1;
						}
						if ($curMis <= $bestMis) {
							$good=1;
						}
					} elsif ($p =~ /^(\d+)([FR])(.+)$/) {
						$pos = $1;
						$dir = 0;
						if ($2 eq 'R') {
							$dir = 1;
							$pos = $pos + $seqLen-1;
						}
						my @matches = split(/[ACGT]+/,$3);
						my @mis = split(/\d+/,$3);
						my $curLen = 0;
						my $curMis = 0;
						my $misIndex=0;
#print "------------------\n";
#print "@matches\n";
#print "@mis\n";
#print "-----------------\n";
						for (my $k=0;$k<@matches;$k++) {
							if ($matches[$k] ne '') {
								$curLen += $matches[$k];
#print "\t+ $matches[$k] Index: $k\t$misIndex\n";
								$misIndex++ if ($k==0);
							}
							if ($curLen >= $curMaxLength) {
#print "P $curLen\t$curMis\n";
								if ($curMis <= $bestMis) {
									$good=1;
									$goodMis = $curMis;
								}
								last;
							}
							next if ($misIndex >= scalar(@mis));
							my $L = length($mis[$misIndex]);
							$misIndex++;

							my $curDiff = $curMaxLength - $curLen;
							if ($curDiff < $L) {
#print "\t+ $L MM ($curDiff) Index:$k\t$misIndex\n";
								$curMis += $curDiff;
								$curLen += $curDiff;
							} else {
#print "\t+ $L MM ($curDiff) Index:$k\t$misIndex\n";
								$curMis += $L;
								$curLen += $L;
							}
							
							if ($curLen >= $curMaxLength) {
#print "M $curLen\t$curMis\n";
								if ($curMis <= $bestMis) {
									$good=1;
									$goodMis = $curMis;
								}
								last;
							}
						}	
					}
					if ($good == 1) {
						$chr =~ s/mm_ref_//;
						$chr =~ s/\.fa\.masked//;
						$chr =~ s/\.fa//;
						print "$id\t$chr\t$pos\t$dir\t1";
						print "\t$seq" if ($seqFlag);
						print "\n";
						#print "GOOD(mode:$mode $bestMis $goodMis $pcount): @line\n";
					} else {
						#print "BAD(mode:$mode $bestMis $goodMis $pcount): @line\n";
					}
				}
				next;
			} elsif (@line > 20) {
				#export format
				my $chr = $line[10];
				if ($chr eq 'QC') {
					$stats{'QC'}++;
					next;
				} elsif ($chr eq 'NM') {
					$stats{'NM'}++;
					next;
				} elsif ($chr =~ /\d+\:\d+\:\d+/) {
					$stats{'R012'}++;
					next;
				}
				$stats{'U012'}++;
				my $seq = $line[8];
				$chr =~ s/\.fa\.masked//;
				$chr =~ s/\.fa//;
				my $p = $line[12];
				if ($line[13] eq 'R') {
					$dir = 1;
					my $L = length($seq) - 1;
					$p+=$L;
				}
				print "$id\t$chr\t$p\t$dir\t1\n";
			} elsif (@line > 13 && ($line[13] eq 'F' || $line[13] eq 'R')) {
				my $chr = $line[10];
				$chr =~ s/\.fa.*$//;
				next if ($chr eq '');
				my $seq = $line[8];
				my $p = $line[12];
				my $id = $line[0];
				my $dir= 0;
				if ($line[13] eq 'R') {
					$dir = 1;
					my $L = length($seq) - 1;
					$p+=$L;
				}
				print "$id\t$chr\t$p\t$dir\t1\n";
			} else {
				#traditional eland format
				$stats{$line[2]}++;
				next if (@line < 9);
				next if ($line[2] eq 'U1' && $maxMisMatches < 1);
				next if ($line[2] eq 'U2' && $maxMisMatches < 2);
	
				my $chr = $line[6];
				my $seq = $line[1];
				$chr =~ s/^.*\///;
				$chr =~ s/\.fa\.masked//;
				$chr =~ s/\.fa//;
	#$chr = "chr" . $chr;	
				my $p = $line[7];
				if ($line[8] eq 'R') {
					$dir = 1;
					my $L = length($seq) - 1;
					$p+=$L;
				}
				print "$id\t$chr\t$p\t$dir\t1";
				print "\t$seq" if ($seqFlag);
				print "\n";
			}
		}
	}
	close IN;
}

foreach(keys %stats) {
	print STDERR "$_\t$stats{$_}\n";
}
