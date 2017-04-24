#!/usr/bin/env perl
use warnings;
use lib "/home/carlos/tools/homer/.//bin";
my $homeDir = "/home/carlos/tools/homer/./";
my $promoterSeqOffset = -2000;


use POSIX;

use Statistics;

my $analysisType = "GETPVALUE";

my $knownFile  = $homeDir . "/data/knownTFs/known.motifs";
my $pvalueThresh = 0.01;
my $cpus = 1;
my $cache = 500;
my $strand = "both";
my $stat = "hypergeo";
my $seqFile = '';
my $groupFile = '';
my $directory = '';
my $homer2Flag =0;
my $bits = "-S";
my $dbview = 0;


if (@ARGV < 3) {
	printCMD();
}

sub printCMD {
	print STDERR "\n\tUsage: findKnownMotifs.pl -s <seq file> -g <group file> -o <Output Directory> [options]\n";
	print STDERR "\n\t!!! Normally this is called by the findMotifs.pl/findMotifsGenome.pl programs\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-s <seq file> (tsv sequence file)\n";
	print STDERR "\t\t-g <group file> (tsv group file)\n";
	print STDERR "\t\t-o <output directory>\n";
	print STDERR "\t\t-m <motif file> (Known motif file, default: $knownFile)\n";
	print STDERR "\t\t-strand <both|+|->\n";
	print STDERR "\t\t-stat <hypergeo|binomial> (default: hypergeo)\n";
	print STDERR "\t\t-pvalue <#> (p-value cutoff, default: $pvalueThresh)\n";
	print STDERR "\t\t-optimize (Optimize motif degeneracy thresholds to maximize enrichment)\n";
	print STDERR "\t\t-homer (use original homer, default, for now...)\n";
	print STDERR "\t\t-homer2 (use homer2)\n";
	print STDERR "\t\t-p <#> (number of CPUs to use, homer2 only)\n";
	print STDERR "\t\t-cache <#> (size in MB of statistics cache, default: 500)\n";
	print STDERR "\t\t-bits (scale logos by information content)\n";
	print STDERR "\t\t-dbview (internal option)\n";
	print STDERR "\n";
	exit;
}
for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-s') {
		$seqFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-g') {
		$groupFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-dbview') {
		$dbview = 1;
	} elsif ($ARGV[$i] eq '-o') {
		$directory = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-bits') {
		$bits = "";
	} elsif ($ARGV[$i] eq '-m') {
		$knownFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-pvalue') {
		$pvalueThresh = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-optimize') {
		$analysisType = "OPTPVALUE";
	} elsif ($ARGV[$i] eq '-homer') {
		$homer2Flag = 0;
	} elsif ($ARGV[$i] eq '-homer2') {
		$homer2Flag = 1;
	} elsif ($ARGV[$i] eq '-stat') {
		$stat = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-p') {
		$cpus = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-cache') {
		$cache = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-strand') {
		$strand = $ARGV[++$i];
	} else {
		print STDERR "!!! Didn't understand \"$ARGV[$i]\"\n";
		printCMD();
	}
}
if ($seqFile eq '' || $groupFile eq '' || $directory eq '') {
	print STDERR "!!! Missing one or more of -s/-g/-o\n";
	printCMD();
}

my $knownDir = "knownResults/";
$scoreMethod ='absDifference';
$height = 2.0;
$widthFactor = 1.0;

`mkdir -p "$directory"`;
`mkdir -p "$directory/$knownDir"`;

my $tmpFile = "$directory/" . rand() . ".motif.tmp";
my $tmpFile2 = "$directory/" . rand() . ".motif.2.tmp";

my $options = '';
if ($homer2Flag) {

	$options .= " -strand $strand" if ($strand ne 'both');
	$options .= " -opt " if ($analysisType eq 'OPTPVALUE');
	#$options .= " -opt \"$tmpFile2\"" if ($analysisType eq 'OPTPVALUE');
	$options .= " -stat $stat";
	$options .= " -p $cpus";
	
	`homer2 known -s "$seqFile" -g "$groupFile" -m "$knownFile" $options -mout "$tmpFile" -cache $cache > "$tmpFile2"`;

} else {
	if ($strand eq '+') {
		$options .= " -norevopp ";
	}
	`homer -s "$seqFile" -g "$groupFile" -m "$knownFile" -a $analysisType $options > "$tmpFile"`;
}


my $mainpageName = "knownResults.html";
my $txtFile = "knownResults.txt";

my $motifs = readMotifFile($tmpFile);
my $totalSeq=0;
my $totalTargets=0;
my $aAA=0;
my $bBB=0;
my $totalNumTargets=0;
my $totalNumBackground=0;
my $numTargets=0;
my $numBackground=0;
my $percentTargets=0;
my $percentBackground=0;
my $mpvalue=0;

if ($homer2Flag) {
	($numTargets,$percentTargets,$numBackground,$percentBackground, $mpvalue) = parseStr2($motifs->[0]->{'str'});
	$percentTargets =~ s/\%//;
	$percentBackground =~ s/\%//;
	$totalNumTargets = floor($numTargets/($percentTargets*0.01)+0.5);
	$totalNumBackground = floor($numBackground/($percentBackground*0.01)+0.5);
} else {
	($totalSeq, $totalTargets, $aAA, $bBB) = parseStr($motifs->[0]->{'str'});
}

my @motifs = sort {$a->{'logp'} <=> $b->{'logp'}} @$motifs;

open TXT, ">$directory/$txtFile";
print TXT "Motif Name\tConsensus\tP-value\tLog P-value\tq-value (Benjamini)";
if ($homer2Flag) {
	print TXT "\t# of Target Sequences with Motif(of $totalNumTargets)\t% of Target Sequences with Motif";
	print TXT "\t# of Background Sequences with Motif(of $totalNumBackground)\t% of Background Sequences with Motif\n";
} else {
	print TXT "\t# of Genes\t# of regulated Genes\t# of Genes with Motif\t#of regulated genes with Motif\n";
}

open MAIN, ">$directory/$mainpageName";

print MAIN "<HTML><HEAD><TITLE>$directory - Homer Known Motif Enrichment Results</TITLE></HEAD>\n";
print MAIN "<BODY>\n";
	
print MAIN "<H1>Homer Known Motif Enrichment Results ($directory)</H1>\n";
print MAIN "<A HREF=\"homerResults.html\">Homer <i>de novo</i> Motif Results</A><BR/>\n";
print MAIN "<A HREF=\"geneOntology.html\">Gene Ontology Enrichment Results</A><BR/>\n";
print MAIN "<A HREF=\"$txtFile\">Known Motif Enrichment Results (txt file)</A><BR/>\n";
if ($homer2Flag) {
	print MAIN "Total Target Sequences = $totalNumTargets, Total Background Sequences = $totalNumBackground</BR>\n";
} else {
	print MAIN "Total Sequences = $totalSeq, Total Targets = $totalTargets</BR>\n";
}
print MAIN "<TABLE border=\"1\" cellpading=\"0\" cellspacing=\"0\">\n";
print MAIN "<TR><TD>Rank</TD><TD>Motif</TD><TD>Name</TD><TD>P-value</TD><TD>log P-pvalue</TD><TD>q-value (Benjamini)</TD>";
if ($homer2Flag) {
	print MAIN "<TD># Target Sequences with Motif</TD><TD>% of Targets Sequences with Motif</TD>";
	print MAIN "<TD># Background Sequences with Motif</TD><TD>% of Background Sequences with Motif</TD>";
} else {
	print MAIN "<TD>Total Sequences with Motif</TD><TD>Targets Sequences with Motif</TD>";
}
print MAIN "<TD>Motif File</TD>\n";
print MAIN "<TD>PDF</TD></TR>\n";

my $motifNum = 1;
my $ZZ = scalar(@motifs);
print STDERR "\tPreparing HTML output with sequence logos...\n";
foreach(@motifs) {
	my $matrixFile = "known$motifNum.motif";
	my $logoFile = "known$motifNum.logo";
	print TXT "$_->{'name'}\t$_->{'cons'}\t";


	$_->{'pvalue'} = exp($_->{'logp'});
	my $pvalue = sprintf("%.3e", $_->{'pvalue'});
	my $lp = sprintf("%.3e", $_->{'logp'});
	my $fdr = sprintf("%.4f", $_->{'fdr'});

	if ($homer2Flag) {
		($numTargets,$percentTargets,$numBackground,$percentBackground,$mpvalue) = parseStr2($_->{'str'});
		$pvalue = $mpvalue;
		print TXT "$pvalue\t$lp\t$fdr\t$numTargets\t$percentTargets";
		print TXT "\t$numBackground\t$percentBackground\n";

	} else {
		my $cCC =0;
		my $dDD = 0;
		($cCC, $dDD,$totalMotif,$totalTargetMotif) = parseStr($_->{'str'});
		print TXT "$pvalue\t$lp\t$fdr\t$totalSeq\t";
		print TXT "$totalTargets\t$totalMotif\t";
		print TXT "$totalTargetMotif\n";
	}

	next if ($_->{'pvalue'} > $pvalueThresh);

	print STDERR "\t\t$motifNum of $ZZ ($pvalue) $_->{'name'}\n";

	print MAIN "<TR><TD>$motifNum</TD><TD><IMG src=\"$knownDir/$logoFile.png\"/></TD>";
	print MAIN "<TD>$_->{'name'}</TD>";
	print MAIN "<TD>$pvalue</TD><TD>$lp</TD><TD>$fdr</TD>";
	if ($homer2Flag) {
		print MAIN "<TD>$numTargets</TD><TD>$percentTargets</TD>";
		print MAIN "<TD>$numBackground</TD><TD>$percentBackground</TD>";
	} else {
		print MAIN "<TD>$totalMotif</TD><TD>$totalTargetMotif</TD>";
	}
	print MAIN "<TD><A target=\"_blank\" HREF=\"$knownDir/$matrixFile\">motif file (matrix)</A></TD>\n";
	print MAIN "<TD><A HREF=\"$knownDir/$logoFile.pdf\">pdf</A></TD></TR>\n";
	
	printMotif($_, $directory . "/$knownDir/$matrixFile");
	`profile2seq.pl "$directory/$knownDir/$matrixFile" 100 > "$tmpFile"`;
	my $width = $widthFactor * $_->{'len'};
	`seqlogo -a $bits -f "$tmpFile" -F PNG -S -c -o "$directory/$knownDir/$logoFile" -h $height -w $width`;
	`seqlogo -a $bits -f "$tmpFile" -F PDF -S -c -o "$directory/$knownDir/$logoFile" -h $height -w $width 2> /dev/null`;

	$motifNum++;
}
print MAIN "</TABLE></BODY></HTML>\n";
close MAIN;
close TXT;

`rm -f "$tmpFile" "$tmpFile2"`;
exit;

sub parseStr2 {
	my ($str) = @_;
	my @str = split(/\)\,/,$str);
	my @t = split(/[\:\(]/,$str[0]);
	my $numTargets = $t[1];
	my $percentTargets = $t[2];
	my @b = split(/[\:\(]/,$str[1]);
	my $numBackground = $b[1];
	my $percentBackground = $b[2];
	my @c = split(/[\:]/,$str[2]);
	my $pvalue = $c[1];
	return ($numTargets, $percentTargets, $numBackground, $percentBackground, $pvalue);
}
sub parseStr {
	my ($str) = @_;
	my $N = 0;
	my $P = 0;
	my $M = 0;
	my $nm = 0;
	my @str = split /\,/, $str;
	if (@str > 3) {
		$N = $str[0];
		$P = $str[1];
		$M = $str[2];
		$nm = $str[3];
	}
	return ($N, $P, $M, $nm);
}

sub matchMotif {
	my ($motif, $known) = @_;

	my @results = ();
	for (my $i=0;$i<@$known;$i++) {
		my ($score, $offset, $dir) = compareMotifs($motif, $known->[$i]);
		my $res = {motif=>$known->[$i], score=>$score, offset=>$offset, dir=>$dir};
		push(@results, $res);
	}
	@results = sort {$b->{'score'} <=> $a->{'score'}} @results;

	return \@results;
}

sub reduceMotifs {
	my ($motifs, $threshold) = @_;
	my %mask = ();	
	my @m = sort {$a->{'logp'} <=> $b->{'logp'}} @$motifs;
	my $startNum = scalar(@m);
	my @newList = ();
	for (my $i=0;$i<@m;$i++) {
		next if ($m[$i]->{'mask'} ==1);
		push(@newList, $m[$i]);
		for (my $j=$i+1;$j<@m;$j++) {
			next if ($m[$j]->{'mask'} ==1);
			my ($s,$o,$d) = compareMotifs($m[$i], $m[$j]);
			if ($s > $threshold) {
				$m[$j]->{'mask'} = 1;
			}
		}
	}
	my $endNum = scalar(@newList);
	print STDERR "$startNum motifs -> $endNum motifs\n";
	return \@newList;
}

sub printMotif {
	my ($motif, $file) = @_;
	open MOTIF, ">$file";
	print MOTIF ">$motif->{'cons'}\t$motif->{'name'}\t$motif->{'v'}\t$motif->{'logp'}\t$motif->{'gapinfo'}\t$motif->{'str'}\n";
	my $matrix = $motif->{'matrix'};
	foreach(@{$matrix}) {
		for (my $i=0;$i<4;$i++) {
			print MOTIF "\t" if ($i>0);
			print MOTIF "$_->[$i]";
		}
		print MOTIF "\n";
	}
}

sub compareMotifs  {
	my ($m1, $m2) = @_;

	my $bestOffset = 0;
	my $bestDirection = 0;
	my $bestScore = -1e10;

	my $rv2 = revoppMotif($m2);

	my $len1 = $m1->{'len'};
	my $len2 = $m2->{'len'};
	my $offset2 = $len2-1;
	my $offset1 = 0;
	for (my $offset2 = $len2-1;$offset2>=-1;$offset2--) {
		my $officialOffset = -1*$offset2+$offset1;

		my $max1 = $len1 - $offset1;
		my $max2 = $len2 - $offset2;
		my $curLen = $max1;
		if ($max2 < $curLen) {
			$curLen = $max2;
		}
		my @mm1 = ();
		for (my $j=$offset1;$j<$offset1+$curLen;$j++) {
			push(@mm1, $m1->{'matrix'}->[$j]);
		}
		my @mm2 = ();
		for (my $j=$offset2;$j<$offset2+$curLen;$j++) {
			push(@mm2, $m2->{'matrix'}->[$j]);
		}
		my @mm2r = ();
		for (my $j=$offset2;$j<$offset2+$curLen;$j++) {
			push(@mm2r, $rv2->{'matrix'}->[$j]);
		}
		my $score = scoreComparison(\@mm1, \@mm2);
		my $rvScore = scoreComparison(\@mm1, \@mm2r);
#print STDERR "$offset1 ($len1)\t$offset2 ($len2)\t$officialOffset\t$score\t$rvScore\n";
		if ($rvScore > $score) {
			if ($rvScore > $bestScore) {
				$bestScore = $rvScore;
				$bestOffset = $officialOffset;
				$bestDirection = 1;
			}
		} else {
			if ($score > $bestScore) {
				$bestScore = $score;
				$bestOffset = $officialOffset;
				$bestDirection = 0;
			}
		}
		if ($offset2 == 0) {
			$offset2++;
			$offset1++;
			if ($offset1 >= $len1) {
				last;
			}
		}
	}
	return ($bestScore, $bestOffset, $bestDirection);
}

sub scoreComparison {
	my ($m1, $m2) = @_;
	my $len = @$m1;
	my $score = 0;
	if ($scoreMethod eq 'absDifference') {
		for (my $i=0;$i<$len;$i++) {
			for (my $j=0;$j<4;$j++){ 
				$score += $m1->[$i][$j] * $m2->[$i][$j] - 0.0625;
			}
		}
	}
	return $score;
}

sub revoppMotif {
	my ($m) = @_;
	my @a = ();
	my $rv = {matrix=>\@a,name=>$m->{'name'},cons=>$m->{'cons'},
					v=>$m->{'v'},p=>$m->{'p'},logp=>$m->{'logp'},str=>$m->{'str'},
					gapinfo=>$m->{'gapinfo'}, len=>$m->{'len'},mask=>$m->{'mask'}};
	for (my $i=-1+$m->{'len'};$i>=0;$i--){ 
		my @b = ();
		for (my $j=3;$j>=0;$j--) {
			push(@b, $m->{'matrix'}->[$i][$j]);
		}
		push(@{$rv->{'matrix'}}, \@b);
	}

	return $rv;
}


sub readMotifFile {
	my ($file) = @_;
	my @motifs = ();

	open IN, $file;
	my $m='';
	my $count = 0;
	while (<IN>) {
		$count++;
		chomp;
		s/\r//g;
		my @line= split /\t/;
		if ($line[0] =~ s/^>//) {
			push(@motifs, $m) if ($count > 1);
			my @matrix = ();
			$m={name=>'',cons=>'',v=>0,p=>1,logp=>0,str=>'',gapinfo=>0, matrix=>\@matrix,len=>0,mask=>0,fdr=>1};
			$m->{'cons'} = $line[0];
			$m->{'name'} = $line[1] if (@line > 1);
			$m->{'v'} = $line[2] if (@line > 2);
			$m->{'logp'} = $line[3] if (@line > 3);
			$m->{'gapinfo'} = $line[4] if (@line > 4);
			$m->{'str'} = $line[5] if (@line > 5);
			next;
		} else {
			if (@line != 4) {
				print STDERR "Wrong file format\n";
				exit;
			}
			push(@{$m->{'matrix'}}, \@line);
			$m->{'len'}++;
		}
	}
	if (scalar($m->{'matrix'}) > 0) {
		push(@motifs, $m);
	}
	close IN;

	if ($homer2Flag == 0) {
		for (my $i=0;$i<@motifs;$i++) {
			cleanUpMotif($motifs[$i]);
		}
	}

	my @pvalues = ();
	for (my $i=0;$i<@motifs;$i++) {
		push(@pvalues, exp($motifs[$i]->{'logp'}));
	}
	my $fdr = Statistics::benjaminiFDR(\@pvalues);
	for (my $i=0;$i<@motifs;$i++) {
		$motifs[$i]->{'fdr'} = $fdr->[$i];
	}

	return \@motifs;
}


sub cleanUpMotif {
	my ($motif) = @_;
	my $consensus = '';
	for(my $i=0;$i<$motif->{'len'};$i++) {
		my $sum = 0;
		my $code = 'A';
		my $best = 0;
		for (my $j=0;$j<4;$j++) {
			my $v = $motif->{'matrix'}->[$i][$j];
			$sum+= $v;
			if ($v > $best) {
				$best = $v;
				$code = 'A' if ($j==0);
				$code = 'C' if ($j==1);
				$code = 'G' if ($j==2);
				$code = 'T' if ($j==3);
			}
			if ($best < 0.4) {
				$code = 'N';
			}
		}
		for (my $j=0;$j<4;$j++) {
			$motif->{'matrix'}->[$i][$j] /= $sum;
		}
		$consensus .= $code;
	}
	$motif->{'cons'} = $consensus;
	if ($motif->{'name'} eq '') {
		$motif->{'name'} = $consensus;
	}
}
