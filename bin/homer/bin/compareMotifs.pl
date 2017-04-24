#!/usr/bin/env perl
use warnings;
use lib "/home/carlos/tools/homer/.//bin";
my $homeDir = "/home/carlos/tools/homer/./";

my $promoterSeqOffset = -2000;


use POSIX;
use Statistics;

my $motifDir = $homeDir . "/data/knownTFs/";
my $chuckFacts = $homeDir . "/data/misc/ChuckFacts.tsv";

my $fpThresh = 1e-12;
$TopN = 10;
my $styleFlag = "";
$styleFlag = "-S";

#$scoreMethod = 'absDifference';
$scoreMethod = 'correlation';

my $reduceThresh = 0.6;
my $matchThresh = 0.6;
my $minOverlap = 6;
$homer2Flag = 0;
my $forceMatrix = 1;
my $rnaFlag = "";
my $norevopp = 0;
my $reduceMotifFile = "";
my $pvalueThresh = 1e10;
my $foldThresh = -1e10;
my $cmpMatrixFile = "";
my $nofacts = 0;
my $dbview = 0;
my $dbTableFile = '';
my $maxCPUs = 1;
$basicFlag = 0;
$backgroundMinimum = -1;
$minT = 0;
$rehliFlag = 1;
$infoThresh = -1e10;

my $knownFile = $motifDir . "/all/" . "all.motifs";
my $defaultFile = $knownFile;



sub printCMD() {
	print STDERR "\n\tUsage: compareMotifs.pl <motifs file> <output directory> [options]\n";
	print STDERR "\n\tProgram for compares collection of motifs, removing similar ones, and generating HTML output\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-reduceThresh <#> (Similarity Threshold used to remove similar motifs, default: $reduceThresh)\n";
	print STDERR "\t\t-matchThresh <#|T#> (Similarity Threshold to report alignments with known motifs.  Use T#\n";
	print STDERR "\t\t\t[i.e. T10] to match the top # known motifs regardless of similarity, default: T10)\n";
	print STDERR "\t\t-known <known motifs filename> (Library of known motifs to compare to input motifs, default:\n";
	print STDERR "\t\t\t\"$knownFile\"\n";
	print STDERR "\t\t-pvalue <#> (Remove Motifs with a p-value higher than #, default: keep all)\n";
	print STDERR "\t\t-F <#> (Remove Motifs with fold enrichment [target%/bg%] less than <#>, default: keep all)\n";
	print STDERR "\t\t-info <#> (Remove Motifs with information content less than #, default: not used)\n";
	print STDERR "\t\t-b <#> (Remove Motifs with background percentage less than #, default: not used)\n";
	print STDERR "\t\t-minT <#> (Remove Motifs with less than # number of target instances, default: not used)\n";
	print STDERR "\t\t-siteReduce (If homer2 known with -siteReduce was used, motifs will be reduced by\n";
	print STDERR "\t\t\ttheir calculated co-occurence.\n";
	print STDERR "\t\t-stat <freqError|absDifference|correlation> (Stat used to compute matrix similarity.\n";
	print STDERR "\t\t\tdefault: correlation)\n";
	print STDERR "\t\t-bits (scale logos to bit content, default present nucleotide percentage)\n";
	print STDERR "\t\t-rna (output RNA motifs, use RNA motif/miRNA database for comparison)\n";
	print STDERR "\t\t-norevopp (do not check for matches on reverse strand)\n";
	print STDERR "\t\t-reducedMotifs <outputfile> (will remove redundant motifs, output unique ones, and then quit)\n";
	print STDERR "\t\t-matrix <outputfile> (will compute all pairwise similarity scores to matrix, then quit)\n";
	print STDERR "\t\t-nofacts (leave out the humor)\n";
	print STDERR "\t\t-dbview (internal option )\n";
	print STDERR "\t\t-dbTable <table.txt filename> (information to use in dbview, internal option)\n";
	print STDERR "\t\t-cpu <#> (number of CPUs to use)\n";
	print STDERR "\t\t-basic (don't compare to known motifs or print similar motifs)\n";
	print STDERR "\n";
	exit;
}

if (@ARGV < 2) {
	printCMD();
}

my $motifInputFile = $ARGV[0];
my $directory = $ARGV[1];

for (my $i=2;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-known') {
		$knownFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-reduceThresh') {
		$reduceThresh = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-matchThresh') {
		my $option = $ARGV[++$i];
		if ($option =~ s/^T//) {
			$TopN = $option;
		} else {
			$matchThresh = $option;
			$TopN = 0;
		}
	} elsif ($ARGV[$i] eq '-nofacts') {
		$nofacts = 1;
	} elsif ($ARGV[$i] eq '-minT') {
		$minT = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-forceMatrix') {
		$forceMatrix = 1;
	} elsif ($ARGV[$i] eq '-b') {
		$backgroundMinimum = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-basic') {
		$basicFlag = 1;
	} elsif ($ARGV[$i] eq '-matrix') {
		$cmpMatrixFile = $ARGV[++$i];
		print STDERR "\tCreating pairwise matrix file\n";
	} elsif ($ARGV[$i] eq '-bits') {
		$styleFlag = "";
	} elsif ($ARGV[$i] eq '-cpu') {
		$maxCPUs = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-F') {
		$foldThresh = $ARGV[++$i];
		print STDERR "\tFiltering out motifs with fold enrichment less than $foldThresh\n";
	} elsif ($ARGV[$i] eq '-info') {
		$infoThresh = $ARGV[++$i];
		print STDERR "\tFiltering out motifs with fold avg information content less than $infoThresh\n";
	} elsif ($ARGV[$i] eq '-pvalue') {
		$pvalueThresh = $ARGV[++$i];
		print STDERR "\tFiltering out motifs with p-values less than $pvalueThresh\n";
	} elsif ($ARGV[$i] eq '-norevopp') {
		$norevopp = 1;
	} elsif ($ARGV[$i] eq '-dbview') {
		$dbview = 1;
	} elsif ($ARGV[$i] eq '-dbTable') {
		$dbTableFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-reduceMotifs') {
		$reduceMotifFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-rna') {
		$rnaFlag = "-rna";
		if ($knownFile =~ /all\.motifs/) {
			$knownFile = $motifDir . "/" . "all.rna.motifs";
		}
	} elsif ($ARGV[$i] eq '-stat') {
		$scoreMethod = $ARGV[++$i];
	} else {
		print STDERR "!!! Didn't recognize \"$ARGV[$i]\" !!!\n";
		printCMD();
	}
}


@facts = ();
if ($nofacts == 0 && open IN, $chuckFacts) {
	while (<IN>) {
		chomp;
		my @line = split /\t/;
		push(@facts, $line[1]);
	}
	close IN;
}
sub getRandomFact {
	if (@facts < 1) {
		return '';
	}
	my $rand = floor(rand()*scalar(@facts));
	return $facts[$rand];
}

$rand = rand();
my $tmpFile = "$directory/" . $rand . ".motif.tmp";

$height = 2.0;
$widthFactor = 0.75;
$widthOffset = 3.0;

my $mainpageName = "homerResults.html";

my $dbTable = '';
if ($dbTableFile ne '') {
	$dbTable = readDBTable($dbTableFile);
}

my $rawMotifs = readMotifFile($motifInputFile,$pvalueThresh,$foldThresh);
if (@$rawMotifs < 1) {
	print STDERR "\t!!! Filtered out all motifs!!!\n";
	exit;
}

if ($rawMotifs->[0]->{'str'} =~ /^T\:/) {
	$homer2Flag = 1;
	print STDERR "\t(Motifs in homer2 format)\n";
	my $i = 0;
	$percentBackground = 0;
	$numBackground = 0;
	while ($i < @$rawMotifs && ($numBackground < 5 || $percentBackground < 0.01)) {
    	($numTargets,$percentTargets,$numBackground,$percentBackground, $mpvalue,$mfdr) = parseStr2($rawMotifs->[$i]->{'str'});
		if ($mfdr ne 'NA') {
			$homer2Flag = 2; # 2=FDR enabled
		}
    	$percentTargets =~ s/\%//;
    	$percentBackground =~ s/\%//;
		$i++;
	}
    $totalNumTargets = floor($numTargets/($percentTargets*0.01)+0.5);
    $totalNumBackground = floor($numBackground/($percentBackground*0.01)+0.5);
}


if ($dbview == 1) {
	$uniqMotifs = $rawMotifs;
} else {
	print STDERR "\tDetermining similar motifs...";
	if ($homer2Flag && $forceMatrix == 0) {
		$uniqMotifs = reduceMotifsHomer2($rawMotifs);
	} else {
		$uniqMotifs = reduceMotifs($rawMotifs, $reduceThresh,$cmpMatrixFile);
	}
}

if ($reduceMotifFile ne '') {
	printMotifs($uniqMotifs,$reduceMotifFile);
	exit;
}

my $knownMotifs = readMotifFile($knownFile,1e10,-1e10);
printMainPage($uniqMotifs,$knownMotifs, $directory);


#my $n1 = 0;
#my $n2 = 3;
#my ($s,$offset,$dir) = compareMotifs($rawMotifs->[$n1], $rawMotifs->[$n2]);
#printMotif($rawMotifs->[$n1]);
#printMotif($rawMotifs->[$n2]);
#print STDERR "\n\n$rawMotifs->[$n1]->{'name'}\t$rawMotifs->[$n2]->{'name'}\tscore = $s\toffset = $offset\tdir = $dir\n";
`rm -f "$tmpFile"`;
exit;

sub printMainPage {
	my ($motifs, $known, $directory) =@_;
	print STDERR "\tOutputing HTML and sequence logos for motif comparison...\n";
	`mkdir -p "$directory"`;
	`mkdir -p "$directory/homerResults"`;
	open MAIN, ">$directory/$mainpageName";
	if ($dbview == 1) {
		print MAIN "<HTML><HEAD><TITLE>Homer Motifs</TITLE></HEAD>\n";
	} else {
		print MAIN "<HTML><HEAD><TITLE>$directory - Homer de novo Motif Results</TITLE></HEAD>\n";
	}
	print MAIN "<BODY>\n";
	
	if ($dbview == 1) {
		print MAIN "<H1>Homer Motifs</H1>\n";
		print MAIN "<TABLE border=\"1\" cellpading=\"0\" cellspacing=\"0\">\n";
		print MAIN "<TR><TD>Motif</TD><TD>Best Match/Details</TD><TD>Motif File</TD></TR>\n";
	} else {
		print MAIN "<H1>Homer <i>de novo</i> Motif Results ($directory)</H1>\n";
		print MAIN "<A HREF=\"knownResults.html\">Known Motif Enrichment Results</A><BR/>\n";
		print MAIN "<A HREF=\"geneOntology.html\">Gene Ontology Enrichment Results</A><BR/>\n";
		print MAIN "If Homer is having trouble matching a motif to a known motif, try copy/pasting the matrix file into \n";
		print MAIN "<A HREF=\"http://www.benoslab.pitt.edu/stamp/\">STAMP</A><BR/>\n";
		print MAIN "More information on motif finding results: ";
		print MAIN "<A HREF=\"http://biowhat.ucsd.edu/homer/\">HOMER</A>\n";
		print MAIN " | <A HREF=\"http://biowhat.ucsd.edu/homer/motif/index.html\">Description of Results</A>\n";
		print MAIN " | <A HREF=\"http://biowhat.ucsd.edu/homer/motif/practicalTips.html\">Tips</A>\n";
		print MAIN "<BR/>\n";
		if ($homer2Flag) {
			print MAIN "Total target sequences = $totalNumTargets<BR/>\n";
			print MAIN "Total background sequences = $totalNumBackground<BR/>\n";
		}
		print MAIN "<FONT color=\"red\">* - possible false positive</FONT><BR/>\n";
		print MAIN "<TABLE border=\"1\" cellpading=\"0\" cellspacing=\"0\">\n";

		print MAIN "<TR><TD>Rank</TD><TD>Motif</TD><TD>P-value</TD><TD>log P-pvalue</TD>";
		if ($homer2Flag == 2) {
			print MAIN "<TD>q-value/FDR</TD>";
		}
		if ($homer2Flag) {
			print MAIN "<TD>% of Targets</TD><TD>% of Background</TD>\n";
			print MAIN "<TD>STD(Bg STD)</TD>\n";
		}
		print MAIN "<TD>Best Match/Details</TD><TD>Motif File</TD></TR>\n";
	}


	print STDERR "\tChecking de novo motifs against known motifs...\n";
	my $motifNum = 1;
	my $ZZ = scalar(@$motifs);

	my $cpus = 0;
	foreach(@$motifs) {
		my $m = $_;
		my $pid = fork();
		$cpus++;
		if ($pid == 0) {
			my $focusPage = "motif$motifNum.info.html";
			my $similarPage = "motif$motifNum.similar.html";
			my $matrixFile = "motif$motifNum.motif";
			my $logoFile = "motif$motifNum.logo";
			my $mtmpFile = $tmpFile . ".motif" . $motifNum;
			my $pvalue = sprintf("%.3e", exp($m->{'logp'}));
			my $mfdr = '';
			if ($homer2Flag) {
	    		($numTargetsM,$percentTargetsM,$numBackgroundM,$percentBackgroundM, $mpvalue,$mfdr) = parseStr2($m->{'str'});
				$pvalue = $mpvalue;
			}
			my $lp = sprintf("%.3e", $m->{'logp'});
	
			my $bestMatch = printFocusPage($m, $known, $directory . "/homerResults/$focusPage", $motifNum, $directory,
														$dbview, $dbTable);
			if ($dbview && exists($dbTable->{$m->{'name'}})) {
				foreach(@{$dbTable->{$m->{'name'}}}) {
					my $g = $_;
					`cp "$directory/homerResults/$focusPage" "$directory/homerResults/$g.html"`;
				}
			}
				
	
			printMotif($m, $directory . "/homerResults/$matrixFile");
			`profile2seq.pl "$directory/homerResults/$matrixFile" 100 $rnaFlag > "$mtmpFile"`;
			my $width = $widthFactor * ($widthOffset+$m->{'len'});
		
			`seqlogo -a -f "$mtmpFile" -F PNG $styleFlag -c -o "$directory/homerResults/$logoFile" -h $height -w $width`;
			`seqlogo -a -f "$mtmpFile" -F PDF $styleFlag -c -o "$directory/homerResults/$logoFile" -h $height -w $width 2> /dev/null`;
			`rm -f "$mtmpFile"`;
	
			#print STDERR "\t\t$motifNum of $ZZ ($pvalue) similar to $bestMatch\n";
			if ($basicFlag==0 && $dbview == 0) {
				printSimilarPage($m, $directory . "/homerResults/$similarPage", $motifNum, $directory);
			}
			exit(0);
		}
		if ($cpus >= $maxCPUs) {
			my $id = wait();
			$cpus--;
		}
		$motifNum++;
		$numTargetsM = 0;
		$numBackgroundM = 0;
	}
	my $id = 0;
	while ($id >= 0) {	
		$id = wait();
		if ($id == -1) {
		} else {
		}
	}


	print STDERR "\tFormatting HTML page...\n";
	$motifNum = 1;
	foreach(@$motifs) {
		$m = $_;
		my $focusPage = "motif$motifNum.info.html";
		my $similarPage = "motif$motifNum.similar.html";
		my $matrixFile = "motif$motifNum.motif";
		my $logoFile = "motif$motifNum.logo";
		my $pvalue = sprintf("%.3e", exp($m->{'logp'}));
		my $mfdr = '';
		if ($homer2Flag) {
    		($numTargetsM,$percentTargetsM,$numBackgroundM,$percentBackgroundM, $mpvalue,$mfdr) = parseStr2($_->{'str'});
			$pvalue = $mpvalue;
		}
		my $lp = sprintf("%.3e", $m->{'logp'});
		print MAIN "<TR>";
		if ($dbview == 1) {
		} else {
			print MAIN "<TD>$motifNum\n";
			if ($homer2Flag == 2) {
				if ($mfdr > 0.10) {
					print MAIN "<FONT color=\"red\">*</FONT>";
				}
			} elsif ($pvalue > $fpThresh) {
				print MAIN "<FONT color=\"red\">*</FONT>";
			}
		}

		my $bestMatch = "NA";
		open BM, $directory . "/homerResults/$matrixFile";
		while (<BM>) {
			if (/BestGuess:(.*?)\t/) {
				$bestMatch = $1;
				last;
			}
		}
		close BM;
		my $width = $widthFactor * ($widthOffset+$m->{'len'});
	
		print STDERR "\t\t$motifNum of $ZZ ($pvalue) similar to $bestMatch\n";

		print MAIN "</TD><TD><IMG src=\"homerResults/$logoFile.png\"/></TD>";
		if ($dbview == 0) {
			print MAIN "<TD>$pvalue</TD><TD>$lp</TD>";
			if ($homer2Flag == 2) {
				print MAIN "<TD>$mfdr</TD>";
			}
			if ($homer2Flag) {
				print MAIN "<TD>$percentTargetsM</TD><TD>$percentBackgroundM</TD>";
				print MAIN "<TD>$m->{'Tstd'}bp ($m->{'Bstd'}bp)</TD>";
			}
		}
		print MAIN "<TD>$bestMatch<BR/><A target=\"_blank\" HREF=\"homerResults/$focusPage\">More Information</A>";
		if ($dbview && exists($dbTable->{$m->{'name'}})) {
			foreach(@{$dbTable->{$m->{'name'}}}) {
				my $g = $_;
				print MAIN " | <A target=\"blank\" HREF=\"http://www.genecards.org/cgi-bin/carddisp.pl?gene=$g\">$g info</A>";
			}
		}
		if ($basicFlag==0 && $dbview == 0) {
			print MAIN " | <A target=\"blank\" HREF=\"homerResults/$similarPage\">Similar Motifs Found</A>";
		}
		print MAIN "</TD><TD><A target=\"_blank\" HREF=\"homerResults/$matrixFile\">motif file (matrix)</A></TD></TR>\n";
	
		$motifNum++;
		$numTargetsM = 0;
		$numBackgroundM = 0;
	}
	print MAIN "</TABLE>\n";
	print MAIN "<P>" . getRandomFact() . "</P>\n";
	print MAIN "</BODY></HTML>\n";
	close MAIN;
}

sub printSimilarPage {
	my ($motif, $page, $motifNum, $directory) = @_;
	open FOCUS, ">$page";
	print FOCUS "<HTML><HEAD><TITLE>motif$motifNum</TITLE></HEAD>\n";
	print FOCUS "<H2>Information for motif$motifNum</H2>\n";

	my $logoFile = "motif$motifNum.logo.png";
	my $matrixFile = "motif$motifNum.motif";
	print FOCUS "<IMG src=\"$logoFile\"/><BR>\n";
	my $rvMatrixFile = "motif$motifNum" . "RV.motif";
	my $rvLogoFile = "motif$motifNum.rvlogo";
	my $mtmpFile = $tmpFile . ".motif" . $motifNum;
	print FOCUS "Reverse Opposite:<BR/><IMG src=\"$rvLogoFile.png\"/><BR>\n";

	my $lp = sprintf("%.3e", $motif->{'logp'});
	my $pvalue = sprintf("%.3e", exp($lp));
	if ($homer2Flag) {
    	($numTargets,$percentTargets,$numBackground,$percentBackground, $mpvalue,$mfdr) = parseStr2($motif->{'str'});
		$pvalue = $mpvalue;
	} else {
		($totalG, $totalP, $totalM, $totalNP) = parseStr($motif->{'str'});
	}
	my $info = $motif->{'info'};
	print FOCUS "<TABLE border=\"1\" cellpading=\"0\" cellspacing=\"0\">\n";
	print FOCUS "<TR><TD>p-value:</TD><TD>$pvalue</TD></TR>\n";
	print FOCUS "<TR><TD>log p-value:</TD><TD>$lp</TD></TR>\n";
	if ($homer2Flag == 2) {
		print FOCUS "<TR><TD>FDR:</TD><TD>$mfdr</TD></TR>\n";
	}
	print FOCUS "<TR><TD>Information Content per bp:</TD><TD>$info</TD></TR>\n";
	if ($homer2Flag) {
		print FOCUS "<TR><TD>Number of Target Sequences with motif</TD><TD>$numTargets</TD></TR>\n";
		print FOCUS "<TR><TD>Percentage of Target Sequences with motif</TD><TD>$percentTargets</TD></TR>\n";
		print FOCUS "<TR><TD>Number of Background Sequences with motif</TD><TD>$numBackground</TD></TR>\n";
		print FOCUS "<TR><TD>Percentage of Background Sequences with motif</TD><TD>$percentBackground</TD></TR>\n";
		print FOCUS "<TR><TD>Average Position of motif in Targets</TD><TD>$motif->{'Tpos'} +/- $motif->{'Tstd'}bp</TD></TR>\n";
		print FOCUS "<TR><TD>Average Position of motif in Background</TD><TD>$motif->{'Bpos'} +/- $motif->{'Bstd'}bp</TD></TR>\n";
		print FOCUS "<TR><TD>Strand Bias (log2 ratio + to - strand density)</TD><TD>$motif->{'StrandBias'}</TD></TR>\n";
		print FOCUS "<TR><TD>Multiplicity (# of sites on avg that occur together)</TD><TD>$motif->{'Multiplicity'}</TD></TR>\n";
	} else {
		print FOCUS "<TR><TD>Total Number of Sequences:</TD><TD>$totalG</TD></TR>\n";
		print FOCUS "<TR><TD>Total Number of Target Sequences:</TD><TD>$totalP</TD></TR>\n";
		print FOCUS "<TR><TD>Total Instances of Motif:</TD><TD>$totalM</TD></TR>\n";
		print FOCUS "<TR><TD>Total Instances of Motif in Targets:</TD><TD>$totalNP</TD></TR>\n";
	}
	print FOCUS "<TR><TD>Motif File:</TD><TD><A target=\"_blank\" HREF=\"$matrixFile\">file (matrix)</A><BR/>"
						. "<A target=\"_blank\" HREF=\"$rvMatrixFile\">reverse opposite</A></TD></TR>\n";
	print FOCUS "</TABLE>\n";


	print FOCUS "<H3>Similar de novo motifs found</H3>\n";
	print FOCUS "<TABLE border=\"1\" cellpading=\"0\" cellspacing=\"0\">\n";
	print FOCUS "<TR><TD>Rank</TD>";
	if ($homer2Flag && $forceMatrix == 0) {
		print FOCUS "<TD>% of shared sites<BR/>Relative to primary motif</TD><TD>% shared sites<BR/>Relative to redundant motif</TD>";
	} else {
		print FOCUS "<TD>Match Score</TD>";
	}
	print FOCUS "<TD>Redundant Motif</TD><TD>P-value</TD><TD>log P-value</TD>";
	if ($homer2Flag==2) {
		print FOCUS "<TD>FDR</TD>";
	}
	if ($homer2Flag) {
		print FOCUS "<TD>% of Targets</TD><TD>% of Background</TD>";
	}
	print FOCUS "<TD>Motif file</TD></TR>\n";


	my $tmpFile1 = "$directory/" . $rand . ".1.motif.tmp";
	my $tmpFile2 = "$directory/" . $rand . ".2.motif.tmp";
 	my $bestMatch = '';

	my $rank = 0;
	foreach(@{$motif->{'similar'}}) { 
		$rank++;
		my $logp = $_->{'logp'};
		my $pvalue = sprintf("%.3e", exp($logp));
		if ($_->{'homer2'}) {
    		($numTargets,$percentTargets,$numBackground,$percentBackground, $mpvalue,$mfdr) = parseStr2($_->{'str'});
			$pvalue = $mpvalue;
		}
		my $logoFile = "motif$motifNum.similar$rank.logo";
		my $matrixFile = "motif$motifNum.similar$rank.motif";
		printMotif($_, $directory . "/homerResults/$matrixFile");
		`profile2seq.pl "$directory/homerResults/$matrixFile" 100 $rnaFlag > "$mtmpFile"`;
		my $width = $widthFactor * ($_->{'len'}+$widthOffset);
		`seqlogo -a -f "$mtmpFile" -F PNG $styleFlag -c -o "$directory/homerResults/$logoFile" -h $height -w $width`;
		`rm -f "$mtmpFile"`;
	
		print FOCUS "<TR><TD>$rank</TD>";
		if ($homer2Flag && $forceMatrix == 0) {
			my $p1 = $_->{'similar1'};
			my $p2 = $_->{'similar2'};
			print FOCUS "<TD>$p1</TD><TD>$p2</TD>";
		} else {
			my $matchScore = sprintf("%.3f", $_->{'ms'});
			print FOCUS "<TD>$matchScore</TD>";
		}
		print FOCUS "<TD><IMG src=\"$logoFile.png\"/></TD><TD>$pvalue</TD><TD>$logp</TD>";
		if ($homer2Flag==2) {
			print FOCUS "<TD>$mfdr</TD>";
		}
		if ($homer2Flag) {
			print FOCUS "<TD>$_->{'targetP'}</TD><TD>$_->{'backP'}</TD>";
		}
		print FOCUS "<TD><A target=\"_blank\" HREF=\"$matrixFile\">motif file (matrix)</A></TD></TR>\n";
	}
	print FOCUS "</TABLE>\n";
	print FOCUS "<P>" . getRandomFact() . "</P>\n";
	print FOCUS "</BODY></HTML>\n";
	close FOCUS;
}





sub printFocusPage {
	my ($motif, $known, $page,$motifNum,$directory, $dbview, $dbTable) = @_;
	open FOCUS, ">$page";
	print FOCUS "<HTML><HEAD><TITLE>Motif $motifNum</TITLE></HEAD>\n";
	print FOCUS "<H2>Information for $motif->{'name'} (Motif $motifNum)</H2>\n";


	my $logoFile = "motif$motifNum.logo.png";
	my $pdfFile = "motif$motifNum.logo.pdf";
	my $matrixFile = "motif$motifNum.motif";
	my $mtmpFile = $tmpFile . ".motif" . $motifNum;
	print FOCUS "<IMG src=\"$logoFile\"/><BR>\n";

	my $rvMatrixFile = "motif$motifNum" . "RV.motif";
	my $rvMOTIF = revoppMotif($motif);
	printMotif($rvMOTIF,"$directory/homerResults/$rvMatrixFile");
	#`revoppMotif.pl "$directory/homerResults/$matrixFile" > "$directory/homerResults/$rvMatrixFile"`;

	my $rvLogoFile = "motif$motifNum.rvlogo";
	my $rvPdfFile = "motif$motifNum.rvlogo.pdf";
	`profile2seq.pl "$directory/homerResults/$rvMatrixFile" 100 $rnaFlag > "$mtmpFile"`;
	my $width = $widthFactor * ($widthOffset+$motif->{'len'});
	`seqlogo -a -f "$mtmpFile" -F PNG $styleFlag -c -o "$directory/homerResults/$rvLogoFile" -h $height -w $width`;
	`seqlogo -a -f "$mtmpFile" -F PDF $styleFlag -c -o "$directory/homerResults/$rvLogoFile" -h $height -w $width 2> /dev/null`;
	`rm "$mtmpFile"`;
	print FOCUS "Reverse Opposite:<BR/><IMG src=\"$rvLogoFile.png\"/><BR>\n";

	my $lp = sprintf("%.3e", $motif->{'logp'});
	my $pvalue = sprintf("%.3e", exp($lp));
	if ($homer2Flag) {
    	($numTargets,$percentTargets,$numBackground,$percentBackground, $mpvalue,$mfdr) = parseStr2($motif->{'str'});
		$pvalue = $mpvalue;
	} else {
		($totalG, $totalP, $totalM, $totalNP) = parseStr($motif->{'str'});
	}
	my $info = $motif->{'info'};
	print FOCUS "<TABLE border=\"1\" cellpading=\"0\" cellspacing=\"0\">\n";
	print FOCUS "<TR><TD>p-value:</TD><TD>$pvalue</TD></TR>\n";
	print FOCUS "<TR><TD>log p-value:</TD><TD>$lp</TD></TR>\n";
	if ($homer2Flag == 2) {
		print FOCUS "<TR><TD>FDR:</TD><TD>$mfdr</TD></TR>\n";
	}
	print FOCUS "<TR><TD>Information Content per bp:</TD><TD>$info</TD></TR>\n";
	if ($homer2Flag) {
		print FOCUS "<TR><TD>Number of Target Sequences with motif</TD><TD>$numTargets</TD></TR>\n";
		print FOCUS "<TR><TD>Percentage of Target Sequences with motif</TD><TD>$percentTargets</TD></TR>\n";
		print FOCUS "<TR><TD>Number of Background Sequences with motif</TD><TD>$numBackground</TD></TR>\n";
		print FOCUS "<TR><TD>Percentage of Background Sequences with motif</TD><TD>$percentBackground</TD></TR>\n";
		print FOCUS "<TR><TD>Average Position of motif in Targets</TD><TD>$motif->{'Tpos'} +/- $motif->{'Tstd'}bp</TD></TR>\n";
		print FOCUS "<TR><TD>Average Position of motif in Background</TD><TD>$motif->{'Bpos'} +/- $motif->{'Bstd'}bp</TD></TR>\n";
		print FOCUS "<TR><TD>Strand Bias (log2 ratio + to - strand density)</TD><TD>$motif->{'StrandBias'}</TD></TR>\n";
		print FOCUS "<TR><TD>Multiplicity (# of sites on avg that occur together)</TD><TD>$motif->{'Multiplicity'}</TD></TR>\n";
	} else {
		print FOCUS "<TR><TD>Total Number of Sequences:</TD><TD>$totalG</TD></TR>\n";
		print FOCUS "<TR><TD>Total Number of Target Sequences:</TD><TD>$totalP</TD></TR>\n";
		print FOCUS "<TR><TD>Total Instances of Motif:</TD><TD>$totalM</TD></TR>\n";
		print FOCUS "<TR><TD>Total Instances of Motif in Targets:</TD><TD>$totalNP</TD></TR>\n";
	}
	print FOCUS "<TR><TD>Motif File:</TD><TD><A target=\"_blank\" HREF=\"$matrixFile\">file (matrix)</A><BR/>"
						. "<A target=\"_blank\" HREF=\"$rvMatrixFile\">reverse opposite</A></TD></TR>\n";
	print FOCUS "<TR><TD>PDF Format Logos:</TD><TD><A HREF=\"$pdfFile\">forward logo</A><BR/>"
						. "<A HREF=\"$rvPdfFile\">reverse opposite</A></TD></TR>\n";
	print FOCUS "</TABLE>\n";

 	my $bestMatch = '';

	if ($basicFlag == 0) {

		my $matches = matchMotif($motif, $known);

		print FOCUS "<H3>Matches to Known Motifs</H3>\n";
		print FOCUS "<TABLE border=\"1\" cellpading=\"0\" cellspacing=\"0\">\n";
	
	
		my $tmpFile1 = "$directory/" . $rand . ".1.motif" . $motifNum . ".tmp";
		my $tmpFile2 = "$directory/" . $rand . ".2.motif" . $motifNum . ".tmp";
	
		for (my $i=0;$i<@$matches;$i++) { 
			if ($TopN > 0) {
				last if ($i >= $TopN);
			} else {
				last if ($matches->[$i]->{'score'} < $matchThresh);
			}
			my $rank = $i+1;
			if ($i==0) {
				$bestMatch = $matches->[$i]->{'motif'}->{'name'};
			}
			print FOCUS "<TR><TD>\n";
			print FOCUS "<H4>$matches->[$i]->{'motif'}->{'name'}</H4>\n";
			print FOCUS "<TABLE><TR><TD>\n";
			print FOCUS "<TABLE>\n";
			print FOCUS "<TR><TD>Match Rank:</TD><TD>$rank</TD></TR>\n";

			my $score = sprintf("%.2f", $matches->[$i]->{'score'});
			print FOCUS "<TR><TD>Score:</TD><TD>$score</TD</TR>\n";
			print FOCUS "<TR><TD>Offset:</TD><TD>$matches->[$i]->{'offset'}</TD</TR>\n";
			my $orientation = "forward strand";
			my $offset = $matches->[$i]->{'offset'};
			my $visMotif = $matches->[$i]->{'motif'};
	
	
			if ($matches->[$i]->{'dir'} == 1) {
				$orientation = "reverse strand" ;
				$visMotif = revoppMotif($visMotif);
				cleanUpMotif($visMotif);
			}
			my $a1 = '';
			my $a2 = '';
			open M1, ">$tmpFile1";
			open M2, ">$tmpFile2";
			print M1 ">M1\n";
			print M2 ">M2\n";
			if ($offset < 0) {
				for (my $j=0; $j>$offset;$j--) {
					$a1 .= '-';
					print M1 "0.25\t0.25\t0.25\t0.25\n";
				}
			} else {
				for (my $j=0; $j<$offset;$j++) {
					$a2 .= '-';
					print M2 "0.25\t0.25\t0.25\t0.25\n";
				}
			}
			$a1 .= $motif->{'cons'};
			foreach(@{$motif->{'matrix'}}) {
				my $cc = 0;
				foreach(@$_) {
					print M1 "\t" if ($cc>0);
					$cc++;
					print M1 $_;
				}
				print M1 "\n";
			}
			$a2 .= $visMotif->{'cons'};
			foreach(@{$visMotif->{'matrix'}}) {
				my $cc = 0;
				foreach(@$_) {
					print M2 "\t" if ($cc>0);
					$cc++;
					print M2 $_;
				}
				print M2 "\n";
			}
			while (length($a1) > length($a2)) {
				$a2 .= '-';
				print M2 "0.25\t0.25\t0.25\t0.25\n";
			}
			while (length($a1) < length($a2)) {
				$a1 .= '-';
				print M1 "0.25\t0.25\t0.25\t0.25\n";
			}
			close M1;
			close M2;
			my $img1 = "motif$motifNum-$rank-homer";
			my $img2 = "motif$motifNum-$rank-known";
			`profile2seq.pl "$tmpFile1" 100 $rnaFlag > "$mtmpFile"`;
			$width = $widthFactor * ($widthOffset+length($a1));
			`seqlogo -a -f "$mtmpFile" -F PNG $styleFlag -c -o "$directory/homerResults/$img1" -h $height -w $width`;
			`profile2seq.pl "$tmpFile2" 100 $rnaFlag > "$mtmpFile"`;
			$width = $widthFactor * ($widthOffset+length($a2));
			`seqlogo -a -f "$mtmpFile" -F PNG $styleFlag -c -o "$directory/homerResults/$img2" -h $height -w $width`;
			`rm -f "$mtmpFile"`;
		
			print FOCUS "<TR><TD>Orientation:</TD><TD>$orientation</TD></TR>\n";
			print FOCUS "<TR><TD>Alignment:</TD><TD>";
			print FOCUS "<FONT face=\"courier\">$a1</FONT><BR>";
			print FOCUS "<FONT face=\"courier\">$a2</FONT>";
			print FOCUS "</TD></TR></TABLE>\n";
			print FOCUS "</TD><TD>\n";
			print FOCUS "<IMG src=\"$img1.png\"/><BR/><IMG src=\"$img2.png\"/>\n";
			print FOCUS "</TD></TR></TABLE>\n";
			print FOCUS "</TD></TR>\n";
			
			`rm -f "$tmpFile1" "$tmpFile2"`;
		}
		print FOCUS "</TABLE>\n";

	}
	print FOCUS "<P>" . getRandomFact() . "</P>\n";
	print FOCUS "</BODY></HTML>\n";
	close FOCUS;
	return $bestMatch;
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
sub parseStr2 {
    my ($str) = @_;
	$str =~ s/\)//g;
    my @str = split(/\,/,$str);
    my @t = split(/[\:\(]/,$str[0]);
    my $numTargets = $t[1];
    my $percentTargets = $t[2];
    my @b = split(/[\:\(]/,$str[1]);
    my $numBackground = $b[1];
    my $percentBackground = $b[2];
    my @c = split(/[\:]/,$str[2]);
    my $pvalue = $c[1];
	my $fdr = 'NA';
	if (@str > 3) {
    	my @c = split(/[\:]/,$str[3]);
    	$fdr = sprintf("%.3lf",$c[1]);
	}
    return ($numTargets, $percentTargets, $numBackground, $percentBackground, $pvalue,$fdr);
}



sub matchMotif {
	my ($motif, $known) = @_;

	my @results = ();
	for (my $i=0;$i<@$known;$i++) {
		my ($score, $offset, $dir, $L) = (0,0,0,0);
		if ($rehliFlag) {
			($score, $offset, $dir, $L) = compareMotifsRehli($motif, $known->[$i]);
		} else {
			($score, $offset, $dir, $L) = compareMotifs($motif, $known->[$i]);
		}
		if ($scoreMethod ne 'correlation') {
			$score = $score/$motif->{'len'};
		}
		my $res = {motif=>$known->[$i], score=>$score, offset=>$offset, dir=>$dir};
		push(@results, $res);
	}
	@results = sort {$b->{'score'} <=> $a->{'score'}} @results;
	if (@results > 0) {
		$motif->{'bestGuess'} = $results[0]->{'motif'}->{'name'} . "(" . sprintf("%.3f",$results[0]->{'score'}) . ")";
	}

	return \@results;
}

sub reduceMotifsHomer2 {
	my ($motifs) = @_;

	my %list = ();
	foreach(@$motifs) {
		$list{$_->{'name'}} = $_;
	}
	my @new = ();
	foreach(@$motifs) {
		if ($_->{'similarName'} ne '') {
			my $m = $list{$_->{'similarName'}};
			push(@{$m->{'similar'}}, $_);
		} else {
			push(@new, $_);
		}
	}
	return \@new;
}

sub reduceMotifs {
	my ($motifs, $threshold, $cmpMatrixFile) = @_;
	my %mask = ();	
	my @m = sort {$a->{'logp'} <=> $b->{'logp'}} @$motifs;
	my $startNum = scalar(@m);
	my @newList = ();
	my $lastNUM = -1;

	my @cmpMatrix = ();
	if ($cmpMatrixFile ne '') {
		for (my $i=0;$i<@m;$i++) {
			my @a = ();
			for (my $j=0;$j<@m;$j++) {
				if ($i==$j) {
					push(@a, 1);
				} else {
					push(@a, 0);
				}
			}
			push(@cmpMatrix, \@a);
		}
	}

	for (my $i=0;$i<@m;$i++) {
		if ($i==$lastNUM) {
			print STDERR "BLAH = $i\n";
		}
		$lastNUM = $i;
		next if ($m[$i]->{'mask'} ==1 && $cmpMatrixFile eq '');
		#print STDERR "$m[$i]->{'name'} $i\n";
		push(@newList, $m[$i]);
		for (my $j=$i+1;$j<@m;$j++) {
			next if ($m[$j]->{'mask'} ==1 && $cmpMatrixFile eq '');
			my ($s,$o,$d,$L) = (0,0,0,0);
			if ($rehliFlag) {
				($s,$o,$d,$L) = compareMotifsRehli($m[$i], $m[$j]);
			} else {
				($s,$o,$d,$L) = compareMotifs($m[$i], $m[$j]);
			}
			if ($scoreMethod ne 'correlation') {
				$s = $s/$L;
			}
			if ($cmpMatrixFile ne '') {
				$cmpMatrix[$i][$j] = $s;
				$cmpMatrix[$j][$i] = $s;
			}
			if ($s > $threshold) {
				$m[$j]->{'mask'} = 1;
				$m[$j]->{'ms'} = $s;
				push(@{$m[$i]->{'similar'}}, $m[$j]);
			}
		}
	}

	if ($cmpMatrixFile ne '') {
		open OUT, ">$cmpMatrixFile" or die "Couldn't open matrix output file $cmpMatrixFile\n";
		print OUT "Motif Comparison";
		for (my $i=0;$i<@m;$i++) {
			my $n = $m[$i]->{'name'};
			print OUT "\t$n";
		}
		print OUT "\n";
		for (my $i=0;$i<@m;$i++) {
			my $n = $m[$i]->{'name'};
			print OUT "$n";
			for (my $j=0;$j<@m;$j++) {
				print OUT "\t$cmpMatrix[$i][$j]";
			}
			print OUT "\n";
		}
		close OUT;
		print STDERR "\n\n";
		exit;
	}
		


	my $endNum = scalar(@newList);
	print STDERR " $startNum reduced to $endNum motifs\n";
	return \@newList;
}

sub printMotif {
	my ($motif, $file) = @_;

	my $name = $motif->{'name'};
	if (defined($motif->{'bestGuess'})) {
		$name .= ",BestGuess:" . $motif->{'bestGuess'};
	}

	open MOTIF, ">$file";
	print MOTIF ">$motif->{'cons'}\t$name\t$motif->{'v'}\t$motif->{'logp'}\t$motif->{'gapinfo'}\t$motif->{'str'}\n";
	my $matrix = $motif->{'matrix'};
	foreach(@{$matrix}) {
		for (my $i=0;$i<4;$i++) {
			print MOTIF "\t" if ($i>0);
			print MOTIF "$_->[$i]";
		}
		print MOTIF "\n";
	}
	close MOTIF;
}
sub printMotifs {
	my ($motifs, $file) = @_;
	open MOTIF, ">$file";
	foreach(@$motifs) {
		my $motif = $_;
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
	close MOTIF;
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
	my $bestLength = 0;
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
		my $curLength = @mm1;

		if (@mm1 < $minOverlap) {
			next;
		}
		if (@mm2 < $minOverlap) {
			next;
		}

		my $score = scoreComparison(\@mm1, \@mm2);
		my $rvScore = scoreComparison(\@mm1, \@mm2r);
#print STDERR "$offset1 ($len1)\t$offset2 ($len2)\t$officialOffset\t$score\t$rvScore\n";
		if ($norevopp == 0 && $rvScore > $score) {
			if ($rvScore > $bestScore) {
				$bestScore = $rvScore;
				$bestOffset = $officialOffset;
				$bestDirection = 1;
				$bestLength = $curLength;
			}
		} else {
			if ($score > $bestScore) {
				$bestScore = $score;
				$bestOffset = $officialOffset;
				$bestDirection = 0;
				$bestLength = $curLength;
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
	return ($bestScore, $bestOffset, $bestDirection,$bestLength);
}

sub compareMotifsRehli  {
	my ($m1, $m2) = @_;

	my $bestOffset = 0;
	my $bestDirection = 0;
	my $bestScore = -1e10;

	my $rv2 = revoppMotif($m2);

	my @default = (0.25,0.25,0.25,0.25);
	my $len1 = $m1->{'len'};
	my $len2 = $m2->{'len'};
	my $offset2 = $len2-1;
	my $offset1 = 0;
	my $bestLength = 0;
	for (my $offset2 = $len1-1;$offset2>-1*$len2;$offset2--) {
		my $officialOffset = -1*$offset2+$offset1;

		my $max1 = $len1 - $offset1;
		my $max2 = $len2 - $offset2;
		my $curLen = $max1;
		if ($max2 < $curLen) {
			$curLen = $max2;
		}
		my @mm1 = ();
		for (my $i=$offset2;$i<0;$i++) {
			push(@mm1, \@default);
		}
		for (my $i=0;$i<$len1;$i++) {
			push(@mm1, $m1->{'matrix'}->[$i]);
		}
		for (my $i=0;$i<$offset2+$len2-$len1;$i++) {
			push(@mm1, \@default);
		}
		my @mm2 = ();
		my @mm2r = ();
		for (my $i=0;$i<$offset2;$i++) {
			push(@mm2, \@default);
			push(@mm2r, \@default);
		}
		for (my $i=0;$i<$len2;$i++) {
			push(@mm2, $m2->{'matrix'}->[$i]);
			push(@mm2r, $rv2->{'matrix'}->[$i]);
		}
		while (scalar(@mm2) < scalar(@mm1)) {
			push(@mm2, \@default);
			push(@mm2r, \@default);
		}
		my $curLength = @mm1;


		my $n2 = scalar(@mm2);
		my $n1 = scalar(@mm1);
		#print STDERR "$offset2  $len1 $len2 $n1 $n2\n";
		#printMatrix(\@mm1);
		#printMatrix(\@mm2);

		my $score = scoreComparison(\@mm1, \@mm2);
		my $rvScore = scoreComparison(\@mm1, \@mm2r);
		#print STDERR "$offset2 ($len1)\t$offset2 ($len2)\t$officialOffset\t$score\t$rvScore\n";
		if ($norevopp == 0 && $rvScore > $score) {
			if ($rvScore > $bestScore) {
				$bestScore = $rvScore;
				$bestOffset = $officialOffset;
				$bestDirection = 1;
				$bestLength = $curLength;
			}
		} else {
			if ($score > $bestScore) {
				$bestScore = $score;
				$bestOffset = $officialOffset;
				$bestDirection = 0;
				$bestLength = $curLength;
			}
		}

	}
	$bestOffset *= -1;
	#print STDERR "$bestScore\t$bestOffset\t$bestDirection\t$bestLength\n";
	return ($bestScore, $bestOffset, $bestDirection,$bestLength);
}

sub printMatrix {
	my ($m) = @_;
	print STDERR "Matrix:\n";
	foreach(@$m) {
		foreach(@$_) {
			print STDERR "\t$_";
		}
		print STDERR "\n";
	} 
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
	} elsif ($scoreMethod eq 'freqError') {
		for (my $i=0;$i<$len;$i++) {
			my $curScore = 0;
			my $expectedScore = 0;
			for (my $j=0;$j<4;$j++) { 
				my $diff = $m1->[$i][$j]-$m2->[$i][$j];
				$curScore -= $diff*$diff;
				for (my $k=0;$k<4;$k++) { 
					my $diff2 = $m1->[$i][$j]-$m2->[$i][$k];
					$expectedScore -= $diff2*$diff2/4;
				}
			}
			$score +=  -1*($curScore-$expectedScore)/$expectedScore;
		}
	} elsif ($scoreMethod eq 'correlation') {
		my @a1 = ();
		my @a2 = ();
		foreach(@$m1) {
			foreach(@$_) {
				push(@a1,$_);
			}
		}
		foreach(@$m2) {
			foreach(@$_) {
				push(@a2,$_);
			}
		}
		my $lp = 0;
		($score,$lp) = Statistics::correlation(\@a1,\@a2);
	}
	#print STDERR "$score\n";
	return $score;
}

sub revoppMotif {
	my ($m) = @_;
	my @a = ();
	my $rv = {matrix=>\@a,name=>$m->{'name'},cons=>$m->{'cons'},
					v=>$m->{'v'},p=>$m->{'p'},logp=>$m->{'logp'},str=>$m->{'str'},
					gapinfo=>$m->{'gapinfo'}, len=>$m->{'len'},mask=>$m->{'mask'},
					homer2=>$m->{'homer2'}};
	for (my $i=-1+$m->{'len'};$i>=0;$i--){ 
		my @b = ();
		for (my $j=3;$j>=0;$j--) {
			push(@b, $m->{'matrix'}->[$i][$j]);
		}
		push(@{$rv->{'matrix'}}, \@b);
	}
	revoppConsensus($rv);

	return $rv;
}
sub revoppConsensus {
	my ($motif) = @_;
	my $new = reverse($motif->{'cons'});
	$new =~ s/A/X/g;
	$new =~ s/T/A/g;
	$new =~ s/X/T/g;

	$new =~ s/C/X/g;
	$new =~ s/G/C/g;
	$new =~ s/X/G/g;

	$new =~ s/R/X/g;
	$new =~ s/Y/R/g;
	$new =~ s/X/Y/g;

	$new =~ s/M/X/g;
	$new =~ s/K/M/g;
	$new =~ s/X/K/g;

	$new =~ s/B/X/g;
	$new =~ s/V/B/g;
	$new =~ s/X/V/g;

	$new =~ s/D/X/g;
	$new =~ s/H/D/g;
	$new =~ s/X/H/g;
	$motif->{'cons'} = $new;
}


sub readMotifFile {
	my ($file,$pvalueThresh,$foldThresh) = @_;
	my @motifs = ();


	my $logpThresh = log($pvalueThresh);
	open IN, $file or die "Could not open file: \"$file\"\n";
	my $m='';
	my $count = 0;
	my $badFlag = 0;
	while (<IN>) {
		$count++;
		chomp;
		s/\r//g;
		my @line= split /\t/;
		if ($line[0] =~ s/^>//) {
			if ($count > 1 && $m->{'logp'} < $logpThresh && $badFlag == 0
					&& $m->{'fold'} > $foldThresh) {
				
				calcInfoContent($m);
				if ($m->{'info'} > $infoThresh) {
					push(@motifs, $m) 
				}
			}
			$badFlag= 0;
			my @matrix = ();
			my @similar = ();
			$m={name=>'',cons=>'',v=>0,p=>1,logp=>0,str=>'',gapinfo=>0, matrix=>\@matrix,len=>0,mask=>0,
					similar=>\@similar,
					homer2=>0,
					targetN=>0,
					targetP=>0,
					backN=>0,
					backP=>0,
					Tpos=>0,
					Tstd=>0,
					Bpos=>0,
					Bstd=>0,
					StrandBias=>0,
					Multiplicity=>1,
					similar1=>0,
					similar2=>0,
					similarName=>"",
					logp=>1,
					fold=>10000,
					FDR=>'NA'
				};
			$m->{'cons'} = $line[0];
			$m->{'name'} = $line[1] if (@line > 1);
			$m->{'v'} = $line[2] if (@line > 2);
			$m->{'logp'} = $line[3] if (@line > 3);
			$m->{'gapinfo'} = $line[4] if (@line > 4);
			if (@line > 5) {
				$m->{'str'} = $line[5];
				if ($m->{'str'} =~ /^T\:/) {
					$m->{'homer2'} = 1;
					#$homer2Flag = 1;
    				my ($numTargets, $percentTargets, $numBackground, $percentBackground, $pvalue,$fdr) = parseStr2($line[5]);
					$m->{'targetN'} = $numTargets;
					$m->{'targetP'} = $percentTargets;
					$m->{'backN'} = $numBackground;
					$m->{'backP'} = $percentBackground;
					$m->{'p'} = $pvalue;
					$m->{'FDR'} = $fdr;
					my $tp = $m->{'targetP'};
					my $tb = $m->{'backP'};
					$tp =~ s/\%//;
					$tb =~ s/\%//;
					if ($tb > 0) {
						$m->{'fold'} = $tp/$tb;
					}
					if ($tb < $backgroundMinimum*100) {
						$badFlag = 1;
					}
					if ($numTargets < $minT) {
						$badFlag = 1;
					}
				} else {
					$m->{'homer2'} = 0;
				}
			}
			if (@line > 6) {
				my @stats = split /\,/, $line[6];
				for (my $i=0;$i<@stats;$i++) {
					my @a = split /\:/, $stats[$i];
					$m->{$a[0]} = $a[1];
				}
			}
			if (@line > 7) {
				$line[7] =~ /Shares (.+?) of sites with (.*?)\((.+?)\)/;
				$m->{'similar1'} = $1;
				$m->{'similarName'} = $2;
				$m->{'similar2'} = $3;
			}
			next;
		} else {
			if (@line != 4) {
				print STDERR "Wrong file format\n";
				exit;
			}
			push(@{$m->{'matrix'}}, \@line);
			$m->{'len'}++;
			for (my $k=0;$k<@line;$k++) {
				$badFlag = 1 if ($line[$k] eq '-nan');
			}
		}
	}
	if (scalar($m->{'matrix'}) > 0) {
		if ($m->{'logp'} < $logpThresh && $badFlag == 0
					&& $m->{'fold'} > $foldThresh) {
			calcInfoContent($m);
			if ($m->{'info'} > $infoThresh) {
				push(@motifs, $m);
			}
		}
	}
	close IN;

	for (my $i=0;$i<@motifs;$i++) {
		cleanUpMotif($motifs[$i]);
	}
	return \@motifs;
}
sub calcInfoContent {
	my ($motif) = @_;
	my $info = 0;
	my $L = @{$motif->{'matrix'}};
	for (my $i=0;$i<$L;$i++) {
		$info += 2;
		for (my $j=0;$j<4;$j++) {
			my $v = $motif->{'matrix'}->[$i]->[$j];
			$v = 0.0001 if ($v < 0.0001);
			$info += $v*log($v)/2;
		}	
	}
	$info /= $L;
	$motif->{'info'} = sprintf("%.3f",$info);
}

sub cleanUpMotif {
	my ($motif) = @_;
	return if ($motif->{'homer2'});
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
			if ($sum < 1e-100) {
				$motif->{'matrix'}->[$i][$j] = 0.25;
			} else {
				$motif->{'matrix'}->[$i][$j] /= $sum;
			}
		}
		$consensus .= $code;
	}
	$motif->{'cons'} = $consensus;
	if ($motif->{'name'} eq '') {
		$motif->{'name'} = $consensus;
	}
}

sub readDBTable {
	my ($filename) = @_;
	open IN, $filename;
	my %table = ();
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		my $name = $line[1];
		my $gene = "";
		if (@line > 11) {
			$gene = $line[11];
		}
		if (!exists($table{$name})) {
			my @a = ();
			$table{$name} = \@a;
		}
		push(@{$table{$name}},$gene)
		#print STDERR "$name=$gene\n";
	}
	close IN;
	return \%table;
}
