#!/usr/bin/env perl
use warnings;
use lib "/home/carlos/tools/homer/.//bin";
my $homeDir = "/home/carlos/tools/homer/./";



use POSIX;
use HomerConfig;
use Statistics;


my $config = HomerConfig::loadConfigFile();

sub printCMD {
	print STDERR "\n\tUsage: analyzeChIP-Seq.pl <exp tag directory> <genome> [global options] [specific options]\n";
	print STDERR "\tIf alignment files have not been parsed yet use:\n";
	print STDERR "\t       analyzeChIP-Seq.pl <exp directory> <genome> [global options] -A s_1_eland_result.txt ...\n";
	print STDERR "\n\tAutomates analysis for a single experiment (creates an index.html file):\n";
	print STDERR "\t\t(A) Alignment Processing..................makeTagDirectory\n";
	print STDERR "\t\t(B) Peak Finding / UCSC Visualization.....findPeaks\n";
	print STDERR "\t\t(C) Motif Finding.........................findMotifsGenome.pl\n";
	print STDERR "\t\t(D) Peak Annotation.......................annotatePeaks.pl\n";

	print STDERR "\n\tAvailable Genomes (required argument): (name,org,directory,default promoter set)\n";
	foreach(keys %{$config->{'GENOMES'}}) {
		print STDERR "\t\t$_\t$config->{'GENOMES'}->{$_}->{'org'}\t$config->{'GENOMES'}->{$_}->{'directory'}"
							. "\t$config->{'GENOMES'}->{$_}->{'promoters'}\n";
	}

	print STDERR "\n\tGlobal Options:\n";
	print STDERR "\t\t-i <Input/IgG experiment tag directory> (or)\n";
	print STDERR "\t\t\t-iraw <input alignment file> [input alignment file2] ...(creates temporary input directory)\n";
	print STDERR "\t\t-style <factor|histone|...> (findPeaks style for peak finding, default: factor)\n";
	print STDERR "\t\t\t-msize (size of reagion to perform motif finding on, default: factor[200], histone[1000])\n";
	print STDERR "\t\t\t-focus (2ndary motif analysis on \"focused\" TF peaks, using 50bp regions, or\n";
	print STDERR "\t\t\t\tanalysis of NFR regions of \"histone\" peaks, using 200bp regions)\n";
	print STDERR "\t\t-p <#> (number of CPUs to run motif finding with, default: 1)\n";
	print STDERR "\t\t-enhancer (when performing analysis, limit motif finding to peaks >3kb from TSS)\n";
	print STDERR "\t\t-force (forces all steps)\n";
	print STDERR "\t\t-mask (Motif finding with repeat masked genome, or add \"r\" to end of genome name)\n";
	print STDERR "\t\t-skipFreq (skips nucleotide frequency and GC quality control calculations)\n";
	print STDERR "\t\t-cpg (For motif finding, use CpG% sequence bias correction, default: GC%)\n";
	print STDERR "\t\t-tagGO (perform Genome Ontology Analysis on tags - need ~ 3Gbs of memory)\n";

	print STDERR "\n\tThis program will attempt to detect previous analysis\n";
	print STDERR "\t\tTo skip analysis:  -s <A|B|C|D> (i.e. \"-s D\")\n";
	print STDERR "\t\tTo force(overwrite) analysis: -f <A|B|C|D> (i.e. \"-f C\" or \"-f B C D\")\n";

	print STDERR "\n\tProgram Specific Options (will be passed to individual programs):\n";
	print STDERR "\tUse will override default options!!\n";
	print STDERR "\t\t-A <makeTagDirectory specific options> (normally -A experiment1.sam  ...)\n";
	print STDERR "\t\t-B <findPeaks specific options>\n";
	print STDERR "\t\t-C <findMotifsGenome.pl specific options>\n";
	print STDERR "\t\t-D <annotatePeaks.pl specific options>\n";


	print STDERR "\n\tDefault Options: (INPUT_DIRECTORY, used if provided as global options)\n";
	print STDERR "\t\t-A\n";
	print STDERR "\t\t-B -i INPUT_DIRECTORY -style factor\n";
	print STDERR "\t\t-C -len 8,10,12 -S 25 -size 200\n";
	print STDERR "\t\t-D -d EXP_DIRECTORY INPUT_DIRECTORY -go EXP_DIRECTORY/GOanalysis\n";

	print STDERR "\n";
	exit;
}

if (@ARGV < 2) {
	printCMD();
}

my $size = '';
my $skipFreq = 0;
my $inputDir = '';
my $cpgFlag = 0;
my $focusFlag = 0;
my $skipGO = 1;
my $distalFlag = 0;
my $enhancerFlag = 0;
my $expDir = $ARGV[0];
my $style = 'factor';
my $msize = "";
my $numCPUs = 1;
my @irawFiles = ();


my $genome = $ARGV[1];
my $maskFlag = "";
if ($genome =~ s/r$//) {
	$maskFlag = " -mask ";
}
my $genomeDir = "";
my $genomeParseDir = "";
my $customGenome = 0;
if (!exists($config->{'GENOMES'}->{$genome})) {
	$customGenome = 1;
	($genome,$genomeDir,$genomeParseDir) = HomerConfig::parseCustomGenome($genome);
} else {
	$genomeDir = $genome;
}

my $Aoptions = "";
my $Boptions = "";
my $Coptions = "";
my $Doptions = "";

my %force = ();
my %skip = ();

for (my $i=2;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-size') {
		print STDERR "\tThe option \"-size #\" is not used anymore.\n";
		#$size = $ARGV[++$i];
		#print STDERR "\tPeak size set to: $size\n";
		next;
	} elsif ($ARGV[$i] eq '-msize') {
		$msize = $ARGV[++$i];
		print STDERR "\tWill analyze peaks for motif enrichment with a size of $msize bp\n";
	} elsif ($ARGV[$i] eq '-mask') {
		$maskFlag = " -mask ";
	} elsif ($ARGV[$i] eq '-p') {
		$numCPUs = $ARGV[++$i];
		print STDERR "\tWill run motif finding with $numCPUs CPUs\n";
	} elsif ($ARGV[$i] eq '-cpg') {
		$cpgFlag = 1;
		print STDERR "\tCpG bias will be corrected instead of GC\n";
	} elsif ($ARGV[$i] eq '-tagGO') {
		$skipGO = 0;
		print STDERR "\tPerforming Genome Ontology analysis on tags\n";
		next;
	} elsif ($ARGV[$i] eq '-style') {
		$style = $ARGV[++$i];
		print STDERR "\tFinding peaks with style \"$style\"\n";
		
	} elsif ($ARGV[$i] eq '-skipFreq') {
		$skipFreq = 1;
		print STDERR "\tSkipping frequencing calculations\n";
		next;
	} elsif ($ARGV[$i] eq '-focus') {
		print STDERR "\tWill perform analysis on 90% focused peaks\n";
		$focusFlag = 1;
		next;
	} elsif ($ARGV[$i] eq '-iraw') {
        print STDERR "\tInput alignment files:\n";
        my $bail = 0;
        while ($ARGV[++$i] !~ /^\-/) {
            push(@irawFiles, $ARGV[$i]);
            print STDERR "\t\t$ARGV[$i]\n";
            if ($i>=@ARGV-1) {
                $bail = 1;
                last;
            }
        }
        last if ($bail == 1);
        $i--;
	} elsif ($ARGV[$i] eq '-enhancer') {
		$size = 1000;
		$Coptions .= "";
		$distalFlag = 1;
		$enhancerFlag = 1;
		print STDERR "\tAdjusting parameters for histone mark analysis\n";
		next;
	} elsif ($ARGV[$i] eq '-i') {
		$inputDir = $ARGV[++$i];
		print STDERR "\tInput Directory: $inputDir\n";
		next;
	} elsif ($ARGV[$i] eq '-force') {
		$force{'A'} = 1;
		$force{'B'} = 1;
		$force{'C'} = 1;
		$force{'D'} = 1;
		next;
	} elsif ($ARGV[$i] eq '-f' || $ARGV[$i] eq '-force') {
		for ($i++;$i<@ARGV;$i++) {
			if ($ARGV[$i] =~ /\-/) {
				$i--;
				last;
			}
			$force{$ARGV[$i]}=1;
		}	
		next;
	} elsif ($ARGV[$i] eq '-s') {
		for ($i++;$i<@ARGV;$i++) {
			if ($ARGV[$i] =~ /\-/) {
				$i--;
				last;
			}
			$skip{$ARGV[$i]}=1;
		}	
		next;
	} elsif ($ARGV[$i] eq '-A') {
		$i++;
		$Aoptions .= " ";
		for (;$i<@ARGV;$i++) {
			if ($ARGV[$i] eq '-A' || $ARGV[$i] eq '-B' || $ARGV[$i] eq '-C' || $ARGV[$i] eq '-D') {
				$i--;
				last;
			}
			$Aoptions .= " \"" . $ARGV[$i] . "\"";
		}
	} elsif ($ARGV[$i] eq '-B') {
		$i++;
		$Boptions .= " ";
		for (;$i<@ARGV;$i++) {
			if ($ARGV[$i] eq '-A' || $ARGV[$i] eq '-B' || $ARGV[$i] eq '-C' || $ARGV[$i] eq '-D') {
				$i--;
				last;
			}
			$Boptions .= " \"" . $ARGV[$i] . "\"";
		}
	} elsif ($ARGV[$i] eq '-C') {
		$i++;
		$Coptions .= " ";
		for (;$i<@ARGV;$i++) {
			if ($ARGV[$i] eq '-A' || $ARGV[$i] eq '-B' || $ARGV[$i] eq '-C' || $ARGV[$i] eq '-D') {
				$i--;
				last;
			}
			$Coptions .= " \"" . $ARGV[$i] . "\"";
		}
	} elsif ($ARGV[$i] eq '-D') {
		$i++;
		$Doptions .= " ";
		for (;$i<@ARGV;$i++) {
			if ($ARGV[$i] eq '-A' || $ARGV[$i] eq '-B' || $ARGV[$i] eq '-C' || $ARGV[$i] eq '-D') {
				$i--;
				last;
			}
			$Doptions .= " \"" . $ARGV[$i] . "\"";
		}
	} else {
		printCMD();
	}
}

if ($Boptions =~ /\"\-style\"\s+\"(.+?)\"/) {
	my $end = $1;
	$style = $end;
	print STDERR "\tDetected Boption style: $style\n";
}

print STDERR "\t------------------------------\n";


my $tmpFile = rand() . ".tmp";
my $tmpInput = '';


if (@irawFiles > 0) {
	print STDERR "\tCreating a temporary input directory:\n";
	$tmpInput = rand() . ".inputDir";
	my $str = '';
	foreach(@irawFiles) {
		$str .= " \"$_\"";
	}
	my $opt = "";
	$opt .= " -checkGC " if ($skipFreq!=1);
	`makeTagDirectory "$tmpInput" $str -genome $genome $opt`;
	$inputDir = $tmpInput;
}

print STDERR "\tPart A: Parsing Alignment\n"; 
if ($Aoptions eq '' || exists($skip{'A'})) {
	print STDERR "\t\tSkipping...\n";
} else {

	my $opt = "";
	$opt .= " -checkGC " if ($skipFreq!=1);
	`makeTagDirectory "$expDir" $Aoptions -genome "$genomeDir" $opt`;
}
if ($skipGO != 1) {
	my $goOptions = "";
	if ($inputDir ne '') {
		$goOptions = " -bg \"$inputDir\"";
	}
	if ($customGenome == 0) {
		`GenomeOntology.pl "$expDir" $genome "$expDir"/GenomeOntology-Tags/ $goOptions`;
	}
}

my %files;
`ls -1 "$expDir/" > $tmpFile`;
open IN, $tmpFile;
while (<IN>) {
	chomp;
	$files{$_}=$_;
}
close IN;
`rm $tmpFile`;

if (!exists($files{"tagInfo.txt"})) {
	print STDERR "!!! tagInfo.txt could not be found $expDir!!\n";
	print STDERR "Have you parsed your alignment files yet?  Try adding \"-A alignedReads.bed\"\n";
	print STDERR "Stopping.\n";
	exit;
}


open HTML, ">$expDir/index.html";
print HTML "<HTML><HEAD><TITLE>Results for $expDir</TITLE></HEAD>\n<BODY>\n";
print HTML "<H1>$expDir</H1>\n";
print HTML "<UL>\n";
print HTML "\t<LI><A HREF=\"tagInfo.txt\">Mapped Tag Information</A></LI>\n";
print HTML "\t<LI><A HREF=\"tagAutocorrelation.txt\">Tag Autocorrelation (open with Excel)</A></LI>\n";
print HTML "\t<LI><A HREF=\"tagCountDistribution.txt\">Clonal tag count distribution (open with Excel)</A></LI>\n";
print HTML "\t<LI><A HREF=\"tagLengthDistribution.txt\">Tag length distribution (open with Excel)</A></LI>\n";
if ($skipFreq ==0) {
	print HTML "\t<LI><A HREF=\"tagFreq.txt\">Nucleotide Frequency relative to tag position (open with Excel)</A></LI>\n";
	print HTML "\t<LI><A HREF=\"tagFreqUniq.txt\">Nucleotide Frequency relative to tag position [limit on tag per bp] (open with Excel)</A></LI>\n";
	print HTML "\t<LI><A HREF=\"tagGCcontent.txt\">Tag GC-content distribution (open with Excel)</A></LI>\n";
	print HTML "\t<LI><A HREF=\"genomeGCcontent.txt\">Expected Tag GC-content distribution [genome] (open with Excel)</A></LI>\n";
}
if ($skipGO == 0) {
	print HTML "\t<LI><A HREF=\"GenomeOntology-Tags/GenomeOntology.html\">Genome Ontology Results (based on raw tags)</A></LI>\n";
}




my $peakFile = "peaks.txt";
if ($style eq 'histone') {
	$peakFile = "regions.txt";
} elsif ($style eq 'tss') {
	$peakFile = "tss.txt";
} elsif ($style eq 'groseq') {
	$peakFile = "transcripts.txt";
} else {
	$peakFile = "peaks.txt";
}
my $annotatedPeaks = $peakFile;
$annotatedPeaks =~ s/\.txt/\.annotated\.txt/;
my $enhancerPeaks = $peakFile;
$enhancerPeaks =~ s/\.txt/\.enhancers\.txt/;
my $annotatedEnhancers = $enhancerPeaks;
$annotatedEnhancers =~ s/\.txt/\.annotated\.txt/;
my $focusedPeaks = $peakFile;
if ($style eq 'histone') {
	$focusedPeaks =~ s/\.txt/\.nfr\.txt/;
} else {
	$focusedPeaks =~ s/\.txt/\.focused\.txt/;
}


print STDERR "\n\tPart B: Finding Peaks and generating UCSC Bedgraph file\n";
if (exists($skip{'B'})) {
	print STDERR "\t\tSkipping...\n";
} else {
	$Boptions = " -style $style " . $Boptions;
	if ($inputDir ne '') {
		$Boptions = " -i \"$inputDir\" " . $Boptions;
	}
	print STDERR "`findPeaks $expDir $Boptions -o auto`;\n";
	`findPeaks "$expDir" $Boptions -o auto`;
	`makeUCSCfile "$expDir" -o auto`;

	my $peakInfoFile = "$expDir/peakInfo.txt";
	open OUT, ">$peakInfoFile";
	my $pfile = "$expDir/$peakFile";
	open IN, $pfile;
	while (<IN>) {
		if (/^#/) {
			print OUT $_;
		} else {
			last;
		}
	}
	close IN;
	close OUT;

	my $gVersion = $genome;
	$gVersion =~ s/r$//;
	
	`ls "$expDir"/*ucsc*gz > .ls`;
	my $ucscfile = '';
	open IN, ".ls";
	while (<IN>) {
		chomp;
		s/^.*\///;
		$ucscfile = $_;
		last;
	}
	close IN;
	`rm .ls`;

	#print HTML "\t<LI><A HREF=\"peakInfo.txt\">Peak Information</A> | ";
	print HTML "\t<A HREF=\"$peakFile\">Peak File (open with Excel)</A></LI>\n";
	print HTML "\t<A HREF=\"peakInfo.txt\">Region Info/Stats</A></LI>\n";

	print HTML "\t<LI><A HREF=\"http://genome.ucsc.edu/cgi-bin/hgGateway?db=$gVersion\">UCSC Genome Browser (genome = $gVersion)</A>\n";
	print HTML "\t<UL>\n";
	print HTML "\t<LI><A HREF=\"$_\">$ucscfile</A></LI>";
	print HTML "\t</UL></LI>\n";

	$focusInputPeaks = $peakFile;	
	if ($distalFlag && $customGenome == 0) {
		`getDistalPeaks.pl "$expDir/$peakFile" $genome > "$expDir/$enhancerPeaks"`;
		print HTML "\t<A HREF=\"$enhancerPeaks\">Enhancer/Distal Peaks File (more than 3kb from TSS) (Excel)</A></LI>\n";
		$focusInputPeaks = $enhancerPeaks;
	}
	if ($focusFlag) {
		if ($style eq 'histone') {
			`getPeakTags "$expDir/$focusInputPeaks" "$expDir/" -center -nfr > "$expDir/$focusedPeaks"`;
			print HTML "\t<A HREF=\"$focusedPeaks\">Nucleosome Free Region centered peaks (Excel)</A></LI>\n";
		} else {
			`getFocalPeaks.pl "$expDir/$focusInputPeaks" 0.9 > "$expDir/$focusedPeaks"`;
			print HTML "\t<A HREF=\"$focusedPeaks\">Focal Peak File (more than 90% focused) (Excel)</A></LI>\n";
		}
	}
}

my $peaksForMotifs = $peakFile;
if ($distalFlag) {
	$peaksForMotifs = $enhancerPeaks;
}


print STDERR "\n\tPart C: Motif Finding\n";
if (!exists($skip{'C'})) {
	if (exists($files{"homerResults.html"}) && !exists($force{'C'}) && $Coptions eq '') {
		print STDERR "\t\tSkipping... (found homerResults.html)\n";
		print HTML "\t<LI><A HREF=\"homerResults.html\"><i>De novo</i> Motif Results</A> | \n";
		print HTML "<A HREF=\"knownResults.html\">Known Motif Enrichment</A></LI>\n";
	} else {
		if ($msize eq '') {
			if ($style eq 'factor')  {
				$msize = 200;
			} elsif ($style eq 'histone') {
				$msize = 1000;
			} elsif ($style eq 'tss') {
				$msize = "-150,50";
			} else {
				$msize = 200;
			}
		}
		my $dirname = "Motifs-" . $msize . "bp";
		if ($Coptions eq '') {
			$Coptions = " -len 8,10,12 ";
		}
		if ($cpgFlag) {
			$Coptions .= " -cpg";
		}

		`findMotifsGenome.pl "$expDir/$peaksForMotifs" "$genomeDir" "$expDir/$dirname" $maskFlag -size $msize -p $numCPUs $Coptions`;
		print HTML "\t<LI><A HREF=\"$dirname/homerResults.html\"><i>De novo</i> Motif Results</A> | \n";
		print HTML "<A HREF=\"$dirname/knownResults.html\">Known Motif Enrichment (Distal)</A></LI>\n";
		
		if ($focusFlag ==1) {
			if ($style eq 'histone') {
				my $fsize = 200;
				$dirname = "Motifs-NFR-" . $fsize . "bp";
				`findMotifsGenome.pl "$expDir/$focusedPeaks" "$genomeDir" "$expDir/$dirname" $maskFlag -p $numCPUs $Coptions -size $fsize`;
				print HTML "\t<LI><A HREF=\"$dirname/homerResults.html\">Nucleosome Free Region <i>De novo</i> Motif Results</A> | \n";
				print HTML "<A HREF=\"$dirname/knownResults.html\">Nucleosome Free Region Known Motif Enrichment</A></LI>\n";
			} else {
				my $fsize = 50;
				$dirname = "Motifs-Focused-" . $fsize . "bp";
				`findMotifsGenome.pl "$expDir/$focusedPeaks" "$genomeDir" "$expDir/$dirname" $maskFlag -p $numCPUs $Coptions -len 6,8,10,12,15,20,30 -S 5 -size $fsize`;
				print HTML "\t<LI><A HREF=\"$dirname/homerResults.html\">Focused <i>De novo</i> Motif Results</A> | \n";
				print HTML "<A HREF=\"$dirname/knownResults.html\">Focused Known Motif Enrichment</A></LI>\n";
			}
		}
	}
}


my $inputOption = "";
if ($inputDir ne '') {
	$inputOption = "\"$inputDir\"";
}

print STDERR "\n\tPart D: Annotate Peaks\n";
if ((exists($files{$annotatedPeaks}) && !exists($force{'D'})) || exists($skip{'D'})) {
	print STDERR "\t\tSkipping...\n";
} else {

	`annotatePeaks.pl "$expDir/$peakFile" "$genomeDir" -d "$expDir" $inputOption -go "$expDir/GOanalysis" -genomeOntology "$expDir/GenomeOntology-Peaks/" $Doptions > "$expDir/$annotatedPeaks"`;
	print HTML "\t<LI><A HREF=\"$annotatedPeaks\">Annotated Peak File (Excel)</A></LI>\n";
	print HTML "\t<LI><A HREF=\"GenomeOntology-Peaks/GenomeOntology.html\">Genome Ontology based on Peaks</A></LI>\n";
	print HTML "\t<LI><A HREF=\"GOanalysis/geneOntology.html\">GO analysis of bound genes</A></LI>\n";

	if ($distalFlag) {
		`annotatePeaks.pl "$expDir/$enhancerPeaks" "$genomeDir" -d "$expDir" $inputOption -go "$expDir/GOanalysis-Enhancers" -genomeOntology "$expDir/GenomeOntology-EnhancerPeaks/" $Doptions > "$expDir/$annotatedEnhancers"`;
		print HTML "\t<LI><A HREF=\"$annotatedEnhancers\">Enhancer/Distal Annotated Peak File (Excel)</A></LI>\n";
		print HTML "\t<LI><A HREF=\"GenomeOntology-EnhancerPeaks/GenomeOntology.html\">Genome Ontology based on Enhancer Peaks</A></LI>\n";
		print HTML "\t<LI><A HREF=\"GOanalysis-Enhancers/geneOntology.html\">GO analysis of genes with Enhancer peaks</A></LI>\n";
	}
}

print HTML "</UL></BODY></HTML>\n";
close HTML;


if ($tmpInput ne '') {
	`rm -r "$tmpInput"`;
}

