#!/usr/bin/env perl
use warnings;

sub printCMD() {
	print STDERR "\n\tUsage: map-bowtie2.pl [options] <FASTQ file1> [FASTQ file2]...\n";
	print STDERR "\n\t\t\tpaired end: <FASTQ end1>,<FASTQ end2> [FASTQ file2 end1],[end2]...\n";
	print STDERR "\n\tRequired Options:\n";
	print STDERR "\t\t-index <path-to-bt2 index> (path to bowtie2 index to use for mapping)\n";
	print STDERR "\t\t-x <path...> (also works same as -index)\n";
	print STDERR "\t\t\ti.e. -index /bioinf/bowtie2/bowtie2-2.0.0-beta6/indexes/hg19\n";
	print STDERR "\n\tAlignment Type:\n";
	print STDERR "\t\t-bowtie2 (map with bowtie2,default)\n";
	print STDERR "\t\t-tophat2 (map with tophat2)\n";
	print STDERR "\t\t\t-path <path-to-program-file> (executable file to run if not in path/diff name)\n";
	print STDERR "\n\tCPU options:\n";
	print STDERR "\t\t-cpu <#> (Number of instances to run at once, default:1)\n";
	print STDERR "\t\t-p <#> (Number of cpus per instance, default: 1)\n";
	print STDERR "\t\t-un (will output unaligned reads)\n";
	print STDERR "\t\tRecommended examples (for 8 cores):\n";
	print STDERR "\t\t\tFor bowtie2: -cpu 1 -p 8 (align each file quickly 1 at a time)\n";
	print STDERR "\t\t\tFor tophat2: -cpu 8 -p 1 (because tophat has several non-parallel steps...)\n";
	print STDERR "\n\tBowtie options:\n";
	print STDERR "\t\t--local (local alignment, default: global/end-to-end)\n";
	print STDERR "\t\t-bam (convert output file to sorted bam file, default: sam)\n";
	print STDERR "\n\tTophat options:\n";
	print STDERR "\t\t--library-type <type> (library type for tophat, default: fr-firststrand\n";
	print STDERR "\t\t\t\tOther optoins: fr-unstranded, fr-secondstrand)\n";
	print STDERR "\t\t-G <GTF file> (Use as guide for splice mapping)\n";
	#print STDERR "\t\t-mis <#> (Max number of mismatches, default: 2)\n";
	print STDERR "\n\tGeneral Options to pass along to alignment program:\n";
	print STDERR "\t\t-f (Input is FASTA files, default expects FASTQ)\n";
	print STDERR "\t\t-pass \"...\" (need to include quotes)\n";
	#print STDERR "\t\t-remap (remap unaligned reads with bowtie afterwards, returning random position)\n";
	print STDERR "\n";
	exit;
}

if (@ARGV < 2) {
	printCMD();
}

my $bt2Index = "";
my $maxCPUs = 1;
my $pCPUs = 1;
my $remapFlag = 0;
my @files = ();
my $libraryType = "fr-firststrand";
my $maxMultihits = 20;
my $guideGTF = "";
my $program = 'bowtie2';
my $exe = '';
my $local = '';
my $pass = "";
my $unFlag = 0;
my $n = 3;
my $maxMisMatches = " -n $n --genome-read-mismatches $n --read-mismatches $n ";
my $bamFlag = 0;

for (my $i=0;$i<@ARGV;$i++) {
	#print STDERR "$ARGV[$i] $i\n";
	if ($ARGV[$i] eq '-index' || $ARGV[$i] eq '-x') {
		$bt2Index = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-bam') {
		$bamFlag= 1;
	} elsif ($ARGV[$i] eq '-mis') {
		my $n = $ARGV[++$i];
		$maxMisMatches = " -n $n --genome-read-mismatches $n --read-mismatches $n ";
	} elsif ($ARGV[$i] eq '-bowtie2') {
		$program = 'bowtie2';
	} elsif ($ARGV[$i] eq '-tophat2') {
		$program = 'tophat2';
	} elsif ($ARGV[$i] eq '-f') {
		$pass .= " -f ";
	} elsif ($ARGV[$i] eq '--local') {
		$local = '--local';
	} elsif ($ARGV[$i] eq '-path') {
		$exe = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-cpu') {
		$maxCPUs = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-un') {
		$unFlag = 1;
	} elsif ($ARGV[$i] eq '-p') {
		$pCPUs = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-remap') {
		$remapFlag = 1;
	} elsif ($ARGV[$i] eq '--library-type') {
		$libraryType = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-pass') {
		$pass = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-G') {
		$guideGTF = " -G $ARGV[++$i]";
	} elsif ($ARGV[$i] =~ /^\-/) {
		print STDERR "\n!!! \"$ARGV[$i]\" not recognized...\n";
		printCMD();
	} else {
		push(@files, $ARGV[$i]);
	}
}
#$program .= " -f ";
if ($exe eq '') {
	$exe = $program;
	print STDERR "\tWill run $exe...\n";
}
if ($bt2Index eq '') {
	print STDERR "!!! Need to specify bowtie2 index and/or CPUs !!!\n";
	printCMD();	
}

my $genomeName = $bt2Index;
$genomeName =~ s/^.+\///;

print STDERR "\tbt2Index = $bt2Index ($genomeName)\n";
print STDERR "\tNumber of instances at once: $maxCPUs\n";
print STDERR "\tNumber of cpus per instance: $pCPUs\n";
if ($program eq 'bowtie2') {
	print STDERR "\n\tbowtie2 will be used to align the following files:\n";
} elsif ($program eq 'tophat2') {
	print STDERR "\tlibrary type: $libraryType\n";
	print STDERR "\n\tTophat2 will be used to align the following files:\n";
	$maxMisMatches = '';
}
foreach(@files) {
	print STDERR "\t\t$_\n";
}
print STDERR "\t\n";

my @pids = ();
my $cpus = 0;
for (my $j=0;$j<@files;$j++) {
	my $file = $files[$j];
	my $pid = fork();
	$cpus++;
	if ($pid == 0) {
		#child process

		my $peflag = 0;
		my %delete = ();
		my $file1 = "";
		my $file2 = "";
		if ($file =~ /\,/) {
			my @f = split /\,/, $file;
			$file1 = $f[0];
			$file2 = $f[1];
			print STDERR "\t\tPaired end files: $file1 and $file2\n";
			$peflag = 1;
			if ($file1 =~ /\.gz$/) {
				my $newfile = $file1;
				$newfile  =~ s/\.gz$//;
				`gunzip -c "$file1" > "$newfile"`;
				$delete{$newfile} = 1;
				$file1 = $newfile;
			}
			if ($file2 =~ /\.gz$/) {
				my $newfile = $file2;
				$newfile  =~ s/\.gz$//;
				`gunzip -c "$file2" > "$newfile"`;
				$delete{$newfile} = 1;
				$file2 = $newfile;
			}
			if ($program eq 'tophat2') {
				$file = "\"$file1\" \"$file2\"";
			} elsif ($program eq 'bowtie2') {
				$file = "-1 \"$file1\" -2 \"$file2\"";
			}
		} else {
			if ($file =~ /\.gz$/) {
				my $newfile = $file;
				$newfile  =~ s/\.gz$//;
				`gunzip -c "$file" > $newfile`;
				$delete{$newfile} = 1;
				$file = $newfile;
			}
			$file1 = $file;
			$file = "\"$file\"";
		}
		if ($program eq 'tophat2') {
			my $outputDir = "$file1.$genomeName.tophat2";
			my $outputFile = "$file1.$genomeName.tophat2.bam";
			my $outputJunc = "$file1.$genomeName.tophat2.junc";
			my $logFile = "$file1.$genomeName.tophat2.log";
			`$exe --library-type $libraryType -p $pCPUs -g $maxMultihits $guideGTF -o $outputDir $pass $maxMisMatches "$bt2Index" $file 2> $logFile`;
			`mv "$outputDir/accepted_hits.bam" "$outputFile"`;
        	`parseTophatJunctions.pl "$outputDir/junctions.bed" > "$outputJunc"`;
			`rm -r "$outputDir"`;
			if ($unFlag) {
				my $samFile = "$file1.$genomeName.tophat2.sam";
				`samtools view -h "$outputFile" > "$samFile" 2>> $logFile`;
				my $unalignedFile = "$file1.$genomeName.tophat2.unaligned.fq";
				my $unalignedFile2 = "$file2.$genomeName.tophat2.unaligned.fq";
				`getUnalignedReadsSam.pl $samFile $file1 > $unalignedFile 2>> $logFile`;
				if ($peflag) {
					`getUnalignedReadsSam.pl $samFile $file2 > $unalignedFile2 2>> $logFile`;
				}
				`rm "$samFile"`;
			}
		} else {
			my $outputFile = "$file1.$genomeName.bowtie2.sam";
			my $logFile = "$file1.$genomeName.bowtie2.log";
			my $unFile = " --no-unal";
			if ($unFlag) {
				$unFile = "--un $file1.$genomeName.bowtie2.unaligned.fq";
			}
			`$exe $local $unFile -p $pCPUs $pass -x "$bt2Index" $file > "$outputFile" 2> "$logFile"`;
			if ($bamFlag) {
				my $bamFile = "$file1.$genomeName.bowtie2.tmp.bam";
				my $sortedBamFile = "$file1.$genomeName.bowtie2";
				`samtools view -S -b -@ $pCPUs "$outputFile" > "$bamFile" 2>> "$logFile"`;
				`samtools sort -m 5000000000 -@ $pCPUs "$bamFile" "$sortedBamFile" 2>> $logFile`;
				`rm "$outputFile" "$bamFile"`;
			}
		}

		foreach(keys %delete) {
			`rm "$_"`;
		}
		exit(0);
	}
	push(@pids, $pid);
	if ($cpus >= $maxCPUs) {
		my $id = wait();
		$cpus--;
	}
}
my $id = 0;
while ($id >= 0) {
	$id = wait();
	if ($id == -1) {
	} else {
	}
}
print STDERR "\n";
