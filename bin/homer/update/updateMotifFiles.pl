#!/usr/bin/perl -w

if (@ARGV < 1) {
#	print STDERR "\tNo input required...\n";
}
print STDERR "\tNo input required...\n";
print STDERR "\n\t!!! If wget has trouble downloading JASPAR matrices, you can download them yourself\n";
print STDERR "\tPlace them in the current directory, filenames:\n";
print STDERR "\t\tpfm_vertebrates.txt\n";
print STDERR "\t\tpfm_plants.txt\n";
print STDERR "\t\tpfm_insects.txt\n";
print STDERR "\t\tpfm_nematodes.txt\n";
print STDERR "\t\tpfm_fungi.txt\n";

my $activeMotifs = "../motifs/";
my $homerMotifs = "../data/knownTFs/";

`mkdir -p $homerMotifs`;
`mkdir -p $homerMotifs/all`;
`mkdir -p $homerMotifs/motifs`;
`cp -r $activeMotifs/* $homerMotifs/motifs/`;

my @groups = (
	"vertebrates",
	"insects",
	"worms",
	"plants",
	"yeast"
);


#my $jasparServer = "http://jaspar.binf.ku.dk/";
my $jasparServer = "http://jaspar.genereg.net/";

# Update JASPAR motifs
unless (-e "pfm_vertebrates.txt") {
	`wget -O pfm_vertebrates.txt $jasparServer/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_vertebrates.txt`;
	`./motifs/parseJasparMatrix.pl pfm_vertebrates.txt > motifs/vertebrates/jaspar.motifs`;
	`rm -f pfm_vertebrates.txt`;
} else {
	`./motifs/parseJasparMatrix.pl pfm_vertebrates.txt > motifs/vertebrates/jaspar.motifs`;
}

unless (-e "pfm_plants.txt") {
	`wget -O pfm_plants.txt $jasparServer/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_plants.txt`;
	`./motifs/parseJasparMatrix.pl pfm_plants.txt > motifs/plants/jaspar.motifs`;
	`rm -f pfm_plants.txt`;
} else {
	`./motifs/parseJasparMatrix.pl pfm_plants.txt > motifs/plants/jaspar.motifs`;
}
	

unless (-e "pfm_insects.txt") {
	`wget -O pfm_insects.txt $jasparServer/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_insects.txt`;
	`./motifs/parseJasparMatrix.pl pfm_insects.txt > motifs/insects/jaspar.motifs`;
	`rm -f pfm_insects.txt`;
} else {
	`./motifs/parseJasparMatrix.pl pfm_insects.txt > motifs/insects/jaspar.motifs`;
}

unless (-e "pfm_nematodes.txt") {
	`wget -O pfm_nematodes.txt $jasparServer/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_nematodes.txt`;
	`./motifs/parseJasparMatrix.pl pfm_nematodes.txt > motifs/worms/jaspar.motifs`;
	`rm -f pfm_nematodes.txt`;
} else {
	`./motifs/parseJasparMatrix.pl pfm_nematodes.txt > motifs/worms/jaspar.motifs`;
}

unless (-e "pfm_fungi.txt") {
	`wget -O pfm_fungi.txt $jasparServer/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_fungi.txt`;
	`./motifs/parseJasparMatrix.pl pfm_fungi.txt > motifs/yeast/jaspar.motifs`;
	`rm -f pfm_fungi.txt`;
} else {
	`./motifs/parseJasparMatrix.pl pfm_fungi.txt > motifs/yeast/jaspar.motifs`;
}

# known motifs first (ones that have degeneracy thresholds)
foreach(@groups) {
	my $g = $_;
	`mkdir -p $homerMotifs/$g`;
	if ($g eq 'vertebrates') {
		`cat $homerMotifs/motifs/*.motif > $homerMotifs/$g/known.motifs`;
	} else {
		`cat $homerMotifs/motifs/$g/*.motif > $homerMotifs/$g/known.motifs`;
	}
	`cat $homerMotifs/$g/known.motifs motifs/common.motifs motifs/$g/*.motifs > $homerMotifs/$g/all.motifs`;
}
`cat $homerMotifs/motifs/*.motif $homerMotifs/motifs/*/*.motif > $homerMotifs/all/known.motifs`;


`cat $homerMotifs/all/known.motifs motifs/common.motifs motifs/*/*.motifs > $homerMotifs/all/all.motifs`;
`cp $homerMotifs/all/known.motifs $homerMotifs/known.motifs`;
`cp $homerMotifs/all/all.motifs $homerMotifs/all.motifs`;
