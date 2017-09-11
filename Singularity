Bootstrap: docker
From: ubuntu:16.04

%labels
Maintainer cag104@ucsd.edu
Version 1.0

%environment
export PATH="/opt/anaconda2/bin:$PATH"

%post

	apt-get update

	apt-get install -y wget bzip2 perl gawk

	ln -s /bin/tar /bin/gtar

	wget https://repo.continuum.io/archive/Anaconda2-4.4.0-Linux-x86_64.sh

	bash Anaconda2-4.4.0-Linux-x86_64.sh -b -p /opt/anaconda2

	/opt/anaconda2/bin/conda install -c r --yes r-base=3.3.2 r-essentials=1.5.2 r-devtools=1.12.0

	/opt/anaconda2/bin/conda install -c bioconda --yes bbmap samtools epic jellyfish sambamba deeptools macs2 bedtools bedops multiqc subread stringtie nextflow bowtie2 bwa hisat2 star fastqc gimmemotifs nucleoatac salmon r-spp=1.14 bioconductor-biocinstaller=1.24.0 bioconductor-edger=3.16.5 bioconductor-deseq2=1.14.1 bioconductor-chipseeker=1.10.0

	wget -P /opt/anaconda2/bin http://hartleys.github.io/QoRTs/QoRTs.jar

	/opt/anaconda2/bin/R --slave -e 'install.packages("http://hartleys.github.io/QoRTs/QoRTs_LATEST.tar.gz", repos=NULL, type="source")'
