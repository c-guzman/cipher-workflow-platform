Installation
============

CIPHER requires no manual installation. Just download and unzip the tar file from our GitHub or git clone the repository. The download may be large.

The only required manual software installation on your local computer / cluster are Docker and Nextflow.

By default, CIPHER will run off a Docker container, but this can be optionally turned off by removing the option in the config file. For more information, please read the 'Nextflow Config File' documentation.

Installing Nextflow
===================

For Linux:

1. We recommend installing Nextflow through the Anaconda package manager.
    ::

        wget https://repo.continuum.io/archive/Anaconda2-4.3.1-Linux-x86_64.sh
        bash Anaconda2-4.3.1-Linux-x86_64.sh

2. Install Nextflow
    ::
    
        conda install -c bioconda nextflow


For macOS:

1. We recommend installing Nextflow through the Anaconda package manager.
    ::

        wget https://repo.continuum.io/archive/Anaconda2-4.3.1-MacOSX-x86_64.sh
        bash Anaconda2-4.3.1-Linux-x86_64.sh

2. Install Nextflow
    ::

        conda install -c bioconda nextflow

For Windows:

CIPHER has not been tested on the Windows operating system. However, Windows 10 has introduced the 'Ubuntu Bash Shell' sub-system, which potentially can be used to run CIPHER.

1. Install Ubuntu Bash on your Windows 10 computer. Follow the instructions here: https://www.howtogeek.com/249966/how-to-install-and-use-the-linux-bash-shell-on-windows-10/

2. Using the Ubuntu Bash terminal we recommend installing nextflow through the Anaconda package manager.
    ::

        wget https://repo.continuum.io/archive/Anaconda2-4.3.1-Linux-x86_64.sh
        bash Anaconda2-4.3.1-Linux-x86_64.sh

3. Install Nextflow
    ::

        conda install -c bioconda nextflow

Installing Docker
=================

Installing Docker can be tricky. Because the instructions for correct installation are constantly changing, we provide links to Docker's official installation process. You'll want to use Linux's installation instructions for Windows via the Ubuntu Bash terminal.

For Linux:

https://docs.docker.com/engine/installation/#supported-platforms

For macS:

https://docs.docker.com/docker-for-mac/

Manual Installation
===================

Some universities may not want to install Docker for security reasons. And others may not want to use Docker at all on their local desktop. This section lists how to manually install all the dependencies for CIPHER.

We recommend using the Anaconda package manager.

For Linux:

1. We recommend installing Nextflow through the Anaconda package manager.
    ::

        wget https://repo.continuum.io/archive/Anaconda2-4.3.1-Linux-x86_64.sh
        bash Anaconda2-4.3.1-Linux-x86_64.sh

2. Install bioconda packages
    ::

        conda install -c bioconda nextflow fastqc bbmap star hisat2 bowtie2 bwa multiqc macs2 deeptools epic preseq samtools sambamba bedtools bedops stringtie subread

3. Install R packages
    ::

        R

        install.packages(c("data.table", "ggplot2", "gplots"))
        install.packages("http://hartleys.github.io/QoRTs/QoRTs_LATEST.tar.gz", repos=NULL, type="source")
        source("https://bioconductor.org/biocLite.R")
        biocLite()
        biocLite(c("ChIPQC", "RUVSeq", "ChIPseeker"))


For macOS:

1. We recommend installing Nextflow through the Anaconda package manager.
    ::

        wget https://repo.continuum.io/archive/Anaconda2-4.3.1-MacOSX-x86_64.sh
        bash Anaconda2-4.3.1-Linux-x86_64.sh

2. Install bioconda packages
    ::

        conda install -c bioconda nextflow fastqc bbmap star hisat2 bowtie2 bwa multiqc macs2 deeptools epic preseq samtools sambamba bedtools bedops stringtie subread

3. Install R packages
    ::

        R

        install.packages(c("data.table", "ggplot2", "gplots"))
        install.packages("http://hartleys.github.io/QoRTs/QoRTs_LATEST.tar.gz", repos=NULL, type="source")
        source("https://bioconductor.org/biocLite.R")
        biocLite()
        biocLite(c("ChIPQC", "RUVSeq", "ChIPseeker"))
