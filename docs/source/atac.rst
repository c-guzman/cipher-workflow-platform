ATAC-seq Workflow
=================

All workflows require the following files:

1. A config file (described below)
2. Reference genome FASTA file
3. Reference genome GTF file

*Config File*
A config file is a tab separated text file that includes information regarding the name, location, and input of your experiment.

**Single-End Config**
This is the config file format for single-ended data.
    ::

        sample1     sample1_rep1    /path/to/sample1_rep1.fastq.gz  -    sample1
        sample1     sample1_rep2    /path/to/sample1_rep2.fastq.gz  -    sample1
        sample2     sample2_rep1    /path/to/sample2_rep1.fastq.gz  control2    sample2
        sample2     sample2_rep2    /path/to/sample2_rep2.fastq.gz  control2    sample2

The columns represent:

1. MergeID: The merge ID that will be used should your files be merged together. Should be the same for all replicates.
2. ID: The ID that will be used to name the majority of your files that are not merged. Recommended to be used to differentiate between different technical replicates.
3. Path: The path to the fastq file to be processed. Can be gzipped or not.
4. ControlID: The ID indicating what control file to be used for peak calling and other downstream analysis. Use "-" (without quotes) if there is no control for a particular sample.
5. Mark: The ID that signifies the type of mark or histone being processed. Use "input" if the line refers to a control. If the line is NOT a control, then use the MergeID name.

**Pair-End Config**
This is the config file format for pair-ended data.
    ::

        sample1		sample1_rep1	/path/to/sample1_rep1_R1.fastq.gz /path/to/sample1_rep1_R2.fastq.gz	-	sample1
        sample1		sample1_rep2	/path/to/sample1_rep2_R1.fastq.gz /path/to/sample1_rep2_R2.fastq.gz	-	sample1
        sample2		sample2_rep1	/path/to/sample2_rep1_R1.fastq.gz /path/to/sample2_rep1_R2.fastq.gz	-	sample2
        sample2		sample2_rep2	/path/to/sample2_rep2_R1.fastq.gz /path/to/sample2_R2_rep1.fastq.gz	-	sample2

The columns represent:

1. MergeID: The merge ID that will be used should your files be merged together. Should be the same for all replicates.
2. ID: The ID that will be used to name the majority of your files that are not merged. Recommended to be used to differentiate between different technical replicates.
3. Path1: The path to the first fastq file to be processed. Can be gzipped or not.
4. Path2: The path to the second fastq file to be processed. Can be gzipped or not.
5. ControlID: The ID indicating what control file to be used for peak calling and other downstream analysis. Use "-" (without quotes) if there is no control for a particular sample.
6. Mark: The ID that signifies the type of mark or histone being processed. Use "input" if the line refers to a control. If the line is NOT a control, then use the MergeID name.

**Simple MNase-seq Tutorial (single-end, 75 length reads)**
    ::

        nextflow run /path/to/main.nf --mode atac --config /path/to/config.txt --fasta /path/to/fasta.fa --gtf /path/to/gtf.gtf --lib s --readLen 75

**Simple MNase-seq Tutorial (pair-end, 75 length reads)**
    ::

        nextflow run /path/to/main.nf --mode atac --config /path/to/config.txt --fasta /path/to/fasta.fa --gtf /path/to/gtf.gtf --lib p --readLen 75

**Simple MNase-seq Tutorial (single-end, 75 length reads, use bowtie2 aligner instead of default bbmap, use 5 threads)**
    ::
    
        nextflow run /path/to/main.nf --mode atac --config /path/to/config.txt --fasta /path/to/fasta.fa --gtf /path/to/gtf.gtf --lib s --readLen 75 --aligner bowtie2 --threads 5