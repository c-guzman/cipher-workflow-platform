#!/usr/bin/env Rscript
# Script to analyze GROseq data using groHMM workflow

args <- commandArgs(trailingOnly = F)
# args[1] = bam file ${bam}
# args[2] = name of output file mergeID
# args[3] = ${params.threads}

# Load groHMM package and set up parallel processing
library(groHMM)
options(mc.cores=getCores(args[3]))

# Read in merged BAM files
bam_file <- as(readGAlignments(args[1]), "GRanges")

# Call transcripts from merged BAM
Sall <- sort(bam_file)
hmmResult <- detectTranscripts(Sall, LtProbB=-200, UTS=5, threshold=1)
txHMM <- as.data.frame(hmmResult$transcripts)

write.table(x = txHMM, file=paste(args[2],"_transcripts.txt", sep=""),
            quote = F, sep = '\t', row.names = F)