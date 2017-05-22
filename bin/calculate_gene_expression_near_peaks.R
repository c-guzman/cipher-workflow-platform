#!/usr/bin/env Rscript
# A script to calculate expression levels of genes near ChIPseq peaks

# Load packages
library(GenomicFeatures)
library(ChIPseeker)

# args[1] GTF File
# args[2] Peak File
# args[3] Stringtie output (gene_abundance.txt file)
txdb = makeTxDbFromGFF(args[1], format = "gtf", dbxrefTag = "GeneID")
peaks = readPeakFile(args[2], header = F)
expression_data = read.table(args[3], header = T, sep = '\t')
colnames(expression_data) = c("geneId", "GeneName", "Chr", "Strand", "Start", "End", "Coverage", "FPKM", "TPM")

# Peak annotation
peakAnno = annotatePeak(peaks, tssRegion = c(-3000, 3000),
                        TxDb=txdb)

# Conver peakAnno to data frame object
peaks_df = as.data.frame(peakAnno)

# Merge data together
merged = merge(peaks_df, expression_data, by = "geneId")

# Write out results
write.table(merged, "geneExpression_near_ChIPseq_peaks.txt", quote = F, row.names = F, sep = '\t')