#!/usr/bin/env Rscript
# Script to conduct differential gene expression analysis using RUVSeq

# Set up command line arguments
args <- commandArgs(trailingOnly=TRUE)

exp_file <- read.table(args[1], header=T, stringsAsFactors = F)
args <- args[-1]

# Load libraries
library(RUVSeq)
library(edgeR)
library(DESeq2)
library(data.table)
library(RColorBrewer)

# Read in count data from all files into a df
# Row names are geneids
temp <- lapply(args, fread, skip="Geneid", header=TRUE, colClasses=c(NA, rep("NULL", 5), NA))

# Merge all count data into a single df
merge.all <- function(x, y) {
  merge(x, y, all=TRUE, by="Geneid")
}

data <- data.frame(Reduce(merge.all, temp))

# Clean sample name headers
colnames(data) <- gsub(".sorted.mapped.bam", "", colnames(data))

# Set GeneID as rowname
rownames(data) <- data[,1]
data[,1] <- NULL

sampleTable <- (condition = factor(exp_file$condition))

# Filter out non-expressed genes, require more than 5 reads in at least two samples
filter <- apply(data, 1, function(x) length(x[x>5])>=2)
filtered <- data[filter,]
genes <- rownames(filtered)

x <- sampleTable
set <- newSeqExpressionSet(as.matrix(filtered),
                           phenoData = data.frame(x, row.names=colnames(filtered)))
set <- betweenLaneNormalization(set, which="upper")

# Use RUVs method to estimate the factors of undwanted variation using replicate/negative control samples
differences = makeGroups(x)
set2 = RUVs(set, genes, k=1, differences) # do QC

# Differential gene expression analysis using edgeR
design = model.matrix(~x + W_1, data = pData(set_ruvs))
y = DGEList(counts=counts(set_ruvs), group=x)
y = calcNormFactors(y, method="upperquartile")
y = estimateGLMCommonDisp(y, design)
y = estimateGLMTagwiseDisp(y, design)

fit = glmFit(y, design)
lrt = glmLRT(fit, coef=2)

pdf('RLE_plot.pdf')
rleplot <- plotRLE(set2, outline=F, ylim=c(4,-4), col=colors[x])
dev.off()

pdf('PCA_plot.pdf')
pcaplot <- plotPCA(set2, col=colors[x], cex=1.2)
dev.off()

# DGE edgeR
design <- model.matrix(~x + W_1, data=pData(set2))
y <- DGEList(counts=counts(set2), group=x)
y <- calcNormFactors(y, method = "upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)

# Select all top genes dge with a p-value of 1
top.df <- topTags(lrt, n=Inf, sort.by="logFC", p.value = 1)

# DGE DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts(set2), colData = pData(set2), design = ~ W_1 + x)
dds <- DESeq(dds)
res <- results(dds)
res_df <- as.data.frame(results)

# Plot MA figure for DESeq2
pdf('DESeq2_MAplot.pdf')
plotMA(dds, pvalCutoff=0.05)
dev.off()

write.table(res_df, "DESeq2.dgeList.txt", row.names = T, quote = F, sep = '\t')
write.table(top.df, "edgeR.dgeList.txt", row.names = T, quote = F, sep = '\t')