#!/usr/bin/env Rscript
# Script to conduct differential gene expression analysis using RUVSeq

# Set up command line arguments
args <- commandArgs(trailingOnly=TRUE)

exp_file <- read.table(args[1], header=T, stringsAsFactors = F)
args <- args[-1]

# Load libraries
library(RUVSeq)
library(edgeR)
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
colnames(data) <- gsub(".mapped.bam", "", colnames(data))

# Set GeneID as rowname
rownames(data) <- data[,1]
data[,1] <- NULL

sampleTable <- (condition = factor(exp_file$condition))

# Filter out non-expressed genes, require more than 5 reads in at least two samples
filter <- apply(data, 1, function(x) length(x[x>5])>=2)
filtered <- data[filter,]

x <- sampleTable
set <- newSeqExpressionSet(as.matrix(filtered),
                           phenoData = data.frame(x, row.names=colnames(filtered)))
set <- betweenLaneNormalization(set, which="upper")

# Obtain a set of empirical negative controls (i.e. the least significantly DE genes based on first pass DGE)
design <- model.matrix(~x, data=pData(set))
y <- DGEList(counts=counts(set), group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)

top <- topTags(lrt, n=nrow(set))$table
empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:5000]))] # By default the top 5000 least DEG are used as empirical controls, change this number

set2 <- RUVg(set, empirical, k=1) # Alter k = number of variances also

# Plot some PCA and RLE datasets for QC
colors <- brewer.pal(3, "Set2")

pdf('RLE_plot.pdf')
rleplot <- plotRLE(set2, outline=F, ylim=c(4,-4), col=colors[x])
dev.off()

pdf('PCA_plot.pdf')
pcaplot <- plotPCA(set2, col=colors[x], cex=1.2)
dev.off()

# DGE
design <- model.matrix(~x + W_1, data=pData(set2))
y <- DGEList(counts=counts(set2), group=x)
y <- calcNormFactors(y, method = "upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)

# Select top genes dge with a p-value of 0.05
top.df <- topTags(lrt, n=Inf, sort.by="logFC", p.value = 0.05)

write.table(top.df, "edgeR.dgeList.txt", row.names = T, quote = F, sep = '\t')