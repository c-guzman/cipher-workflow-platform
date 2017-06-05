Analysis-Mode
=============

Analysis mode provides several functions that either combine/integrate various NGS datasets (typically files output already output by CIPHER's processing workflows) or performs complex bioinformatic tasks.

You can think of analysis mode as an entirely independent set of 'functions' that each have their own flags/parameters. For specific requests please submit a issue to our GitHub with a "[FUNCTION REQUEST]" title.

**CURRENTLY IMPLEMENTED FUNCTIONS:**

1. Enhancer prediction ("--function predictEnhancers") - predicts potential enhancers using a combination of ChIP-seq (H3K4me3, H3K4me1, H3K27Ac) and DNase-seq data.

2. Calculate expression of genes near a set of peaks or features ("--function geneExpressionNearPeaks") - calculates the gene expression level of genes (FPKM and TPM) near a set of features (e.g. peaks, enhancers, etc.) from RNA-seq data and a BED formatted list of features.