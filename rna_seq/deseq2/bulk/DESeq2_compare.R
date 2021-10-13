#!/usr/bin/Rscript

message("Loading required libraries...\n")
library('DESeq2')
library('impute')
library('pcaExplorer')

#import raw counts and coldata matrices
message("Reading in data...\n")
cts_raw <- as.matrix(read.table('../counts.tsv', header=TRUE, sep="\t", row.names=1))
coldata <- as.matrix(read.table('../coldata.txt', header=TRUE, sep="\t", row.names=1))
samples <- readLines('./samples.txt')

#reduce lists to only samples we are interested in
message("Reducing data to only provided sample list...\n")
cts_raw <- cts_raw[, samples]
coldata <- coldata[samples, ]
coldata <- data.frame(coldata)
colnames(coldata) <- c('Disease')

#read in arguments for plot titles and comparison
args <- commandArgs(TRUE)
comparison = c(args[3],args[1],args[2])

#sanity checks
message("Performing sanity checks...\n")
if(! all(rownames(coldata) %in% colnames(cts_raw))) {
    stop('!!! Not all row names of coldata are in counts.tsv !!! Aborting...')
} else if (! all(rownames(coldata) == colnames(cts_raw))) {
    stop('!!! Row names of coldata are not in the same order as counts.tsv !!! Aborting...')
}


#begin DESeq2 comparisons
message("Running DESeq2 comparison on raw data...\n")
#head(cts_raw)
dds <- DESeqDataSetFromMatrix(countData = cts_raw, colData = coldata, design = ~ Disease)
dds <- DESeq(dds)
res <- results(dds)
rlt <- DESeq2::rlogTransformation(dds)
res <- results(dds, contrast = comparison)
write.table(counts(dds, normalized=TRUE), 'raw_norm_counts.txt', sep="\t", quote=FALSE)
write.table(res, 'raw.results.txt', sep="\t", quote=FALSE)
pcaplot(rlt, intgroup = c(args[3]), ntop = 1000, pcX = 1, pcY = 2, title = paste("Samples PCA: ", args[1]," vs. ", args[2], " (All Transcripts)", sep = ""), ellipse = TRUE, text_labels=TRUE, point_size = 4)
pcaplot(rlt, intgroup = c(args[3]), ntop = 1000, pcX = 1, pcY = 2, title = paste("Samples PCA: ", args[1]," vs. ", args[2]," vs. ", args[4], " (All Transcripts)", sep = ""), ellipse = FALSE, text_labels=FALSE, point_size = 4)
