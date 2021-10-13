#!/usr/bin/Rscript

# Authored by Altan Turkoglu (turkoglu.12@osu.edu) and Alex Pan (alexander.pan@osumc.edu)
# Usage: Rscript DESeq2_covariateShrink.R [COLDATA_COL_NAME] [GROUP1] [GROUP2]

message("Loading required libraries...\n")
library('DESeq2')
library('impute')
#library('pcaExplorer')
library('apeglm')

#import raw counts and coldata matrices
message("Reading in data...\n")
##################################################
# MAKE SURE SAMPLES ARE IN CORRECT LOCATIONS
cts_raw <- as.matrix(read.table('../counts.tsv', header=TRUE, sep="\t", row.names=1))
cts_clear <- as.matrix(read.table('./clear.counts.tsv', header=TRUE, sep="\t", row.names=1))
coldata <- as.matrix(read.table('./coldata.txt', header=TRUE, sep="\t", row.names=1))
samples <- readLines('./samples.txt')
##################################################

#reduce lists to only samples we are interested in
message("Reducing data to only provided sample list...\n")
cts_raw <- cts_raw[, samples]
cts_clear <- cts_clear[, samples]
coldata <- coldata[samples, ]
# Replaces colnames with custom names
rownames(coldata) <- coldata[,2]
colnames(cts_raw) <- coldata[,2]
colnames(cts_clear) <- coldata[,2]

coldata <- data.frame(coldata)
##################################################
# MAKE SURE TO PASS THE SAME COLUMN NAMES AS THE FIRST LINE OF COLDATA.TXT
colnames(coldata) <- c('Treatment', 'Patient', 'SampleName')
coldata
##################################################

#read in arguments for plot titles and comparison
args <- commandArgs(TRUE)
comparision_name <- args[1]
comparison = c(comparision_name,args[2],args[3])
comparison

#sanity checks
message("Performing sanity checks...\n")
if(! all(rownames(coldata) %in% colnames(cts_raw))) {
    stop('!!! Not all row names of coldata are in counts.tsv !!! Aborting...')
} else if (! all(rownames(coldata) == colnames(cts_raw))) {
    stop('!!! Row names of coldata are not in the same order as counts.tsv !!! Aborting...')
}

if(! all(rownames(coldata) %in% colnames(cts_clear))) {
    stop('!!! Not all row names of coldata are in clear.counts.tsv !!! Aborting...')
} else if (! all(rownames(coldata) == colnames(cts_clear))) {
    stop('!!! Row names of coldata are not in the same order as clear.counts.tsv !!! Aborting...')
}

#create list with any genes with NAs omitted
message("Omitting NAs...\n")
cts_passed <- cts_clear[complete.cases(cts_clear[,1:ncol(cts_clear)]), ]

#function for deleting rows with more than "n" NAs
delete.na <- function(DF, n = 0) {
    DF[rowSums(is.na(DF)) <= n,]
}

#begin DESeq2 comparisons
message("Running DESeq2 comparison on raw data...\n")
head(cts_raw)
dds <- DESeqDataSetFromMatrix(countData = cts_raw, colData = coldata, design = ~ Treatment + Patient)
dds <- DESeq(dds)
res <- results(dds)
rlt <- DESeq2::rlogTransformation(dds)
res <- results(dds, contrast = comparison)
write.table(counts(dds, normalized=TRUE), 'raw_norm_counts.txt', sep="\t", quote=FALSE)
write.table(res, 'raw.results.txt', sep="\t", quote=FALSE)
pcaplot(rlt, intgroup = comparision_name, ntop = 1000, pcX = 1, pcY = 2, ellipse = TRUE, text_labels=TRUE, point_size = 4)
pcaplot(rlt, intgroup = comparision_name, ntop = 1000, pcX = 1, pcY = 2, ellipse = TRUE, text_labels=FALSE, point_size = 4)
pcaplot(rlt, intgroup = c(args[3]), ntop = 1000, pcX = 1, pcY = 2, title = paste("Samples PCA: ", args[1]," vs. ", args[2]," vs. ", args[4], " (All Transcripts)", sep = ""), ellipse = FALSE, text_labels=FALSE, point_size = 4)

message("Running DESeq2 comparison on CLEAR-only data...\n")
dds <- DESeqDataSetFromMatrix(countData = cts_passed, colData = coldata, design = ~ Treatment + Patient)
dds <- DESeq(dds)
res1 <- lfcShrink(dds, coef=2, type="apeglm", svalue=TRUE)
res2 <- lfcShrink(dds, coef=2, type="apeglm", svalue=FALSE)
#res <- results(dds)
rlt <- DESeq2::rlogTransformation(dds)
#res <- results(dds, contrast = comparison)
write.table(res1, 'clear.results.svalue.txt', sep="\t", quote=FALSE)
write.table(res2, 'clear.results.pvalue.txt', sep="\t", quote=FALSE)
write.table(counts(dds, normalized=TRUE), 'clear_norm_counts.txt', sep="\t", quote=FALSE)
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
plotDispEsts(dds)
pcaplot(rlt, intgroup = comparision_name, ntop = 1000, pcX = 1, pcY = 2, ellipse = TRUE, text_labels=TRUE, point_size = 4)
pcaplot(rlt, intgroup = comparision_name, ntop = 1000, pcX = 1, pcY = 2, ellipse = TRUE, text_labels=FALSE, point_size = 4)
#pcaplot(rlt, intgroup = comparision_name, ntop = 1000, pcX = 2, pcY = 3, ellipse = TRUE, text_labels=TRUE, point_size = 4)
#pcaplot(rlt, intgroup = comparision_name, ntop = 1000, pcX = 3, pcY = 4, ellipse = TRUE, text_labels=TRUE, point_size = 4)
pcascree(prcomp(t(assay(rlt))), type="pev")
pcaplot(rlt, intgroup = c(args[3]), ntop = 1000, pcX = 1, pcY = 2, title = paste("Samples PCA: ", args[1]," vs. ", args[2]," vs. ", args[4], " (CLEAR-Only Transcripts)", sep = ""), ellipse = FALSE, text_labels=FALSE, point_size = 4)
