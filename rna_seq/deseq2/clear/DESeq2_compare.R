#!/usr/bin/Rscript

# Authored by Altan Turkoglu (turkoglu.12@osu.edu) and Alex Pan (alexander.pan@osumc.edu)
# Wednesday, March 12th, 2019

message("Loading required libraries...\n")
library('DESeq2')
library('impute')
library('pcaExplorer')

#import raw counts and coldata matrices
message("Reading in data...\n")
cts_raw <- as.matrix(read.table('./counts.tsv', header=TRUE, sep="\t", row.names=1))
cts_clear <- as.matrix(read.table('./clear.counts.tsv', header=TRUE, sep="\t", row.names=1))
coldata <- as.matrix(read.table('./coldata.txt', header=TRUE, sep="\t", row.names=1))
samples <- readLines('./samples.txt')

#reduce lists to only samples we are interested in
message("Reducing data to only provided sample list...\n")
cts_raw <- cts_raw[, samples]
cts_clear <- cts_clear[, samples]
coldata <- coldata[samples, ]
coldata <- data.frame(coldata)
colnames(coldata) <- c("Treatment", "Patient", "SampleName")

#replace rownames with custom names
message("Replacing rownames with custom names...\n")
rownames(coldata) <- coldata[,"SampleName"]
colnames(cts_raw) <- coldata[,"SampleName"]
colnames(cts_clear) <- coldata[,"SampleName"]

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

#create expression matrix with a maximum value of 25% NAs
#max_NAs = length(samples) * 0.25 #calculate how many NAs are permissable
max_NAs = 2 #allowing for 2 of the 8 lcRNA samples to have NAs
cts_imputed <- delete.na(cts_clear, max_NAs)

#impute
message("Imputing data...\n")
cts_imputed <- round(as.matrix(impute.knn(as.matrix( cts_imputed ) )$data))

#more sanity checks
message("Performing more sanity checks...\n")
if(! all(rownames(coldata) %in% colnames(cts_passed))) {
    stop('!!! Not all row names of coldata are in object cts_passed !!! Aborting...')
} else if (! all(rownames(coldata) == colnames(cts_passed))) {
    stop('!!! Row names of coldata are not in the same order as object cts_passed !!! Aborting...')
}

if(! all(rownames(coldata) %in% colnames(cts_imputed))) {
    stop('!!! Not all row names of coldata are in object cts_imputed !!! Aborting...')
} else if (! all(rownames(coldata) == colnames(cts_imputed))) {
    stop('!!! Row names of coldata are not in the same order as object cts_imputed !!! Aborting...')
}

#begin DESeq2 comparisons
message("Running DESeq2 comparison on raw data...\n")
#head(cts_raw)
coldata
dds <- DESeqDataSetFromMatrix(countData = cts_raw, colData = coldata, design = ~ Treatment + Patient)
dds <- DESeq(dds)
res <- results(dds)
rlt <- DESeq2::rlogTransformation(dds)
res <- results(dds, contrast = comparison)
write.table(counts(dds, normalized=TRUE), 'raw_norm_counts.txt', sep="\t", quote=FALSE)
write.table(res, 'raw.results.txt', sep="\t", quote=FALSE)
pcaplot(rlt, intgroup = c(args[3]), ntop = 1000, pcX = 1, pcY = 2, title = paste("Samples PCA: ", args[1]," vs. ", args[2], " (All Transcripts)", sep = ""), ellipse = TRUE, text_labels=TRUE, point_size = 4)
#pcaplot(rlt, intgroup = c(args[3]), ntop = 1000, pcX = 1, pcY = 2, title = paste("Samples PCA: ", args[1]," vs. ", args[2]," vs. ", args[4], " (All Transcripts)", sep = ""), ellipse = FALSE, text_labels=FALSE, point_size = 4)

message("Running DESeq2 comparison on CLEAR-only data...\n")
dds <- DESeqDataSetFromMatrix(countData = cts_passed, colData = coldata, design = ~ Treatment + Patient)
dds <- DESeq(dds)
res <- results(dds)
rlt <- DESeq2::rlogTransformation(dds)
res <- results(dds, contrast = comparison)
write.table(res, 'clear.results.txt', sep="\t", quote=FALSE)
write.table(counts(dds, normalized=TRUE), 'clear_norm_counts.txt', sep="\t", quote=FALSE)
pcaplot(rlt, intgroup = c(args[3]), ntop = 1000, pcX = 1, pcY = 2, title = paste("Samples PCA: ", args[1]," vs. ", args[2], " (CLEAR-Only Transcripts)", sep = ""), ellipse = TRUE, text_labels=TRUE, point_size = 4)
#pcaplot(rlt, intgroup = c(args[3]), ntop = 1000, pcX = 1, pcY = 2, title = paste("Samples PCA: ", args[1]," vs. ", args[2]," vs. ", args[4], " (CLEAR-Only Transcripts)", sep = ""), ellipse = FALSE, text_labels=FALSE, point_size = 4)


message("Running DESeq2 comparison on imputed data...\n")
dds <- DESeqDataSetFromMatrix(countData = cts_imputed, colData = coldata, design = ~ Treatment + Patient)
dds <- DESeq(dds)
res <- results(dds)
rlt <- DESeq2::rlogTransformation(dds)
res <- results(dds, contrast = comparison)
write.table(res, 'imputed.results.txt', sep="\t", quote=FALSE)
write.table(counts(dds, normalized=TRUE), 'imputed_norm_counts.txt', sep="\t", quote=FALSE)
pcaplot(rlt, intgroup = c(args[3]), ntop = 1000, pcX = 1, pcY = 2, title = paste("Samples PCA: ", args[1]," vs. ", args[2], " (25% kNN-Imputed Transcripts)", sep = ""), ellipse = TRUE, text_labels=TRUE, point_size = 4)
#pcaplot(rlt, intgroup = c(args[3]), ntop = 1000, pcX = 1, pcY = 2, title = paste("Samples PCA: ", args[1]," vs. ", args[2]," vs. ", args[4], " (25% kNN-Imputed Transcripts)", sep = ""), ellipse = FALSE, text_labels=FALSE, point_size = 4)
