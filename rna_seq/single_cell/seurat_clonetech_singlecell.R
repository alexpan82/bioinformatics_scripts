# Performs DEG analysis on scRNA data from a counts table
# Argument 1: counts table in tsv format with headers
# Argument 2: A tsv that specifies the columns in counts table to analyze
	# 1st column must contain counts table column names (make sure names are in same order relative to counts table)
	# 2nd column contains group id

library(Seurat)
library(ggplot2)
library(cowplot)

#pdf("output.pdf")

args = commandArgs(trailingOnly=TRUE)

# Read in arguments
countFile <- read.table(args[1], header=TRUE, sep='\t', row.names=1)
groups <- read.table(args[2], sep='\t', row.names=1)

# Reduce counts table to only the samples of interest
countFile <- countFile[, rownames(groups)]

# Make Seurat object and give conditions to cells
	# Only include genes expressed in at least 5 cells
	# Only include cells with at least 200 genes

sobject <- CreateSeuratObject(counts = countFile, min.cells = 5, min.features=200)
sobject <- AddMetaData(object = sobject, metadata = groups[,1], col.name = 'groupname')
sobject

# QC
# mitrochondrial gene
sobject[["percent.mt"]] <- PercentageFeatureSet(sobject, pattern = "^MT-")
	# Visualize QC metrics as a violin plot
# sobject@meta.data
# VlnPlot(sobject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
	# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
	# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
FeatureScatter(sobject, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(sobject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
	# We filter cells that have unique feature counts over 10,000 or less than 200
	# We filter cells that have >5% mitochondrial counts
sobject <- subset(sobject, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 5)
sobject

# Normalize
sobject <- NormalizeData(sobject, normalization.method = "LogNormalize", scale.factor = 10000)


# Find highly variable features
sobject <- FindVariableFeatures(sobject, selection.method = "vst", nfeatures = 2000)
	# ID 10 most highly variable genes
top10 <- head(VariableFeatures(sobject), 10)
	# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sobject)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2

# Scales and centers genes in the dataset (idk exactly what this means...)
all.genes <- rownames(sobject)
sobject <- ScaleData(sobject, features = all.genes)

# PCA
numcells <- length(groups[,1])
if(numcells > 51){
pc <- 50
} else {
pc <- numcells-1
}
print(pc)
sobject <- RunPCA(sobject, features = VariableFeatures(object = sobject), npcs=pc)
#print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(sobject, dims = 1:2, reduction = "pca")
DimPlot(sobject, reduction = "pca")
DimHeatmap(sobject, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(sobject)

# Cluster
sobject <- RunUMAP(sobject, reduction = "pca", dims = 1:30)
sobject <- FindNeighbors(sobject, reduction = "pca", dims = 1:30)
tobject <- FindClusters(sobject, resolution = 0.5)
#sobject <- FindNeighbors(sobject, dims = 1:13)
#sobject <- FindClusters(sobject, resolution = 0.5)
#Idents(sobject)
#write.table(as.matrix(sobject[["RNA"]]@data), file = "seurat_lognorm_counts.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

# tSNE
#tobject <- RunTSNE(object = sobject, reduction.type = "pca", dims.use = 1:13, perplexity = 28)
DimPlot(tobject, reduction = "umap")
DimPlot(tobject, reduction="umap", group.by="groupname")
#FeaturePlot(tobject, features = c('CD5', 'CD19', 'BTK', 'BIRC3'))



# Report markers of every cluster
tobject.markers <- FindAllMarkers(tobject, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25)
write.table(tobject.markers, file = "seurat_markers.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
if (FALSE) {
#sobject
#sobject@ident


# Get DEGs 
DEGs <- FindMarkers(sobject, ident.1 = 'Top_Clone', ident.2 = 'W/O_Top_Clone')
write.table(DEGs, file = "seurat_DEGS.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
}


