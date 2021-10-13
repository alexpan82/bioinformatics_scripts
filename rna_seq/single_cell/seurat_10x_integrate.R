# Performs DEG analysis across scRNA data from 10x (cellranger output)
# Needs path to 10x cellranger filtered_feature_bc_matrix/ directory for both samples
# Please give adequate resources to this script
	# Min 12 cores and x hours
# Need to install Seurat
	# To avoid any compling/installation issues, create a conda environment with R
		# create or load your conda environment and type: conda install -c r r
	# Then systematically start downloading your dependancies
		# conda install -c bioconda r-seurat
	# If 'unable to load shared object' errors occur, just delete the pkg that R tries to reference (thats outside the conda env)

# Usage: Rscript seurat_10x_integrate.R [filtered_feature_bc_matrix1/] [filtered_feature_bc_matrix2/] [NAME1] [NAME2]

library(Seurat)
library(ggplot2)
library(cowplot)

pdf("output.pdf")

args = commandArgs(trailingOnly=TRUE)

# Read in arguments
countFile1 <- Read10X(data.dir=args[1])
countFile2 <- Read10X(data.dir=args[2])
name1 <- args[3]
name2 <- args[4]

# Make Seurat object and give conditions to cells
	# Only include genes expressed in at least 5 cells
	# Only include cells with at least 200 genes

sobject1 <- CreateSeuratObject(counts = countFile1, min.cells = 5, min.features=200)
sobject2 <- CreateSeuratObject(counts = countFile2, min.cells = 5, min.features=200)
sobject1
sobject2
name1list <- rep(name1 ,length(Idents(sobject1)))
name2list <- rep(name2 ,length(Idents(sobject2)))

# Add information about groups to seurat objects
sobject1 <- AddMetaData(object = sobject1, metadata = name1list, col.name = 'groupname')
sobject2 <- AddMetaData(object = sobject2, metadata = name2list, col.name = 'groupname')


# QC for each object
	# mitrochondrial gene
sobject1[["percent.mt"]] <- PercentageFeatureSet(sobject1, pattern = "^MT-")
sobject2[["percent.mt"]] <- PercentageFeatureSet(sobject2, pattern = "^MT-")
	# We filter cells that have unique feature counts over 2,500 or less than 200
	# We filter cells that have >5% mitochondrial counts
sobject1 <- subset(sobject1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
sobject2 <- subset(sobject2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
sobject1
sobject2

# Normalize
sobject1 <- NormalizeData(sobject1, normalization.method = "LogNormalize", scale.factor = 10000)
sobject2 <- NormalizeData(sobject2, normalization.method = "LogNormalize", scale.factor = 10000)

# Find highly variable features
sobject1 <- FindVariableFeatures(sobject1, selection.method = "vst", nfeatures = 2000)
sobject2 <- FindVariableFeatures(sobject2, selection.method = "vst", nfeatures = 2000)

# Combine the two objects
immune.anchors <- FindIntegrationAnchors(object.list = list(sobject1, sobject2), dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)

DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "groupname")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
p1
p2

DimPlot(immune.combined, reduction = "umap", split.by = "groupname")

# Identify conserved cell type markers
DefaultAssay(immune.combined) <- "RNA"

markers <- FindConservedMarkers(immune.combined, ident.1 = 1, grouping.var = "groupname", verbose = FALSE)

write.table(as.matrix(markers), file = "conserved_markers.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

# Feature plot
FeaturePlot(immune.combined, features = c("CD5", "CD19"), min.cutoff = "q9", blend=TRUE)
diffmarkers <- FindAllMarkers(immune.combined, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(as.matrix(diffmarkers), file = "diff_markers.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

