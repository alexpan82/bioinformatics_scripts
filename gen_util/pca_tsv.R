# Plots PCA for any tsv file
# Assumes dependent variables (things to plot) are col
# Independent variables (things to turn into PCs) are rows
# Example usage:
# Rscript pca_tsv.R normal [TSV]

# Can also filter by row mean
# Rscript pca_tsv.R rowFilter [TSV] [rowMeanCutoff]

# Can also filter by a column in a different file (with same row names)
# Rscript pca_tsv.R fileFilter [TSV] [TSV2] [columnOfTSV2] [cutoff]

library(ggfortify)
pdf("pca.pdf")

args <- commandArgs(TRUE)
print(args[2])
data <- as.matrix(read.table(args[2], header=TRUE, sep="\t", row.names=1))
filter <- 0

if (args[1] == "rowFilter") {
	filter <- as.numeric(args[3])
	print(paste("Filtering row avgs by:", filter))
	
	avg_row <- apply(data, 1, mean)
	filtered_rows <- Filter(function(x) x > filter, avg_row)
	#filtered_rows <- lapply(avg_row, function(x) x > filter)
	data <- data[names(filtered_rows),]
	print(nrow(data))

} else if (args[1] == "fileFilter") {
	criteria <- as.matrix(read.table(args[3], header=TRUE, sep="\t", row.names=1))
	filterCol <- as.numeric(args[4])
	filter <- as.numeric(args[5])
	print(paste("Filtering rows in", args[3], ">=", filter, "in column", filterCol, sep=" "))
	
	filtered_rows <- criteria[criteria[, filterCol] >= filter,]
	data <- data[rownames(filtered_rows),]
	summary(data)

} else if (args[1] == "normal") {
	print("PCA without any type of filtering expect removing rows with all 0's")
	# print(ncol(data))
	# print(nrow(data[rowSums(data[])==0,]))
	# non_zero <- apply(data[,-1], 1, function(x) !all(x==0))
	non_zero <- apply(data[,1:ncol(data)], 1, function(x) !all(x==0))
	data <- data[non_zero,]
	# print(nrow(data))

} else {
	print("Please enter 'rowFilter' or 'fileFilter' or 'normal'")
	q()
}
head(data)
data <- t(data)
#data <- log(data[,1:ncol(data)]+1)
pca <- prcomp(data, center = TRUE, scale. = TRUE)
#pca <- prcomp(data, center = TRUE, scale. = FALSE)
summary(pca)
write.table(pca$rotation, 'pca.txt', sep="\t", quote=FALSE)
autoplot(pca, label=FALSE, label.size = 4) + geom_text(vjust=-1, label=rownames(data))
