# Takes DESeq2 results and makes volcano plot
library(ggrepel)

pdf("volcano.pdf")
args <- commandArgs(TRUE)
# input_name
input_name <- args[1]
deseq2_table <- as.data.frame(read.table(input_name,
                                        header=TRUE, sep="\t", row.names=NULL))

# Make row names unique and remove NAs
rownames(deseq2_table) <- make.names(deseq2_table[, 1], unique = TRUE)
deseq2_table <- deseq2_table[complete.cases(deseq2_table), ]
# Rename 1st col
colnames(deseq2_table)[1] <- "gene_symbol"

# Add col of significant DEGs
deseq2_table$diffexpressed <- "NO"
deseq2_table$diffexpressed[deseq2_table$log2FoldChange > 1 & deseq2_table$padj < 0.05] <- "UP"
deseq2_table$diffexpressed[deseq2_table$log2FoldChange < -1 & deseq2_table$padj < 0.05] <- "DOWN"
deseq2_table$delabel <- NA
deseq2_table$delabel[deseq2_table$diffexpressed != "NO"] <- deseq2_table$gene_symbol[deseq2_table$diffexpressed != "NO"]
# write.table(deseq2_table, file = "test.txt", quote = FALSE, sep = "\t")
# Plot
# ggplot(data=deseq2_table, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
ggplot(data=deseq2_table, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed)) +
  geom_point() + 
  theme_minimal() +
  # geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylim(-5, 100)
print(range(deseq2_table$padj))
