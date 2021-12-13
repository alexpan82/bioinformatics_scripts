# Takes a counts table and the 2 columns to plot on the x and y axis respectively
# Example, if counts.tsv contains headers named 'A' and 'B'
# Rscript plotExpression_1x1.R counts.tsv A B [OUTPREFIX]

library('ggplot2')
args <- commandArgs(TRUE)

png(paste(args[4], ".png", sep=""))
counts <- read.table(args[1], header=TRUE, sep='\t')
xy <- c(args[2], args[3])

# Subset counts table
counts <- counts[, xy]
counts <- data.frame(counts)

# Only keep shared genes
row_sub = apply(counts, 1, function(row) all(row != 0 ))
counts <- counts[row_sub,]

# remove 0s and add 1 for log transform
counts <- counts[rowSums(counts[])>0,]
print(dim(counts))

# print correlation
print(paste("Correlation of log10(x) and log10(y): ", cor(x=log10(counts[, xy[1]]), y=log10(counts[, xy[2]]))))
# Plot
colfunc <- c("darkblue","green","yellow","red","black")
breaks <- c(0.01, 0.1, 1, 10, 100, 1000, 10000)
p <- ggplot(counts, aes(x=counts[, xy[1]], y=counts[, xy[2]])) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_gradientn(colors=colfunc) +
  scale_x_log10(limits=c(0.01, 10000), breaks = breaks, labels = scales::trans_format("log10",scales::math_format(.x)), expand=c(0,0)) +
  scale_y_log10(limits=c(0.01, 10000), breaks = breaks, labels = scales::trans_format("log10",scales::math_format(.x)), expand=c(0,0)) +
  xlab(xy[1]) + ylab(xy[2])
plot(p)
