#!/usr/bin/Rscript

# Written by Altan Turkoglu (turkoglu.12@osu.edu)
# Wednesday, May 11th, 2020

# This must be run using R/3.6.0
# Requires gglot2, scales, gridExtra, reshape2

# argument 1 should be a directory containing: 10X.txt, bulk.txt, lcRNA.txt, scRNA.txt, 1ng.txt

#/users/PAS0472/osu8725/tech_manuscript/10x_data/VAC05-C24-comb2/filtered_feature_bc_matrix.csv
#/users/PAS0472/osu8725/tech_manuscript/clonetech_lcRNA/vac_5_15/output/L5295_AwanF_Sub5-C24D1-300Bcell-1_V1C_1_S5_L008/L5295_AwanF_Sub5-C24D1-300Bcell-1_V1C_1_S5_L008.genes.results

#/fs/scratch/PAS0472/osu9900/tech_manuscript/clonetech_scRNA/vac5_15/output/sub5-C24_comb/sub5-C24_comb.genes.results

#change to show expression and not percentile rank DONE
#log scales DONE
#create UMIs per million (divide by the sum of all UMIs, multiply by 1,000,000) DONE
#r^2 value (like in CLEAR) almost DONE
#incorporate scRNA tpm_bulk, tpm_bulk RNA, 1 ng DONE

#add # genes with 0 coverage on diagonal

#import arguments
args <- commandArgs(TRUE)
filedir <- args[1] #this should be a directory containing at least one of the following: 10X.txt, bulk.txt, lcRNA.txt, scRNA.txt, 1ng.txt
title <- args[2] #title

#filedir <- '/users/PAS0472/osu8725/tech_manuscript/transcriptome_figures/10X/'
#title <- "10X"

require('ggplot2')
require('scales')
require('grid')
require('gridExtra')
require('reshape2')
require('lattice')
require('gridGraphics')
require('VennDiagram')
require('cowplot')

na.zero <- function (x) {
    x[is.na(x)] <- 0
    return(x)
}

grabFile <- function(name) {
	return(paste(filedir,list.files(filedir)[match(name, list.files(filedir))], sep='/'))
}
	
import <- function(name, id, data, tenx=FALSE) { #
	if ( name %in% list.files(filedir)) { message(paste("Detected ", name, "...", sep="")) } else { return(data) }
	
	if (tenx) {
		dat <- as.data.frame(rowSums(as.data.frame(read.table(grabFile(name), header=TRUE, sep=',', row.names=1))))
		dat <- dat/sum(dat)*1000000 #create UMIs per million transform
		colnames(dat) <- id
	}
	else {
		dat <- as.data.frame(read.table(grabFile(name), header=TRUE, sep="\t", row.names=1))[, "TPM", drop=FALSE]
		rownames(dat) <- sapply(strsplit(row.names(dat), split="_"), "[[", 1)
		colnames(dat) <- id
	}
	
	data <- merge(data, dat, by=0, all=TRUE)
	rownames(data) <- data$Row.names
	data$Row.names <- NULL
	
	return(data)
}

parseFiles <- function() {
	message("Importing data...")

	data = data.frame()
	for (i in list.files(filedir)[grepl(".txt",list.files(filedir))]) {
		data <- import(i, gsub(".txt", "", i), data, FALSE)
	}

	data[data == 0] <- NA

	return(data)
}

data <- parseFiles()

#create label vector
axis_labels <- colnames(data)
axis_labels <- paste(axis_labels, "(TPM,")

#filtered = data[data$umi > 0 & data$ng > 0, ]
#nofilter = data[!(data$umi == 0 & data$lc == 0), ]

#plotting all genes
message("Plotting...")

heat_theme <- theme(axis.line=element_blank(),
					axis.text=element_text(face="bold", colour="black"),
					#axis.ticks=element_blank(),
					#axis.title.x=element_blank(),
					#axis.title.y=element_blank(),
					axis.title=element_blank(),
					legend.position="none",
					panel.background=element_blank(),
					#panel.border=element_blank(),
					#panel.grid.major=element_blank(),
					#panel.grid.minor=element_blank(),
					#plot.background=element_blank()
					plot.margin = unit(c(0,0,0.1,0.1), "in"),
					axis.ticks.length = unit(0,"null"),
				)
					
no_theme <- theme(	axis.line=element_blank(),
					axis.text=element_text(face="bold", colour="white"),
					axis.ticks=element_line(colour="white"),
					axis.title=element_blank(),
					#axis.title.x=element_text(colour="white"),
					#axis.title.y=element_text(colour="white"),
					legend.position="none",
					panel.background=element_blank(),
					panel.border=element_rect(colour = 'black', fill = NA, size=1),
					#panel.grid.major=element_blank(),
					#panel.grid.minor=element_blank(),
					#plot.background=element_blank()
					plot.margin = rep(unit(0,"null"),4),
					axis.ticks.length = unit(0,"null"),
				)

text_theme <- theme(panel.background=element_blank(), 
					plot.background=element_blank(),
					#plot.margin = rep(unit(0,"null"),4),
				)


heatplotlog <- function(data, datax, datay) {
	colfunc <- c("darkblue","green","yellow","red","black")
	
	#function for breaks
	base_breaks <- function(n = 10){
		function(x) {
			axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
		}
	}
	
	#remove zeroes
	cases <- complete.cases(datax, datay)
	data <- data[cases, ]
	datax <- datax[cases]
	datay <- datay[cases]
	
	lowlim <- min(c(min(datax), min(datay)))
	highlim <- max(c(max(datax), max(datay)))
		
	return(
		ggplot(data, aes(x=datax, y=datay)) +
		stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
		scale_fill_gradientn(colors=colfunc) +
		scale_x_log10(limits=c(0.01, 12000), breaks = base_breaks(), labels = scales::trans_format("log10", scales::math_format(.x)), expand=c(0,0)) +
		scale_y_log10(limits=c(0.01, 12000), breaks = base_breaks(), labels = scales::trans_format("log10", scales::math_format(.x)), expand=c(0,0)) +
		theme(legend.position = "none") +
		xlab("") +
		ylab("") +
		heat_theme
		+ annotation_logticks(colour="white", sides="trbl", mid=unit(0.1, "cm"), long=unit(0.2, "cm"))
		# + geom_segment(aes(x = lowlim, xend = highlim, y = lowlim, yend = highlim))
	)
}


plotGrid <- function(data) {
	if ( ! is.null(data$umi) ) { data <- data[ ,2:ncol(data)] } #remove non-upm column
	cols = ncol(data)
	graphs = vector('list', (cols+1)*(cols+1))
	count = 1
	
	#transcript dropout metrics
	total_genes <- nrow(data[rowSums(na.zero(data)) != 0, ])

	for (i in 0:cols) {
		for (j in 0:cols+1) {
			if (i > 0 & j > 0 & j <= cols) {
				#extract data
				dat <- data[, c(i,j)]
				datx <- data[, i]
				daty <- data[, j]
				
				#pull only complete cases from data
				cases = complete.cases(dat)
				dat <- dat[cases, ]
				datx <- datx[cases]
				daty <- daty[cases]
			}
			
			graph = ggplot()
			
			if (i == 0) { #top text
				if (j == cols+1) { #top right blank buffer
					graph <- ggplot() + theme_void()
				} else {
					graph <- ggplot() + xlab(paste(axis_labels[j], "log10)")) + text_theme
				}
			} else if (j > cols) { #right text
				graph <- ggplot() + ylab(paste(axis_labels[i], "log10)")) + text_theme + theme()
			} else if (i == j) { #diagonals
				dropout <- total_genes - length(datx)
				anno <- paste(paste(head(unlist(strsplit(axis_labels[i], " ")), -1), collapse=" "), "\nDropout:\n", comma(dropout), " / ", comma(total_genes), " (", percent(dropout/total_genes), ")\nGenes", sep="")
				
				tdat <- rbind(data.frame(),c(1,1, anno)); colnames(tdat) <- c("x","y", "label")
				
				graph <- ggplot(tdat, aes(x=x,y=y, label=anno)) + geom_point(colour="white") + geom_text(aes(label=label)) + no_theme
			} else if (i > j) { #bottom left
				comp <- paste(paste(head(unlist(strsplit(axis_labels[i], " ")), -1), collapse=" "), "\nvs.\n", paste(head(unlist(strsplit(axis_labels[j], " ")), -1), collapse=" "), sep="")
				stat <- format(cor.test(log10(datx), log10(daty), method="pearson")$estimate[[1]], digits=2) #log transformed
				anno <- comp
				
				tdat <- rbind(data.frame(),c(1,1, anno)); colnames(tdat) <- c("x","y", "label")
				
				graph <- ggplot(tdat, aes(x=x,y=y, label=anno)) + geom_point(colour="white") + geom_text(aes(label=label)) + no_theme

			} else if (i < j) { #heatmaps
				stat <- format(cor.test(log10(datx), log10(daty), method="pearson")$estimate[[1]], digits=2) #log transformed
				anno1 <- data.frame(xpos=0.03162, ypos=3162, annotateText=c(paste("PCC:",stat)), hjustvar=0, vjustvar=1)
				
				tran <- comma(length(datx))
				anno2 <- data.frame(xpos=3162, ypos=0.03162, annotateText=c(paste("Genes:",tran)), hjustvar=1, vjustvar=0)
				
				graph <- heatplotlog(dat, datx, daty) + heat_theme + 
				geom_text(data=anno1,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), color="white") +
				geom_text(data=anno2,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), color="white")
				#if (j==1) { graph <- graph + ylab(paste(axis_labels[i], "(log10)")) } else { graph <- graph + theme(axis.title.y=element_blank())}
				#if (i==cols) { graph <- graph + xlab(paste(axis_labels[j], "(log10)")) } else { graph <- graph + theme(axis.title.x=element_blank())}
			}
			
			graphs[[count]] <- graph
			count = count + 1
		}
	}

	return(graphs)
}

if ( ! is.null(data$umi) ) { cols <- ncol(data) } else { cols <- ncol(data) + 1}

graphs <- plotGrid(data)

png(paste(title,".heatmaps.png",sep=""), width = 12, height=12, units = "in", res = 600)
do.call(grid.arrange, c(graphs, ncol=cols))
dev.off()

#create aggregate of 1ng and bulk data
agg_bulk <- function(data) {
	agg <- na.zero(data$ng) + na.zero(data$bulk)
	rownames(agg) <- rownames(data)

	return(agg)
}
	
plotGenes <- function(data) {
	if ( ! is.null(data$umi) ) { data <- data[ ,2:ncol(data)] } #remove non-upm column
	cols = ncol(data)

	data <- data[complete.cases(data), ]
	
	P <- apply(data, 2, ecdf)
	#calculate gene ranks
	#umiP <- ecdf(data$umi)
	#scP <- ecdf(data$sc)
	#lcP <- ecdf(data$lc)
	#ngP <- ecdf(data$ng)
	#bulkP <- ecdf(data$bulk)

	#data$umiP <- umiP(data$umi)
	#data$scP <- scP(data$sc)
	#data$lcP <- lcP(data$lc)
	#data$ngP <- SCp(data$ng)
	#data$bulkP <- bulkP(data$bulk)
}
#dropout

if (FALSE) {
expressionVenn <- function(data) {
	if ( ! is.null(data$umi) ) { data <- data[ ,2:ncol(data)] }
	cols = ncol(data)
	genes = vector('list', cols)
	
	#create label vector
	axis_labels <- colnames(data)
	axis_labels <- gsub("upm", "10X", axis_labels)
	axis_labels <- gsub("sc", "scRNA", axis_labels)
	axis_labels <- gsub("lc", "lcRNA", axis_labels)
	axis_labels <- gsub("ng", "1ng", axis_labels)
	axis_labels <- gsub("bulk", "Bulk", axis_labels)	

	fill <- c("red", "green", "yellow", "blue","orange")
	alpha <- c(0.5,0.5,0.5,0.5,0.5)
	fill <- fill[1:cols]
	alpha <- alpha[1:cols]

	for (i in 1:cols) {
		cases = complete.cases(data[ , i])
		genes[[i]] <- row.names(data[cases, ][i])
	}
	
	venn.plot <- venn.diagram(genes, NULL, fill=fill, alpha=alpha, cex = 2, cat.fontface=4, category.names=axis_labels, main=paste("Gene Expression Overlap:", title))
	
	grid.newpage()
	grid.draw(venn.plot)
}

dropoutTable <- function(data) {
	if ( ! is.null(data$umi) ) { data <- data[ ,2:ncol(data)] }
	cols = ncol(data)
	genes = vector('list', cols)

	dropout <- data.frame(row.names=colnames(data))
	
}

#unused plots #####################################################################################################


heatplot <- function(data, datax, datay) {
	colfunc <- c("darkblue","green","yellow","red","black")
	
	return(ggplot(data, aes(x=datax, y=datay)) +
		stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
		scale_fill_gradientn(colors=colfunc) +
		scale_x_continuous(expand = c(0, 0)) +
		scale_y_continuous(expand = c(0, 0)) +
		theme(legend.position = "none") +
#		ggtitle(title) +
		xlab("UMIs per Million") +
		ylab("TPM (300 Cells)") +
		geom_abline())
}
heatplotperc <- function(data) {
	colfunc <- c("darkblue","green","yellow","red","black")
	
	return(ggplot(data, aes(x=data$umiP, y=data$lcP)) +
		stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
		scale_fill_gradientn(colors=colfunc) +
		scale_x_continuous(expand = c(0, 0)) +
		scale_y_continuous(expand = c(0, 0)) +
		theme(legend.position = "none") +
#		ggtitle(title) +
		xlab("UMI Percentile Rank") +
		ylab("TPM Percentile Rank") +
		geom_abline())
}
scatterplot <- function(data) {
	return(ggplot(data, ) +
		geom_point(aes(x=upm, y=lc)) + 
#		ggtitle(title) +
		xlab("UMIs per Million") +
		ylab("TPM (300 Cells)") +
		theme(plot.title = element_text(hjust = 0.5), legend.position='none'))
		#scale_colour_manual(values = setNames(c('red','black'),c(T, F))))
}
scatterplotlog <- function(data) {
	return(ggplot(data, ) +
		geom_point(aes(x=upm, y=lc)) + 
		scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x))) +
		scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x))) +		
#		ggtitle(title) +
		xlab("UMIs per Million") +
		ylab("TPM (300 Cells)") +
		theme(plot.title = element_text(hjust = 0.5), legend.position='none'))
		#scale_colour_manual(values = setNames(c('red','black'),c(T, F))))
}

heatplotlogscale <- function(data, datax, datay) {
	#function for breaks
	base_breaks <- function(n = 10){
		function(x) {
			axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
		}
	}

	#remove zeroes
	cases <- complete.cases(datax, datay)
	data <- data[cases, ]
	datax <- datax[cases]
	datay <- datay[cases]
	
	lowlim <- min(c(min(datax), min(datay)))
	highlim <- max(c(max(datax), max(datay)))

	return(
		ggplot(data, aes(x=datax, y=datay)) +
		geom_blank() +
		scale_x_log10(limits=c(min(datax), max(datax)), breaks = base_breaks(), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
		scale_y_log10(limits=c(min(datay), max(datay)), breaks = base_breaks(), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
		theme(legend.position = "none") +
		xlab("") +
		ylab("") +
		#annotation_logticks() +
		theme(
			panel.background=element_blank(),
			axis.ticks.margin=unit(0.2, "cm"),
			axis.ticks.length=unit(-0.1, "cm"),
			axis.line.x = element_line(color="black"),
			axis.line.y = element_line(color="black"),
			plot.margin = unit(c(0,0,0,0), "cm")
			#axis.text.y = element_text(margin=margin(1,-1,1,0, unit="cm")),
			#axis.text.x = element_text(margin=margin(-1,0,0,1, unit="cm"))
		)
		
	)
}
}