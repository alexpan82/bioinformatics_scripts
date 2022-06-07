#!/usr/bin/Rscript

# Written by Altan Turkoglu (turkoglu.12@osu.edu)
# Wednesday, May 11th, 2020

# This must be run using R/3.6.0
# Requires gglot2, scales, gridExtra, reshape2

# argument 1 should be a directory containing: 10X.txt, bulk.txt, lcRNA.txt, scRNA.txt, 1ng.txt

#/users/PAS0472/osu8725/tech_manuscript/10x_data/VAC05-C24-comb2/filtered_feature_bc_matrix.csv
#/users/PAS0472/osu8725/tech_manuscript/clonetech_lcRNA/vac_5_15/output/L5295_AwanF_Sub5-C24D1-300Bcell-1_V1C_1_S5_L008/L5295_AwanF_Sub5-C24D1-300Bcell-1_V1C_1_S5_L008.genes.results

#/fs/scratch/PAS0472/osu9900/tech_manuscript/clonetech_scRNA/vac5_15/output/sub5-C24_comb/sub5-C24_comb.genes.results

#add # genes with 0 coverage on diagonal

#import arguments
args <- commandArgs(TRUE)
filedir <- args[1] #this should be a directory containing at least one of the following: 10X.txt, bulk.txt, lcRNA.txt, scRNA.txt, 1ng.txt
title <- args[2] #title

#filedir <- '/users/PAS0472/osu8725/tech_manuscript/transcriptome_figures/CLL1/'
#title <- "CLL1"

#file_name <- "/users/PAS0472/osu8725/tech_manuscript/transcriptome_figures/CLL1/lcRNA.bed"

require('ggplot2')
require('scales')
require('grid')
require('ggplotify')
require('gridExtra')
require('reshape2')
require('lattice')
require('gridGraphics')
require('VennDiagram')

library('reticulate')
use_python("~osu8725/anaconda2/bin/python")
source_python("/users/PAS0472/osu8725/tools/tech_manuscript/replaceGencode.py")
source_python("/users/PAS0472/osu8725/tools/tech_manuscript/gene-coverage-plotter.py")


na.zero <- function (x) {
    x[is.na(x)] <- 0
    return(x)
}
	
grabFile <- function(name) {
	return(paste(filedir,list.files(filedir)[match(name, list.files(filedir))], sep='/'))
}
	
import <- function(name, id, data, tenx=FALSE) { #
	if ( name %in% list.files(filedir)) { message(paste(title, ": Detected ", name, "...", sep="")) } else { return(data) }
	
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
	message(paste(title,": Importing data...", sep=""))

	data = data.frame()

	data <- import("10X.txt", "upm", data, TRUE)
	data <- import("lcRNA.txt", "lc", data)
	data <- import("1ng.txt", "ng", data)
	data <- import("bulk.txt", "bulk", data)

	data[data == 0] <- NA

	return(data)
}

data <- parseFiles()

#labeling
axis_labels <- colnames(data)
axis_labels <- gsub("upm", "10X-5' (UPM,", axis_labels)
axis_labels <- gsub("sc", "scRNA (TPM,", axis_labels)
axis_labels <- gsub("lc", "300-LC (TPM,", axis_labels)
axis_labels <- gsub("ng", "1ng-LC (TPM,", axis_labels)
axis_labels <- gsub("bulk", "BULK (TPM,", axis_labels)

return_label <- function(lab) {
	lab <- gsub("upm", "10X-5' (UPM,", lab)
	lab <- gsub("sc", "scRNA (TPM,", lab)
	lab <- gsub("lc", "300-LC (TPM,", lab)
	lab <- gsub("ng", "1ng-LC (TPM,", lab)
	lab <- gsub("bulk", "BULK (TPM,", lab)
	
	return(lab)
} 	

sublabel <- function(string) {
	if (string == "upm") {
		return("10X-5'")
	} else if (string == "lc") {
		return("300-LC")
	} else if (string == "ng") {
		return("1ng-LC")
	} else if (string == "bulk") {
		return("BULK")
	}
}

subunit <- function(string) {
	if (string == "upm") { return("UPM") } 
	else { return("TPM") }
}

#filtered = data[data$umi > 0 & data$ng > 0, ]
#nofilter = data[!(data$umi == 0 & data$lc == 0), ]

#plotting all genes
message(paste(title, ": Plotting...", sep=""))

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
				
heat_theme_labeled <- theme(axis.line=element_blank(),
					axis.text=element_text(face="bold", colour="black"),
					#axis.ticks=element_blank(),
					#axis.title.x=element_blank(),
					#axis.title.y=element_blank(),
					axis.title=element_text(),
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

assignlevels <- function(data) {
	comb <- rowMeans(na.zero(data)[ , , drop=FALSE])
	comb[comb == 0] <- NA

	data$level <- vector("character", length(row.names(data)))
	
	quants <- quantile(comb, probs=seq(0,1,1/3), na.rm = TRUE)
	
	
	message(paste(title, ": Assigning gene levels", sep=""))
	for (gene in row.names(data)) {
		if (is.na(comb[gene])) {
			data[gene, ]$level <- "NA"
		}
		else if (comb[gene] >= quants[1] & comb[gene] <= quants[2]) {
			data[gene, ]$level <- "Low"
		}
		else if (comb[gene] > quants[2] & comb[gene] <= quants[3]) {
			data[gene, ]$level <- "Medium"
		}
		else if (comb[gene] > quants[3] & comb[gene] <= quants[4]) {
			data[gene, ]$level <- "High"
		}
	}
	
	return(data)
}

barplotdropout <- function(data, LABELS) {

	message(paste(title,": Creating dropouts matrix"), sep="")
	#matrix of dropouts
	dropout <- data.frame(matrix(ncol=4, nrow=1))
	dropoutlabels <- unlist(lapply(lapply(strsplit(axis_labels, " "), head, -1), paste, collapse=" "))
	colnames(dropout) <- c("Sample", "Level", "Dropouts", "Total")

	for (sample in colnames(data)[1:length(dropoutlabels)]) {
		for (level in c("Low", "Medium", "High")) {
			drops <- length(rownames(data[is.na(data[, sample]) & data$level == level, ]))
			total <- length(rownames(data[data$level == level, ]))
			dropout <- rbind(dropout, c(sample, level, drops, total))
		}
	}
	dropout <- dropout[2:length(dropout$Sample), ] #tidy rownames 
	rownames(dropout) <- 1:length(dropout$Sample)
	dropout$Dropouts <- as.numeric(dropout$Dropouts)
	dropout$Total <- as.numeric(dropout$Total)
	
	y <- length(data[data$level == "Low", "level"])
	
	plot <- ggplot(dropout, aes(x=factor(Sample, levels=colnames(data)[1:length(dropoutlabels)]), y=Dropouts)) + 
			geom_bar(aes(fill=factor(Level, levels=c("Low", "Medium", "High"))), position = "dodge", stat="identity") +
			scale_x_discrete(labels=dropoutlabels) +
			scale_y_continuous(label=comma, expand = c(0, 0), limits=c(0, y*1.05)) +
			labs(fill="Expression Level") +
			ggtitle(paste("Gene Dropout by Sample Type and Expression Level:", title)) +
			xlab("Sample Type") +
			ylab("Dropout Genes") +
			geom_hline(yintercept=y, linetype="dashed", color = "red")
			
	if (LABELS) {
		plot <- plot + 
				geom_text(aes(label=Dropouts), position=position_dodge2(width=0.9), vjust=-0.25) + 
				annotate("text", x=0, y=y, label=as.character(y), vjust=-0.5, hjust=-.25)
	}
	plot
	
	return(plot)
}

plotGeneCovs <- function(data) {
	#reorder our data by row means
	data <- data[ order(rowMeans(na.zero(data[ , 1:length(axis_labels)])), decreasing=T), ]
	images <- matrix("", nrow=4, ncol=length(axis_labels))
	
	i=1
	for (value in c(1, 10, 100, 1000)) {
		#midexp <- median(data[data$level == lev & ! is.na(data[samp]), samp])
		closest <- which.min(abs(value-rowMeans(na.zero(data[ , 1:length(axis_labels)])))) #find closest row to value
		next_closest <- closest
		
		blocklist = c("TMEM63C", "ARL13B") #genes that don't look good.
		
		#test our gene
		test <- data[closest, 1:length(axis_labels)]
		flip <- FALSE
		inc <- 1
		while (any(is.na(test)) | max(test)/min(test) > 2 | startsWith(as.character(mapping[rownames(test)]), "RP")  | toupper(mapping[rownames(test)]) %in% blocklist | next_closest <= 0 ) { #do not allow for NAs OR l2fc > 2 OR RPL OR RPS
			if (flip) {
				next_closest <- closest - inc
				test <- data[next_closest, 1:length(axis_labels)]
				inc <- inc + 1
				flip <- !flip
			}
			else {
				next_closest <- closest + inc
				test <- data[next_closest, 1:length(axis_labels)]
				flip <- !flip
			}
		}
		
		gene <- mapping[rownames(test)]
		
		j=1
		for (samp in colnames(data)[1:length(axis_labels)]) {
			file_name <- ""
			if (samp == "upm") {
				file_name <- "10X.bed"
			} else if (samp == "lc") {
				file_name <- "lcRNA.bed"
			} else if (samp == "ng") {
				file_name <- "1ng.bed"
			} else if (samp == "bulk") {
				file_name <- "bulk.bed"
			}
			
			message(paste(title, ": ", paste("Plotting", as.character(gene),"for", file_name, "at", value, "TPM..."), sep=""))
			
			tpm <- data[next_closest, samp]
			tpm <- comma(tpm, accuracy=0.01)
			suffix <- paste("(", tpm, " ", subunit(samp), ")", sep="")
			plot_cov(grabFile(file_name), as.character(gene), suffix, paste(samp, value, sep="."))
			
			images[i,j] <- paste("./covplots/",paste(samp,value,toupper(as.character(gene)),"covprofile.png", sep="."),sep="")
			
			j=j+1
			
		}
		
		i=i+1
	}
	
	merge_images(images, title)
}

dropoutVenn <- function(data, lev="all") {
	cols = length(axis_labels)
	genes = vector('list', cols)

	fill <- c("red", "green", "yellow", "blue","orange")
	alpha <- c(0.5,0.5,0.5,0.5,0.5)
	fill <- fill[1:cols]
	alpha <- alpha[1:cols]
	
	dropoutlabels <- unlist(lapply(lapply(strsplit(axis_labels, " "), head, -1), paste, collapse=" "))

	for (i in 1:cols) {
		cases = !complete.cases(data[ , i]) # vector of dropouts
		if (lev =="all") {
			genes[[i]] <- row.names(data[cases, ][i])
		}
		else {
			genes[[i]] <- row.names(data[cases & data$level==lev, ][i])
		}
	}
	if (lev == "all") {
		venn.plot <- venn.diagram(genes, NULL, fill=fill, alpha=alpha, cex = 2, cat.fontface=4, category.names=dropoutlabels, main=paste(title, "Gene Dropout"))
	}
	else {
		venn.plot <- venn.diagram(genes, NULL, fill=fill, alpha=alpha, cex = 2, cat.fontface=4, category.names=dropoutlabels, main=paste(title, "Gene Dropout:", lev, "Expressors"))
	}
	
	grid.newpage()
	grid.draw(venn.plot)

}
#draw single plot
plotSingle <- function(data, str1, str2) {

	if (! (str1 %in% colnames(data))) {
		return(paste(title, ": Data for ", str1, " not found.", sep=""))
	} else if (! (str2 %in% colnames(data))) {
		return(paste(title, ": Data for ", str2, " not found.", sep=""))
	}
	
	#extract data
	dat <- data[, c(str1,str2)]
	datx <- data[, str1]
	daty <- data[, str2]
	
	#pull only complete cases from data
	cases = complete.cases(dat)
	dat <- dat[cases, ]
	datx <- datx[cases]
	daty <- daty[cases]

	stat <- format(cor.test(log10(datx), log10(daty), method="pearson")$estimate[[1]], digits=2) #log transformed
	anno1 <- data.frame(xpos=0.03162, ypos=3162, annotateText=c(paste("PCC:",stat)), hjustvar=0, vjustvar=1)
				
	tran <- comma(length(datx))
	anno2 <- data.frame(xpos=3162, ypos=0.03162, annotateText=c(paste("Genes:",tran)), hjustvar=1, vjustvar=0)
	
	graph <- heatplotlog(dat, datx, daty) + heat_theme_labeled +
				xlab(paste(return_label(str1), "log10)")) + 
				ylab(paste(return_label(str2), "log10)")) +  
				geom_text(data=anno1,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), color="white") +
				geom_text(data=anno2,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), color="white")
	
	return(graph)
}

#draw grid
plotGrid <- function(data) {
	if ( ! is.null(data$umi) ) { data <- data[ ,2:ncol(data)] } #remove non-upm column
	if ( ! is.null(data$level) ) { data$level <- NULL } #remove levels

	cols = ncol(data)
	graphs = vector('list', (cols)*(cols))
	count = 1
	
	#transcript dropout metrics
	total_genes <- nrow(data[rowSums(na.zero(data)) != 0, ])

	for (i in 2:(cols+1)) {
		for (j in 0:(cols-1)) {
			if (i > 0 & j > 0 & i <= cols & j <= cols) {
				#extract data
				dat <- data[, c(j,i)]
				datx <- data[, j]
				daty <- data[, i]
				
				#pull only complete cases from data
				cases = complete.cases(dat)
				dat <- dat[cases, ]
				datx <- datx[cases]
				daty <- daty[cases]
			}
			
			graph = ggplot()
			type = "blank"

			if (j == 0) { #left text
				if (i  == cols+1) { #bot left blank buffer
					type = "botleft"
					graph <- ggplot() + theme_void()
				} else {
					type = "left text"
					graph <- ggplot() + ylab(paste(axis_labels[i], "log10)")) + text_theme + scale_y_continuous(position='right') + theme(axis.title.y.right = element_text(angle=90))
				}
			} else if (i > cols) { #bot text
				type = "bot text"
				graph <- ggplot() + xlab(paste(axis_labels[j], "log10)")) + text_theme + scale_x_continuous(position = 'top')
			} else if (i > j) { #bottom left heatmaps
				type = "heatmap"
				stat <- format(cor.test(log10(datx), log10(daty), method="pearson")$estimate[[1]], digits=2) #log transformed
				anno1 <- data.frame(xpos=0.03162, ypos=3162, annotateText=c(paste("PCC:",stat)), hjustvar=0, vjustvar=1)
				
				tran <- comma(length(datx))
				anno2 <- data.frame(xpos=3162, ypos=0.03162, annotateText=c(paste("Genes:",tran)), hjustvar=1, vjustvar=0)
				
				graph <- heatplotlog(dat, datx, daty) + heat_theme + 
				geom_text(data=anno1,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), color="white") +
				geom_text(data=anno2,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), color="white")
				#if (j==1) { graph <- graph + ylab(paste(axis_labels[i], "(log10)")) } else { graph <- graph + theme(axis.title.y=element_blank())}
				#if (i==cols) { graph <- graph + xlab(paste(axis_labels[j], "(log10)")) } else { graph <- graph + theme(axis.title.x=element_blank())}
			} else { #other
				type = "other"
				graph <- ggplot() + theme_void()
			}
			
			#print(paste(i,j,type, count))
			graphs[[count]] <- graph
			count = count + 1
		}
	}

	return(graphs)
}


#heatmaps
graphs <- plotGrid(data)

message(paste(title, ": Drawing heatmaps...", sep=""))
system("if [ ! -d ./heatmaps/ ]; then mkdir ./heatmaps/; fi")
png(paste("./heatmaps/", gsub(" ", "_", title),".heatmaps.png",sep=""), width = 7.5, height=7.5, units = "in", res = 600)
do.call(grid.arrange, c(graphs, ncol=length(axis_labels)))
dev.off()

png(paste("./heatmaps/", gsub(" ", "_", title),".300_vs_10X.heatmaps.png",sep=""), width = 2.5, height=2.5, units = "in", res = 600)
plotSingle(data, "lc", "upm")
dev.off()

png(paste("./heatmaps/", gsub(" ", "_", title),".bulk_vs_ng.heatmaps.png",sep=""), width = 2.5, height=2.5, units = "in", res = 600)
plotSingle(data, "bulk", "ng")
dev.off()

#assign levels, dropout
data <- assignlevels(data)

message(paste(title, ": Drawing barcharts...", sep=""))
system("if [ ! -d ./barcharts/ ]; then mkdir ./barcharts/; fi")
png(paste("./barcharts/", gsub(" ", "_", title),".dropout.png",sep=""), width = 6, height=6, units = "in", res = 600)
barplotdropout(data, LABELS=FALSE)
dev.off()

png(paste("./barcharts/", gsub(" ", "_", title),".dropout.labeled.png",sep=""), width = 6, height=6, units = "in", res = 600)
barplotdropout(data, LABELS=TRUE)
dev.off()

#message(paste(title, ": Drawing Venn Diagrams...", sep=""))
#system("if [ ! -d ./venns/ ]; then mkdir ./venns/; fi")
#png(paste("./venns/", gsub(" ", "_", title),".doVenn.all.png",sep=""), width = 6, height=6, units = "in", res = 600)
#dropoutVenn(data)
#dev.off()

#png(paste("./venns/", gsub(" ", "_", title),".doVenn.low.png",sep=""), width = 6, height=6, units = "in", res = 600)
#dropoutVenn(data, "Low")
#dev.off()

#png(paste("./venns/", gsub(" ", "_", title),".doVenn.med.png",sep=""), width = 6, height=6, units = "in", res = 600)
#dropoutVenn(data, "Medium")
#dev.off()

#png(paste("./venns/", gsub(" ", "_", title),".doVenn.high.png",sep=""), width = 6, height=6, units = "in", res = 600)
#dropoutVenn(data, "High")
#dev.off()
#system("rm VennDiagram*.log")

message(paste(gsub(" ", "_", title), ": Choosing and plotting gene coverages...", sep=""))
plotGeneCovs(data)

message(paste(gsub(" ", "_", title), ": Done.", sep=""))
quit()


#unused plots #####################################################################################################

if (FALSE) {

#create aggregate of 1ng and bulk data
agg_bulk <- function(data) {
	agg <- na.zero(data$ng) + na.zero(data$bulk)
	rownames(agg) <- rownames(data)

	return(agg)
}

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