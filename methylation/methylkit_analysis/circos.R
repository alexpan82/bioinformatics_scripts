### Graphs 3 concentric rings. If you want more rings, add more inputs and code accordingly
### Make sure inputs are in BED format with start and stop locations being a range of numbers
	### For example, a DMC list will only plot if you make stop locations differ from start by 1 bp instead of 0 bp. 
 #   ______________________________________________________________________
    #   <RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO>
    #
    #   This demo draw chromosome ideogram with padding between chromosomes, 
    #   highlights, chromosome names, and histogram. 
    #
    #   Usage:
    #
    #   library(RCircos);
    #   demo("RCircos.Polygon.Demo");
    #   ______________________________________________________________________
    #   <RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO>



    #   Load RCircos package and defined parameters
    #   ________________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

library(RCircos);


    #   Load human cytoband data and gene expression data
    #   ________________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

input <- read.table("Amal_3x3_WTvCF_diffMeth_10p.bed");
input1 <- read.table("Amal_3x3_CFvCFEGCG_diffMeth_10p.bed");
input2 <- read.table("test1.txt");
data(UCSC.Mouse.GRCm38.CytoBandIdeogram);
cyto.info <- UCSC.Mouse.GRCm38.CytoBandIdeogram;
tracks.inside <- 3
tracks.outside <- 0


    #   Setup RCircos core components:
    #   _________________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#RCircos.Set.Core.Components(cyto.info, NULL, 10, 10);
RCircos.Set.Core.Components(cyto.info, chr.exclude = NULL, tracks.inside, tracks.outside)
rcircos.params <- RCircos.Get.Plot.Parameters();
rcircos.params$base.per.unit <- 3000;
rcircos.params$track.height <- 0.3; 
rcircos.params$hist.width <- 50;
RCircos.Reset.Plot.Parameters(rcircos.params);
RCircos.List.Plot.Parameters();


    #   Open the graphic device (here a pdf file)
    #   _________________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

pdf(file="RCircosDemoMouseGenome.pdf", height=8, width=8);
RCircos.Set.Plot.Area();


    #   Draw chromosome ideogram
    #   _________________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

message("Draw chromosome ideogram ...\n");

RCircos.Chromosome.Ideogram.Plot();
title("RCircos Polygon Plot Demo");


    #   Plot histogram Inside of chromosome ideogram
    #   _________________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

polygon.data <- input;
polygon1.data <- input1;
polygon2.data <- input2;

plot.colors <- rep("red", nrow(polygon.data))
plot1.colors <- rep("blue", nrow(polygon1.data))
plot2.colors <- rep("green", nrow(polygon2.data))
polygon.data["PlotColor"] <- plot.colors
polygon1.data["PlotColor"] <- plot1.colors
polygon2.data["PlotColor"] <- plot2.colors
#total <- rbind(polygon.data, polygon1.data)
newdata <- polygon.data[order( polygon.data[,1], polygon.data[,2] ),]
newdata1 <- polygon1.data[order( polygon1.data[,1], polygon1.data[,2] ),]
newdata2 <- polygon2.data[order( polygon2.data[,1], polygon2.data[,2] ),]
#colnames(newdata) <- c("Chromosome", "chromStart", "chromEnd", "Data", "PlotColor")
#head(newdata)
#tail(newdata)

#newdata.color <- as.data.frame(newdata[,5])
#head(newdata.color)
#tail(newdata.color)
#RCircos.Polygon.Plot(newdata, track.num=1, data.col=4, genomic.columns=3, side="in", border.col = "red", is.sorted = TRUE)
#, inside.pos=1.6, outside.pos=1.8, is.sorted = TRUE)
RCircos.Polygon.Plot(newdata, track.num=1, data.col=4, genomic.columns=3, side="in", border.col = "red", is.sorted = TRUE)

RCircos.Polygon.Plot(newdata1, track.num=2, data.col=4, genomic.columns=3, side="in", border.col = "blue", is.sorted = TRUE)
RCircos.Polygon.Plot(newdata2, track.num=3, data.col=4, genomic.columns=3, side="in", border.col = "green", is.sorted = TRUE)




#RCircos.Polygon.Plot(polygon1.data, track.num=3, data.col=4, genomic.columns=3, side="in", border.col="red", polygon.col=NULL, is.sorted = TRUE)



    #   Close the graphic device and clear memory
    #   _________________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

dev.off();
message("RCircos Polygon Plot Demo Done!");
