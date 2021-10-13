library(methylKit)
library(rtracklayer)
library(GenomicRanges)
library(data.table)

pdf("diffmeth_methylKit.pdf")

###### input *_pe_CpG.txt files (comma separated, same directory as this script)
# for i in *CpG.txt; do printf '"'${i}'", '; done
file.list=list("")

###### input the names you wish to label the above files (in the same order)
id.list=list("")

###### Reads the file list, determine treatment vector (same order as id.list). 1 is treated, 0 is control.
###### Can put multiple different ids in the same vector. For example: (1,1,0,0,2,2)
###### Next 2 lines filter each CpG by coverage and normalizes
methylObj=methRead(file.list, sample.id=id.list, assembly="hg38", treatment=c(1,1,1,1,1,0,0,0,0,0,2,2), context="CpG", mincov=5)
#filter=filterByCoverage(methylObj, lo.count = 1)
myobj=normalizeCoverage(filter,method="median")

###### Comment out if sliding window analysis is not desired. Set for 200bp windows
#tiles=tileMethylCounts(myobj,win.size=200,step.size=200,cov.bases=5,mc.cores=12)

###### Plots methylation and coverage stats if desired.
if (FALSE)
{
getMethylationStats(myobj[[1]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj[[1]],plot=TRUE,both.strands=FALSE)
getMethylationStats(myobj[[2]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj[[2]],plot=TRUE,both.strands=FALSE)
getMethylationStats(myobj[[3]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj[[3]],plot=TRUE,both.strands=FALSE)
getMethylationStats(myobj[[4]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj[[4]],plot=TRUE,both.strands=FALSE)
#getMethylationStats(myobj[[5]],plot=TRUE,both.strands=FALSE)
#getCoverageStats(myobj[[5]],plot=TRUE,both.strands=FALSE)
#getMethylationStats(myobj[[6]],plot=TRUE,both.strands=FALSE)
#getCoverageStats(myobj[[6]],plot=TRUE,both.strands=FALSE)
}

###### Establish covariates here if desired
#covariates=data.frame(rep=c('Rep1','Rep23','Rep23','Rep1','Rep23','Rep23'))



#############################CpG Level Analysis####################################

###### input the absolute path of a tab-delimenated ROI/BED file (first 3 columns are #chr \t start \t stop)
#region=import.bed("---------------", extraCols=1)
#class(region)

###### Meat of the methylKit analysis. Finds DMCs/DMRs and outputs tables. Give file name outputs below
###### change qvalue and difference values if desired

# Find shared CpGs
meth=unite(methylObj, destrand=FALSE)
#meth=unite(tiles, destrand=FALSE) #### Uncomment line if DMRs are desired 

# Find DMCs
myDiff=calculateDiffMeth(meth, adjust="fdr", overdispersion="MN", test="Chisq", num.cores=12)
myDiff10p.hyper=getMethylDiff(myDiff,difference=10,qvalue=0.05,type="hyper")
myDiff10p.hypo=getMethylDiff(myDiff,difference=10,qvalue=0.05,type="hypo")
myDiff10p=getMethylDiff(myDiff,difference=10,qvalue=0.05)
per_myDiff=diffMethPerChr(myDiff,plot=TRUE,qvalue.cutoff=0.05, meth.cutoff=10)

# Write results to file
write.table(as(meth, "data.frame"), file="meth_shared.tsv", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)
write.table(as(myDiff, "data.frame"), file="diffMeth.tsv", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)
write.table(as(myDiff10p, "data.frame"), file="dmc.tsv", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)
write.table(as(myDiff10p.hyper, "data.frame"), file="dmc_hyper.tsv", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)
write.table(as(myDiff10p.hypo, "data.frame"), file="dmc_hypo.tsv", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)
write.table(as(per_myDiff, "data.frame"), file="diffMethPerChr.tsv", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)

###### Outputs correlation stats tables
#write.table(as(getCorrelation(meth,plot=FALSE), "data.frame"), file="methCorrelation.tsv", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
PCASamples(meth, screeplot=TRUE)
PCASamples(meth)
pr=PCASamples(meth, obj.return=TRUE)
# Report variances
print(pr$sdev^2/sum(pr$sdev^2))

###### Commented out
if (FALSE)
{
myDiff10p.hyper$strand="*"
myDiff10p.hypo$strand="*"


diff_hyper_gr=as(myDiff10p.hyper, "GRanges")
diff_hypo_gr=as(myDiff10p.hypo, "GRanges")


diff_promo_hyper=subsetByOverlaps(region, diff_hyper_gr)
diff_promo_hyper_cpgs=subsetByOverlaps(diff_hyper_gr, region)
diff_promo_hypo=subsetByOverlaps(region, diff_hypo_gr)
diff_promo_hypo_cpgs=subsetByOverlaps(diff_hypo_gr, region)

export.bed(diff_promo_hyper, "Amal_diff_1kb_hyper.3x3.bed")
export.bed(diff_promo_hypo, "Amal_diff_1kb_hypo.3x3.bed")
export.bed(diff_promo_hyper_cpgs, "Amal_diff_1kb_hyper_cpgs.3x3.bed")
export.bed(diff_promo_hypo_cpgs, "Amal_diff_1kb_hypo_cpgs.3x3.bed")
}
#######

if(FALSE)
{
###### Graphs distance to TSS, percentage of DMCs in various regions, etc.
# Look up the BED argument this takes (an example can be found here: /fs/scratch/ccri0063/amal/methylkit/mkitpermeth/refGene_mm10_readfeatures.bed)
gene.obj=readTranscriptFeatures("--------------")
diffAnn=annotateWithGeneParts(as(myDiff10p, "GRanges"), gene.obj)
head(getAssociationWithTSS(diffAnn))
getTargetAnnotationStats(diffAnn, percentage=FALSE, precedence=TRUE)
plotTargetAnnotation(diffAnn, precedence=TRUE, main="Differential Methylation Annotation (10p)")

cpg.obj=readFeatureFlank("--------------", feature.flank.name=c("Cgi","shores"))
diffCpGann=annotateWithFeatureFlank(as(myDiff10p,"GRanges"), cpg.obj$Cgi,cpg.obj$shores, feature.name="CpGi",flank.name="shores")
plotTargetAnnotation(diffCpGann,col=c("green","gray","white"), main="Differential Methylation Annotation (10p)")

tss=getAssociationWithTSS(diffAnn)
hist(tss$dist.to.feature[abs(tss$dist.to.feature)<=100000],main="Distance to Nearest TSS (10p)", xlab="Distance in bp", breaks=50, col="brown4")

#####bedgraph(myDiff10p, file.name="Amal_3x3_1kb_bedgraph.bed", col.name="meth.diff", unmeth=FALSE)
#####write.table(as(gene.obj, "data.frame"), file="gene_obj_table.txt", sep="\t", col.names=NA, row.names=TRUE)

}
##############################Regional Analysis####################################

##### Commented out because we have no clue how methylKit is calling these DMRs anymore lol
if(FALSE)
{
region=regionCounts(myobj, unique(region))

head(region[[1]])
meth_region=unite(region, destrand=FALSE)

myDiff_region=calculateDiffMeth(meth_region, adjust="fdr", overdispersion="MN", test="Chisq", num.cores=4, weighted.mean=FALSE)
class(myDiff_region)
myDiff_region10p.hyper=getMethylDiff(myDiff_region,difference=10,qvalue=0.05,type="hyper")
myDiff_region10p.hypo=getMethylDiff(myDiff_region,difference=10,qvalue=0.05,type="hypo")
myDiff_region10p=getMethylDiff(myDiff_region,difference=10,qvalue=0.05)
per_myDiff_region=diffMethPerChr(myDiff_region,plot=FALSE,qvalue.cutoff=0.05, meth.cutoff=10)

class(as(myDiff_region, "data.frame"))

write.table(as(myDiff_region, "data.frame"), file="autophagy_1kb_Amal_3x3_diffMeth.txt", sep="\t", col.names=NA, row.names=TRUE)
write.table(as(myDiff_region10p, "data.frame"), file="autophagy_1kb_Amal_3x3_diffMeth_10p.txt", sep="\t", col.names=NA, row.names=TRUE)
write.table(as(myDiff_region10p.hyper, "data.frame"), file="autophagy_1kb_Amal_3x3_diffMeth_10pHyper.txt", sep="\t", col.names=NA, row.names=TRUE)
write.table(as(myDiff_region10p.hypo, "data.frame"), file="autophagy_1kb_Amal_3x3_diffMeth_10pHypo.txt", sep="\t", col.names=NA, row.names=TRUE)
}
