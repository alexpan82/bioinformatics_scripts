library(methylKit)

args <- commandArgs(TRUE)

filename <- args[1]
samplename <- args[2]
genome <- args[3]

print(filename)
print(samplename)
print(genome)

if (is.null(filename) | is.null(samplename)){
	print("Usage: [filename] [samplename]")
}else{
	myobj <- processBismarkAln(filename, sample.id = samplename, assembly = genome, nolap = TRUE, save.context= c("CpG"), mincov = 1, save.folder=getwd())
}
