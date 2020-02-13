

## This script generates example raw count, FPKM data
suppressMessages(library(recount))
## get GTEx brain data from recount package
load(url("http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_brain.Rdata"))
rse = scale_counts(rse_gene)
## get RPKM matrix
fpkm <- getRPKM(rse)

## this is the gene map used in our paper
gene_map = read.csv("https://raw.githubusercontent.com/reese3928/RNAAgeData/master/gene_map_03182019.csv")

tempid = which(!duplicated(vapply(strsplit(rownames(fpkm),split="\\."), '[', 1, FUN.VALUE=character(1))))
fpkm = fpkm[tempid,]
rownames(fpkm) = vapply(strsplit(rownames(fpkm),split="\\."), '[', 1, FUN.VALUE=character(1))
fpkm = fpkm[as.character(gene_map$ENSEMBL),]
rownames(fpkm) = as.character(gene_map$SYMBOL)
fpkm = fpkm[,1:2]  ## subset the samples to make the example size smaller
save(fpkm, file = "data/fpkmExample.RData", compress = "xz")


rawcount = assay(rse)
tempid = which(!duplicated(vapply(strsplit(rownames(rawcount),split="\\."), '[', 1, FUN.VALUE=character(1))))
rawcount = rawcount[tempid,]
rownames(rawcount) = vapply(strsplit(rownames(rawcount),split="\\."), '[', 1, FUN.VALUE=character(1))
rawcount = rawcount[as.character(gene_map$ENSEMBL),]
rownames(rawcount) = as.character(gene_map$SYMBOL)
rawcount = rawcount[,1:2]  ## subset the samples to make the example size smaller
save(rawcount, file = "data/countExample.RData", compress = "xz")




