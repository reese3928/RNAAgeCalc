

## This script generates all the internal data used in the package.
suppressMessages(library(AnnotationDbi))
library(org.Hs.eg.db)
suppressMessages(library(recount))

## get GTEx brain data from recount package
load(url("http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_brain.Rdata"))
lengthDB = rowData(rse_gene)
tempid = which(!duplicated(vapply(strsplit(rownames(lengthDB),split="\\."),
                                  '[', 1, FUN.VALUE=character(1))))
lengthDB = lengthDB[tempid,]
rownames(lengthDB) = vapply(strsplit(rownames(lengthDB),split="\\."), '[', 1,
                            FUN.VALUE=character(1))
lengthDB = lengthDB[,setdiff(colnames(lengthDB), c("gene_id","symbol")),
                    drop = FALSE]
lengthDB = data.frame(lengthDB)


library(glmnet)
library(devtools)
tissues = c("adipose_tissue", "adrenal_gland", "blood", "blood_vessel", "brain",
            "breast", "colon", "esophagus", "heart", "liver", "lung", "muscle",
            "nerve", "ovary", "pancreas", "pituitary", "prostate",
            "salivary_gland", "skin", "small_intestine", "spleen", "stomach",
            "testis", "thyroid", "uterus", "vagina")

signatures = c("DESeq", "Pearson", "Dev", "deMagalhaes", "GenAge",
               "GTExAge", "Peters", "all")

genelist_one = vector("list", length(signatures))
names(genelist_one) = signatures
genelist_all = vector("list", length(tissues))
names(genelist_all) = tissues
for(i in 1:length(genelist_all)) genelist_all[[i]] = genelist_one
genelist_cau = genelist_all

dir = "/Users/renxu/Desktop/Rpackage/RNAAgeData/all_samples"
for(t in tissues){
    message("tissue: ", t)
    for(sig in signatures){
        load(file.path(dir, t, paste0(sig,".RData")))
        coefs = coef(opt_model, s = "lambda.min")[,1]
        genelist_all[[t]][[sig]] = coefs[which(coefs!=0)]
    }
}

for(sig in signatures[-1]){
    load(file.path(dir, "all_tissues", paste0(sig,".RData")))
    coefs = coef(opt_model, s = "lambda.min")[,1]
    genelist_all[["all_tissues"]][[sig]] = coefs[which(coefs!=0)]
}

for(t in tissues) names(genelist_all[[t]])[1] = "DESeq2"

dir2 = "/Users/renxu/Desktop/Rpackage/RNAAgeData/Caucasian"
for(t in tissues){
    message("tissue: ", t)
    for(sig in signatures){
        load(file.path(dir2, t, paste0(sig,".RData")))
        coefs = coef(opt_model, s = "lambda.min")[,1]
        genelist_cau[[t]][[sig]] = coefs[which(coefs!=0)]
    }
}

for(sig in signatures[-1]){
    load(file.path(dir2, "all_tissues", paste0(sig,".RData")))
    coefs = coef(opt_model, s = "lambda.min")[,1]
    genelist_cau[["all_tissues"]][[sig]] = coefs[which(coefs!=0)]
}

for(t in tissues) names(genelist_cau[[t]])[1] = "DESeq2"

use_data(lengthDB, genelist_all, genelist_cau, internal = TRUE, 
         compress = "xz", overwrite = TRUE)


## some tests
load(file.path(dir, "breast", "DESeq.RData"))
load("data/FPKMExample.RData")
geneid = rownames(coef(opt_model, s = "lambda.min"))[-1]
fpkm1 = fpkm[geneid,]
dat.test = t(log2(fpkm1+1))
predict(opt_model, newx = dat.test, s = "lambda.min")

geneid2 = names(genelist_all[["breast"]][["DESeq"]])[-1]
fpkm2 = fpkm[geneid2,]
dat.test2 = t(log2(fpkm2+1))
coefs = genelist_all[["breast"]][["DESeq"]]
apply(dat.test2, 1, function(x) sum(x*coefs[-1])+coefs[1])




