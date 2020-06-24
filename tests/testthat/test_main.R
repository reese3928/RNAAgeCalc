context("Test pred_age, makeplot")

data(fpkmExample)
library(SummarizedExperiment)

test_that("Test invalid input", {
    exprdata1 = 1
    expect_error(predict_age(exprdata = exprdata1, exprtype = "FPKM"))
    
    exprdata2 = matrix(c("a", "b", "c", "d"), 2, 2)
    expect_error(predict_age(exprdata = exprdata2, exprtype = "FPKM"))
    
    exprdata3 = matrix(c(1:8, -1), 3, 3)
    rownames(exprdata3) = c("a", "b", "c")
    colnames(exprdata3) = c("sample1", "sample2", "sample3")
    expect_error(predict_age(exprdata = exprdata3, exprtype = "FPKM"))
    
    exprdata4 = data.frame(sample1 = c(1,2,3), sample2 = c("a", "b", "c"),
                           sample3 = as.factor(c("ctr", "ctr", "trt")))
    expect_error(predict_age(exprdata = exprdata4, exprtype = "FPKM"))
    
    exprdata5 = data.frame(matrix(c(1:8, -1), 3, 3))
    expect_error(predict_age(exprdata = exprdata5, exprtype = "FPKM"))
    
    exprdata6 = matrix(1:9, 3, 3)
    expect_error(predict_age(exprdata = exprdata6))
    rownames(exprdata6) = c("a", "b", "c")
    expect_error(predict_age(exprdata = exprdata6))
    colnames(exprdata6) = c("sample1", "sample2", "sample3")
    expect_error(predict_age(exprdata = exprdata6, exprtype = "unknown"))
    
    expect_error(predict_age(exprdata = fpkm, exprtype = "FPKM", 
                             stype = "unknown"))
    
    expect_error(predict_age(exprdata = fpkm, exprtype = "FPKM", 
                             maxp = "unknown"))
    
    se1 = 1
    expect_error(predict_age_fromse(se = se1, exprtype = "FPKM"))
    
    colData = data.frame(age = c(40, 50))
    se2 = SummarizedExperiment(assays=list(unknown=fpkm), colData=colData)
    expect_error(predict_age_fromse(se = se2, exprtype = "FPKM"))
    
    se3 = SummarizedExperiment(assays=list(FPKM=fpkm), colData=colData)
    expect_error(predict_age_fromse(se = se3, exprtype = "counts"))
    
    se4 = SummarizedExperiment(assays=list(counts=fpkm), colData=colData)
    expect_error(predict_age_fromse(se = se4, exprtype = "FPKM"))
    
    # These are some intermediate tests when we are writing the 
    # predict_age_fromse() function. It no longer holds once we finish writing
    # the function. Thus, they are comment.
    # rowData = DataFrame(bp_length = runif(nrow(fpkm)))
    # se5 = SummarizedExperiment(assays=list(counts=fpkm), colData=colData, 
    #                            rowData = rowData)
    # res = predict_age_fromse(se5, tissue = NULL, exprtype="counts")
    # expect_identical(res, rowData[["bp_length"]])
    
    # se5 = SummarizedExperiment(assays=list(fpkm=fpkm), colData=colData)
    # res = predict_age_fromse(se5, tissue = NULL)
    # expect_identical(res, data.frame(sampleid = colnames(fpkm), 
    #                                 age = c(40, 50)))
    
    se5 = SummarizedExperiment(assays=list(fpkm=fpkm), colData=colData)
    expect_error(predict_age_fromse(se = se5, exprtype = "unknown"))
    
    expect_error(predict_age_fromse(se = se5, exprtype = "FPKM", 
                                    stype = "unknown"))
    
    expect_error(predict_age_fromse(se = se5, exprtype = "FPKM", 
                                    maxp = "unknown"))
    
    res1 = predict_age(exprdata = fpkm, tissue = "brain", signature = "GTExAge", 
                       chronage = data.frame(sampleid = colnames(fpkm), 
                                             age = c(40, 50)))
    res2 = predict_age_fromse(se = se5, tissue = "brain", 
                              signature = "GTExAge")
    expect_identical(res1, res2)
    
    se6 = SummarizedExperiment(assays=list(fpkm=fpkm))
    res1 = predict_age(exprdata = fpkm, tissue = "brain", 
                       signature = "GTExAge")
    res2 = predict_age_fromse(se = se6, tissue = "brain", 
                              signature = "GTExAge")
    expect_identical(res1, res2)
})


test_that("Test tissue, signature, idtype, and chronage arguments", {
    expect_error(predict_age(exprdata = fpkm, exprtype = "FPKM"), NA)
    expect_error(predict_age(exprdata = fpkm, tissue = "breaste",
        exprtype = "FPKM"), NA)
    expect_error(predict_age(exprdata = fpkm, exprtype = "FPKM",
        idtype = "symbole"))
    
    colData = data.frame(age = c(40, 50))
    se1 = SummarizedExperiment(assays=list(fpkm=fpkm), colData=colData)
    expect_error(predict_age_fromse(se = se1, exprtype = "FPKM"), NA)
    expect_error(predict_age_fromse(se = se1, tissue = "breaste",
                                    exprtype = "FPKM"), NA)
    expect_error(predict_age_fromse(se = se1, exprtype = "FPKM",
                                    idtype = "symbole"))
    
    expect_error(predict_age(exprdata = fpkm, idtype = "SYMBOL",
        signature = "unknown"))
    expect_error(predict_age(exprdata = fpkm, idtype = "SYMBOL",
        signature = "DESeq2"),NA)
    expect_error(predict_age(exprdata = fpkm, idtype = "SYMBOL",
        signature = "DESeq"),NA)
    expect_error(predict_age(exprdata = fpkm, idtype = "SYMBOL",
        tissue = "brain"),NA)
    expect_error(predict_age(exprdata = fpkm, idtype = "SYMBOL",
        tissue = "brain", signature = "unknown"))
    
    temp1 = predict_age(exprdata = fpkm, idtype = "SYMBOL",
        signature = "DESeq2")
    temp2 = predict_age(exprdata = fpkm, idtype = "SYMBOL",
        signature = "DESeq")
    expect_identical(temp1, temp2)
    
    expect_error(predict_age_fromse(se = se1, idtype = "SYMBOL",
                                    signature = "unknown"))
    expect_error(predict_age_fromse(se = se1, idtype = "SYMBOL",
                                    signature = "DESeq2"),NA)
    expect_error(predict_age_fromse(se = se1, idtype = "SYMBOL",
                                    signature = "DESeq"),NA)
    expect_error(predict_age_fromse(se = se1, idtype = "SYMBOL",
                                    tissue = "brain"),NA)
    expect_error(predict_age_fromse(se = se1, idtype = "SYMBOL", 
                                    tissue = "brain", signature = "unknown"))
    
    temp1 = predict_age_fromse(se = se1, idtype = "SYMBOL", 
        signature = "DESeq2")
    temp2 = predict_age_fromse(se = se1, idtype = "SYMBOL", 
        signature = "DESeq")
    expect_identical(temp1, temp2)
    

    expect_error(predict_age(exprdata = fpkm, idtype = "SYMBOL",
        tissue = "unknown", signature = "unknown"))
    expect_error(predict_age(exprdata = fpkm, idtype = "SYMBOL",
        tissue = "unknown", signature = "DESeq2"), NA)
    expect_error(predict_age(exprdata = fpkm, idtype = "SYMBOL",
        tissue = "unknown", signature = "DESeq"), NA)
    
    expect_error(predict_age_fromse(se = se1, idtype = "SYMBOL",
        tissue = "unknown", signature = "DESeq2"), NA)
    expect_error(predict_age_fromse(se = se1, idtype = "SYMBOL",
        tissue = "unknown", signature = "DESeq"), NA)

    expect_error(predict_age(exprdata = fpkm, idtype = "SYMBOL",
        tissue = "breast"),NA)
    expect_error(predict_age(exprdata = fpkm, idtype = "SYMBOL",
        tissue = "unknown"),NA)
    
    expect_error(predict_age_fromse(se = se1, idtype = "SYMBOL",
                             tissue = "breast"),NA)
    expect_error(predict_age_fromse(se = se1, idtype = "SYMBOL",
                             tissue = "unknown"),NA)

    ## test the chronage argument
    chronage = data.frame(sampleid = colnames(fpkm), age = c(30,50))
    expect_error(predict_age(exprdata = fpkm, exprtype = "FPKM",
        chronage = chronage),NA)

    res = predict_age(exprdata = fpkm, exprtype = "FPKM")
    expect_error(makeplot(res))

    res = predict_age(exprdata = fpkm, exprtype = "FPKM", chronage = chronage)
    expect_warning(makeplot(res))

    #chronage1 = data.frame(sampleid = c(30,50), age = c(30,50))
    #expect_error(predict_age(exprdata = fpkm, exprtype = "FPKM",
    #    chronage = chronage1))

    chronage2 = data.frame(sampleid = colnames(fpkm),
        age = c("unknown", "unknown"))
    expect_error(predict_age(exprdata = fpkm, exprtype = "FPKM",
        chronage = chronage2))

    chronage3 = data.frame(sampleid = c("unknown1","unknown2"), age = c(30,50))
    expect_error(predict_age(exprdata = fpkm, exprtype = "FPKM",
        chronage = chronage3))

    chronage4 = data.frame(sampleid = colnames(fpkm), age = c(NA,-3))
    expect_error(predict_age(exprdata = fpkm,
        exprtype = "FPKM", chronage = chronage4))

    chronage5 = data.frame(sampleid = colnames(fpkm),
                           age = c(30,50), age2 = c(30,50))
    expect_error(predict_age(exprdata = fpkm, exprtype = "FPKM",
        chronage = chronage5),NA)

    fpkm_large = cbind(fpkm, fpkm, fpkm, fpkm)
    fpkm_large = cbind(fpkm_large, fpkm_large, fpkm_large, fpkm_large)
    colnames(fpkm_large) = paste0("sample",1:32)
    chronage6 = data.frame(sampleid = colnames(fpkm_large), age = 1:32)
    expect_error(predict_age(exprdata = fpkm_large, exprtype = "FPKM",
        chronage = chronage6),NA)
    res = predict_age(exprdata = fpkm_large, exprtype = "FPKM",
                      chronage = chronage6)
    expect_error(makeplot(res), NA)
    
    chronage7 = rbind(chronage6, data.frame(sampleid = "sample2", age = 2))
    expect_error(predict_age(exprdata = fpkm_large, exprtype = "FPKM",
                             chronage = chronage7))

    ## check for imputation
    tissue = "brain"
    signature = "DESeq2"
    genes_required = names(genelist_all[[tissue]][[signature]])[-1]
    fpkm2 = fpkm
    fpkm2 = fpkm2[which(!rownames(fpkm2)%in%genes_required[c(1,10,20,30)]),]
    fpkm2 = fpkm2[which(rownames(fpkm2)%in%genes_required),]
    ## fpkm2 contains missing genes
    expect_error(predict_age(fpkm2, "brain"),NA)

    ## partial row NA
    fpkm3 = fpkm
    fpkm3 = fpkm3[which(rownames(fpkm3)%in%genes_required),]
    fpkm3[1,2] = NA
    fpkm3[3,1] = NA
    expect_error(predict_age(fpkm3, "brain"),NA)

    ## some other tests:
    ## These are some tests for our internal data, it has nothing to do
    ## with the user facing functions. Since some of the code takes too long,
    ## these tests are commented in order to pass the Bioconductor check
    ## time limit

    ## test whether the id type conversion work well
    #suppressMessages(library(recount))
    ## get GTEx brain data from recount package
    #load(url("http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_brain.Rdata"))
    #rse = scale_counts(rse_gene)
    ## get RPKM matrix
    #fpkm <- getRPKM(rse)

    ## this is the gene map used in our paper
    #gene_map = read.csv("https://raw.githubusercontent.com/reese3928/RNAAgeData/master/gene_map_03182019.csv")

    #tempid = which(!duplicated(vapply(strsplit(rownames(fpkm),split="\\."), '[', 1, FUN.VALUE=character(1))))
    #fpkm = fpkm[tempid,]
    #rownames(fpkm) = vapply(strsplit(rownames(fpkm),split="\\."), '[', 1, FUN.VALUE=character(1))
    #fpkm2 = fpkm[sample(nrow(fpkm)),]

    #exprdata = fpkm2
    #idtype = "ENSEMBL"
    #mapid = suppressMessages(select(org.Hs.eg.db, rownames(exprdata),
    #                                columns = "SYMBOL", keytype = idtype))
    #mapid = mapid[which(!is.na(mapid$SYMBOL)),]
    #mapid = mapid[!duplicated(mapid[,idtype]),]
    #ind1 = which(!mapid$SYMBOL %in% mapid$SYMBOL[duplicated(mapid$SYMBOL)])
    #mapid = mapid[ind1,]
    #rownames(mapid) = as.character(mapid[,idtype])
    #exprdata = exprdata[rownames(exprdata)%in%rownames(mapid),]
    #genesymbol = mapid[rownames(exprdata),"SYMBOL"]
    #rownames(exprdata) = genesymbol
    #all(gene_map$SYMBOL%in%rownames(exprdata))

})








