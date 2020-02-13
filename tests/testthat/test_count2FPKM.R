context("Test count2FPKM")

data(countExample)

test_that("Test count2FPKM", {
    rawcount2 = rawcount
    rownames(rawcount2)[c(1000,2000,3000)] = c("ABC1", "ABC2", "ABC3")
    expect_warning(count2FPKM(rawcount2, genelength = NULL))

    genelength1 = rep("unknown", nrow(rawcount))
    expect_error(count2FPKM(rawcount, genelength = genelength1))

    genelength2 = matrix(1:9,3,3)
    expect_error(count2FPKM(rawcount, genelength = genelength2))

    genelength3 = rep(1, nrow(rawcount)-1)
    expect_error(count2FPKM(rawcount, genelength = genelength3))

    rawcount3 = rawcount[1:5,]
    rownames(rawcount3) = c("ENSG00000000003", "ENSG00000000004",
        "ENSG00000000005", "ENSG00000000419", "ENSG00000000457")
    expect_warning(count2FPKM(rawcount3, idtype = "ENSEMBL"))

})
