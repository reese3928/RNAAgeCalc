#' @title An example of FPKM data
#'
#' @description This is an example of FPKM data. It is the gene expression 
#' FPKM data from The genotype-tissue expression (GTEx) project brain samples. 
#' The data was downloaded from 
#' recount2(https://jhubiostatistics.shinyapps.io/recount/). For 
#' illustration purpose, only two samples were included. The R script to 
#' generate this example data can be found in inst/scripts/createExamples.R.
#' @docType data
#' @name fpkm
#' @usage fpkm
#' @format A data frame which contains the fpkm data. The dimension of this 
#' data frame is 24989 by 2. Each row is a gene and each column is a sample. 
#' @keywords datasets
#' @references Collado-Torres, Leonardo, et al. "Reproducible RNA-seq analysis 
#' using recount2." Nature biotechnology 35.4 (2017): 319-321.
#' @references Lonsdale, John, et al. "The genotype-tissue expression (GTEx) 
#' project." Nature genetics 45.6 (2013): 580.
NULL

#' @title An example of RNASeq counts data
#'
#' @description An example of RNASeq counts data. It is the gene expression 
#' counts data from The genotype-tissue expression (GTEx) project brain samples. 
#' The data was downloaded from 
#' recount2(https://jhubiostatistics.shinyapps.io/recount/). For 
#' illustration purpose, only two samples were included. The R script to 
#' generate this example data can be found in inst/scripts/createExamples.R.
#' @docType data
#' @name rawcount
#' @usage rawcount
#' @format A data frame which contains the RNASeq counts data. The dimension of 
#' this data frame is 24989 by 2. Each row is a gene and each column is a 
#' sample. 
#' @keywords datasets
#' @references Collado-Torres, Leonardo, et al. "Reproducible RNA-seq analysis 
#' using recount2." Nature biotechnology 35.4 (2017): 319-321.
#' @references Lonsdale, John, et al. "The genotype-tissue expression (GTEx) 
#' project." Nature genetics 45.6 (2013): 580.
NULL