
#' @title Calculate RNA age using SummarizedExperiment
#' @description This function takes SummarizedExperiment object as input and 
#' calculates RNA age. 
#' @param se a \code{\link{SummarizedExperiment}} object. The assays(se) should
#' contain gene expression data. The name of assays(se) should be either "FPKM"
#' or "counts". Use `exprtype` argument to specify the type of gene expression 
#' data provided. Users are able to provide the chronological age of samples
#' using colData(se). This is optional. If provided, the column name for 
#' chronological age in colData(se) should be "age". If some samples' 
#' chronological age are not available, users are expected to set the 
#' chronological age in colData(se) to NA. If chronological age is not 
#' provided, the age acceleration residual will not be calculated. See package 
#' vignette for the definition of age acceleration residual. In addition, users 
#' are able to provide their own gene length using rowData(se). This is also
#' optional. If using `exprtype = "FPKM"`, the provided gene length will be 
#' ignored. If provided, the column name for gene length in rowData(se) should 
#' be "bp_length". The function will convert raw count to FPKM by the
#' user-supplied gene length. Otherwise, gene length is obtained from the
#' internal database. See below for an example of se object.
#' @param tissue a string indicate which tissue the gene expression data is
#' obtained from. Users are expected to provide one of the following tissues.
#' If the tissue argument is not provide or the provided tissue is not in this
#' list, then the age predictor trained on all tissues will be used to
#' calculate RNA age.
#' \itemize{
#'   \item adipose_tissue
#'   \item adrenal_gland
#'   \item blood
#'   \item blood_vessel
#'   \item brain
#'   \item breast
#'   \item colon
#'   \item esophagus
#'   \item heart
#'   \item liver
#'   \item lung
#'   \item muscle
#'   \item nerve
#'   \item ovary
#'   \item pancreas
#'   \item pituitary
#'   \item prostate
#'   \item salivary_gland
#'   \item skin
#'   \item small_intestine
#'   \item spleen
#'   \item stomach
#'   \item testis
#'   \item thyroid
#'   \item uterus
#'   \item vagina
#' }
#' @param exprtype either "counts" or "FPKM". For RPKM data, please use
#' `exprtype` = "FPKM".
#' @param idtype a string which indicates the gene id type in `exprdata`.
#' It should be one of "SYMBOL", "ENSEMBL", "ENTREZID" or "REFSEQ". Default is
#' "SYMBOL".
#' @param stype a string which specifies which version of pre-trained 
#' calculators to be used. It should be either "all" or "Caucasian". "all" 
#' means samples from all races (American Indian/Alaska Native, Asian, 
#' Black/African American, and Caucasian) are used to obtain the 
#' pre-trained calculator. "Caucasian" means only the Caucasian samples are 
#' used to build up the pre-trained calculator.
#' @param signature a string which indicates the age signature to use when
#' calculating RNA age. This argument is not required. \cr
#' In the case that this argument is not provided, if `tissue` argument is also
#' provided and the tissue is in the list above, the tissue specific age
#' signature given by our DESeq2 analysis result on GTEx data will be used.
#' Otherwise, the across tissue signature "GTExAge" will be used. \cr
#' In the case that this argument is provided, it should be one of the
#' following signatures. A detailed description of the meaning of these
#' signatures is given in the package vignette.
#' \itemize{
#'   \item DESeq2
#'   \item Pearson
#'   \item Dev
#'   \item deMagalhaes
#'   \item GenAge
#'   \item GTExAge
#'   \item Peters
#'   \item all
#' }
#' @param maxp the maxp argument used in \code{\link[impute]{impute.knn}} 
#' function. This is optional. 
#' @return a data frame contains RNA age.
#' @importFrom SummarizedExperiment assayNames assay colData rowData
#' @importFrom methods is
#' @export
#' @examples 
#' library(SummarizedExperiment)
#' data(fpkmExample)
#' colData = data.frame(age = c(40, 50))
#' se = SummarizedExperiment(assays=list(FPKM=fpkm), 
#' colData=colData)
#' res = predict_age_fromse(se = se, exprtype = "FPKM")

predict_age_fromse <- function(se, tissue, exprtype = c("FPKM", "counts"), 
    idtype = c("SYMBOL", "ENSEMBL", "ENTREZID", "REFSEQ"),
    stype = c("all", "caucasian"), signature = NULL, maxp = NULL){
    
    if (!is(se, "SummarizedExperiment")) {
        stop("se should be a SummarizedExperiment object. ")
    }
    
    exprtype = match.arg(exprtype)
    ind = which(toupper(assayNames(se))==toupper(exprtype))
    if(length(ind)==0){
        stop("The assay name should contain ", exprtype)
    }
    
    exprdata = assay(se, ind)
    
    ## process gene length
    genelength = NULL
    if(exprtype=="counts"){
        lengthind = which(names(rowData(se))=="bp_length")
        if(lengthind > 0){
            genelength = rowData(se)[[lengthind]]
        }
    }
    
    ## process chronological age
    chronage = NULL
    ageind = which(toupper(names(colData(se)))=="AGE")
    if(length(ageind)>0){
        chronage = data.frame(sampleid = rownames(colData(se)), 
                              age = colData(se)[[ageind]])
    }
    
    if(missing(tissue)){
        res = predict_age(exprdata, exprtype = exprtype, idtype = idtype,
                          stype = stype, signature = signature, 
                          genelength = genelength, chronage = chronage, 
                          maxp = maxp)
    }else{
        res = predict_age(exprdata, tissue = tissue, 
                          exprtype = exprtype, idtype = idtype,
                          stype = stype, signature = signature, 
                          genelength = genelength, chronage = chronage,
                          maxp = maxp)
    }
    
    res
}





