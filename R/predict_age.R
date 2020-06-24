
#' @title Calculate RNA age
#' @description This function calculates RNA age based on pre-trained
#' calculators.
#' @param exprdata a matrix or data frame which contains gene expression data
#' with each row represents a gene and each column represents a sample. Use
#' the argument `exprtype` to specify raw count or FPKM. The rownames of
#' `exprdata` should be gene ids and colnames of `exprdata` should be sample
#' ids.
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
#' @param genelength a vector which contains gene length in bp. 
#' The size of `genelength` should be
#' equal to the number of rows in `exprdata`. This argument is optional. If
#' using `exprtype = "FPKM"`, `genelength` argument is ignored. If using
#' `exprtype = "counts"`, the raw count will be converted to FPKM. If
#' `genelength` is provided, the function will convert raw count to FPKM by the
#' user-supplied gene length. Otherwise, gene length is obtained from the
#' internal database.
#' @param chronage a data frame which contains the chronological age of each
#' sample. This argument is optional. If provided, it should be a dataframe
#' with 1st column sample id and 2nd column chronological age. The sample order
#' in `chronage` doesn't have to be in the same order as in
#' `exprdata`. However, the samples in `chronage` and `exprdata` should be the
#' same. If some samples' chronological age are not available, users are
#' expected to set the chronological age in `chronage` to NA. If `chronage`
#' contains more than 2 columns, only the first 2 columns will be considered.
#' If this argument is not provided, the age acceleration residual will not be
#' calculated. See package vignette for the definition of age acceleration
#' residual.
#' @param maxp the maxp argument used in \code{\link[impute]{impute.knn}} 
#' function. This is optional. 
#' @importFrom stats lm na.omit residuals
#' @importFrom impute impute.knn
#' @return a data frame contains RNA age.
#' @export
#' @examples
#' data(fpkmExample)
#' res = predict_age(exprdata = fpkm, exprtype = "FPKM")


predict_age <- function(exprdata, tissue, exprtype = c("FPKM","counts"), 
    idtype = c("SYMBOL", "ENSEMBL", "ENTREZID", "REFSEQ"),
    stype = c("all", "caucasian"), signature = NULL, genelength = NULL, 
    chronage = NULL, maxp = NULL){

    ## check input:
    stopifnot(is.matrix(exprdata)|is.data.frame(exprdata))
    if(is.matrix(exprdata)){
        if(!is.numeric(exprdata)){
            stop("Only numeric values are allowed in the exprdata matrix.")
        }
        if(is.null(rownames(exprdata))){
            stop("Row names should be provided in the exprdata matrix.")
        }
        if(is.null(colnames(exprdata))){
            stop("Column names should be provided in the exprdata matrix.")
        }
        if( any(duplicated(rownames(exprdata))) ){
            stop("Duplicated gene names found in the exprdata matrix.")
        }
        exprdata = data.frame(exprdata)
    }else{
        if(any(!vapply(exprdata, is.numeric, TRUE))){
            stop("Only numeric values are allowed in the exprdata data.frame.")
        }
    }
    if(any(exprdata<0, na.rm=TRUE)){
        stop("Gene expression data cannot contain negative value(s).")
    }
    if(any(duplicated(rownames(exprdata)), na.rm = TRUE)){
        stop("Duplicated gene names found in the exprdata matrix.")
    }

    exprtype = match.arg(exprtype)
    idtype = match.arg(idtype)
    stype = match.arg(stype)

    if(!is.null(chronage)){
        stopifnot(is.data.frame(chronage))
        if(ncol(chronage)>2){
            message("More than 2 columns are provided in chronage. Only the
first 2 columns will be used.")
        }
        if(!is.character(chronage[,1]) & !is.factor(chronage[,1]) & 
            !is.numeric(chronage[,1]) ){
            stop("The 1st column in chronage should be sample ids.")
        }
        if(any(duplicated(chronage[,1]))){
            stop("chronage contains duplicated sample ids.")
        }
        if(!is.numeric(chronage[,2])){
            stop("The 2nd column in chronage should be chronological age.")
        }
        if(!setequal(chronage[,1], colnames(exprdata))){
            stop("Samples in chronage and exprdata should be the same.")
        }
        if(any(chronage[,2]<0, na.rm=TRUE)){
            stop("Chronological age contains negative value(s).")
        }
    }

    if(!is.null(maxp)){
        stopifnot(is.numeric(maxp))
        if(length(maxp) > 1){
            message("maxp contains more than 1 element. Only the first one
                    will be used.")
            maxp = maxp[1]
        }
    }
    ## end check input

    ## process tissue and signature argument
    if(missing(tissue)){
        message("tissue is not provided, using the RNA age predictor trained
by all tissues automatically.")
        tissue = "all_tissues"
        if(is.null(signature)){
            message("signature is not provided, using GTExAge signature
automatically.")
            signature = "GTExAge"
        }else{
            signature = match.arg(signature, c("DESeq", "DESeq2", "Pearson", 
                "Dev", "deMagalhaes", "GenAge", "GTExAge", "Peters", "all"))
            if(signature=="DESeq"){
                signature = "DESeq2"  ## replace "DESeq" with "DESeq2"
            }
            if(signature=="DESeq2"){
                message("DESeq2 signature is currently not available for
all tissues, using Pearson signature automatically.")
                signature = "Pearson"
            }

        }
    }else{
        stopifnot(is.character(tissue))
        stopifnot(length(tissue)==1)
        if( tolower(tissue) %in% names(genelist_all)){
            tissue = tolower(tissue)
            if(is.null(signature)){
                message("signature is not provided, using DESeq2 signature
automatically.")
                signature = "DESeq2"
            }else{
                signature = match.arg(signature, c("DESeq", "DESeq2", 
                    "Pearson", "Dev", "deMagalhaes", "GenAge", "GTExAge", 
                    "Peters", "all"))
                if(signature=="DESeq"){
                    signature = "DESeq2"  ## replace "DESeq" with "DESeq2"
                }
            }
        }else{
            message("the provided tissue was not found, using the RNA age
predictor trained by all tissues automatically.")
            tissue = "all_tissues"
            if(is.null(signature)){
                message("signature is not provided, using GTExAge signature
automatically.")
                signature = "GTExAge"
            }else{
                signature = match.arg(signature, c("DESeq", "DESeq2", 
                    "Pearson", "Dev", "deMagalhaes", "GenAge", "GTExAge", 
                    "Peters", "all"))
                if(signature=="DESeq"){
                    signature = "DESeq2"  ## replace "DESeq" with "DESeq2"
                }
                if(signature=="DESeq2"){
                    message("DESeq2 signature is currently not available for
all tissues, using Pearson signature automatically.")
                    signature = "Pearson"
                }

            }

        }
    }
    ## end: process tissue and signature argument

    if(exprtype=="counts"){
        exprdata = count2FPKM(exprdata, genelength = genelength,
            idtype = idtype)
    }

    ##  convert gene id type to gene symbol
    if(idtype!="SYMBOL"){
        mapid = suppressMessages(select(org.Hs.eg.db, rownames(exprdata),
            columns = "SYMBOL", keytype = idtype))
        mapid = mapid[which(!is.na(mapid$SYMBOL)),]
        mapid = mapid[!duplicated(mapid[,idtype]),]
        ind1 = which(!mapid$SYMBOL %in% mapid$SYMBOL[duplicated(mapid$SYMBOL)])
        mapid = mapid[ind1,]
        rownames(mapid) = as.character(mapid[,idtype])
        exprdata = exprdata[rownames(exprdata)%in%rownames(mapid),]
        genesymbol = mapid[rownames(exprdata),"SYMBOL"]
        rownames(exprdata) = genesymbol
    }

    if(stype=="all"){
        genes_required = names(genelist_all[[tissue]][[signature]])[-1]
    }else{
        genes_required = names(genelist_cau[[tissue]][[signature]])[-1]
    }
        
    sig_in_expr = genes_required%in%rownames(exprdata)
    if(sum(!sig_in_expr)!=0){
        message(round(sum(!sig_in_expr)/length(genes_required)*100,4),
                "% genes in the gene signature are not included in the supplied
                gene expression.")
        
        ## impute the gene expression in the log scale
        message("performing imputation...")
        tempmat = data.frame(matrix(NA, sum(!sig_in_expr), ncol(exprdata)))
        rownames(tempmat) = setdiff(genes_required, rownames(exprdata))
        colnames(tempmat) = colnames(exprdata)
        
        exprdata_withNA = rbind(exprdata, tempmat)
        exprdata_log = log2(exprdata_withNA + 1)
        exprdata_log = as.matrix(exprdata_log)
        
        if(is.null(maxp)){
            maxp = nrow(exprdata_log)
        }
        tmp = suppressWarnings(impute.knn(exprdata_log, maxp=maxp))
        exprdata_log_impute = data.frame(tmp$data)
        exprdata_sub = t(exprdata_log_impute[genes_required,])
        
    }else{
        ##  check partial row NA, then applied pre-trained calculator
        if(any(is.na(exprdata[genes_required,]))){
            message("performing imputation...")
            exprdata_log = log2(exprdata + 1)
            exprdata_log = as.matrix(exprdata_log)
            
            if(is.null(maxp)){
                maxp = nrow(exprdata_log)
            }
            tmp = suppressWarnings(impute.knn(exprdata_log, maxp=maxp))
            exprdata_log_impute = data.frame(tmp$data)
            exprdata_sub = t(exprdata_log_impute[genes_required,])
        }else{
            exprdata_sub = t(log2(exprdata[genes_required,] + 1))
        }
    }
    
    ## make prediction using pre-trained calculator
    if(stype=="all"){
        coefs = genelist_all[[tissue]][[signature]]
    }else{
        coefs = genelist_cau[[tissue]][[signature]]
    }
  
    RNAAge = apply(exprdata_sub, 1, function(x) sum(x*coefs[-1])+coefs[1])
    res = data.frame(RNAAge = RNAAge)
    if(!is.null(chronage)){
        matchid = match(colnames(exprdata), chronage[,1])
        res$ChronAge = chronage[matchid,2]
        ## if sample size is too small, age acceleration residual cannot be
        ## calculated
        if(nrow(na.omit(res))>30){
            res$AgeAccelResid = NA
            res[!is.na(res$ChronAge),"AgeAccelResid"] =
                residuals(lm(res$RNAAge~res$ChronAge))
        }
    }

    res

}


