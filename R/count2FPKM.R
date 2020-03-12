
#' @title Converting gene expression data from raw count to FPKM
#' @description This function converts gene expression data from raw count to
#' FPKM by using \code{\link[recount]{getRPKM}} 
#'
#' @param rawcount a matrix or data frame which contains gene expression counts
#' data.
#' @param genelength gene length in bp. The size of `genelength` should be
#' equal to the number of rows in `rawcount`. This argument is optional. If
#' not provided, gene length is obtained from the internal database.
#' @param idtype a string which indicates the gene id type in rawcount matrix.
#' It should be one of "SYMBOL", "ENSEMBL", "ENTREZID" or "REFSEQ".
#' Default is "SYMBOL".
#' @export
#' @import org.Hs.eg.db
#' @importFrom AnnotationDbi select
#' @importFrom recount getRPKM
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @return a data frame contains FPKM.
#' @references Carlson M (2019). org.Hs.eg.db: Genome wide annotation for
#' Human. R package version 3.8.2.
#' @references Collado-Torres, Leonardo, et al. "Reproducible RNA-seq analysis 
#' using recount2." Nature biotechnology 35.4 (2017): 319-321.
#' @examples
#' data(countExample)
#' head(rawcount)
#' fpkm = count2FPKM(rawcount)
#' head(fpkm)

count2FPKM = function(rawcount, genelength = NULL, idtype = "SYMBOL"){
    if(is.null(genelength)){
        if(idtype!="ENSEMBL"){
            ## convert the gene id to ENSEMBL to match the gene id in gene
            ## length database
            temp = suppressMessages(select(org.Hs.eg.db, rownames(rawcount),
                columns = "ENSEMBL", keytype = idtype))
            outcome = temp[temp$ENSEMBL%in%rownames(lengthDB),]
            outcome = outcome[!duplicated(outcome[,idtype]),]
            ## ind is all the one to one map
            ind = which(!outcome$ENSEMBL %in%
                outcome$ENSEMBL[duplicated(outcome$ENSEMBL)])
            outcome1 = outcome[ind,]   ## outcome1 is all the 1-1 map

            matchid = match(rownames(rawcount),outcome1[,idtype])
            genelength = lengthDB[outcome1$ENSEMBL[matchid],"bp_length"]
        }else{
            genelength = lengthDB[rownames(rawcount), "bp_length"]
        }
        
    }else{
        stopifnot(is.vector(genelength))
        stopifnot(is.numeric(genelength))
        if(any(genelength<0, na.rm = TRUE)){
            stop("genelength cannot contain negative value(s).")
        }
        if(length(genelength)!=nrow(rawcount)){
            stop("The size of genelength vector should be equal to the number 
of rows in rawcount.")
        }
    }
    
    if(any(is.na(genelength))){
        warning("Can't find gene length for ",
                round(sum(is.na(genelength))/nrow(rawcount)*100,4),
                "% genes when converting raw count to FPKM.")
    }

    rowData = data.frame(bp_length = genelength)
    se = SummarizedExperiment(assays=list(counts=rawcount), rowData=rowData)
    df = data.frame(getRPKM(se))
    df
}
