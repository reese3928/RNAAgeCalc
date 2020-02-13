

#' @title Make plot to visualize RNA age
#'
#' @description This function makes plots to visualize the relationship
#' between chronological age and RNA age.
#' @param res a data frame returned by `predict_age` function. If the
#' chronological age is not provided when using `predict_age` function,
#' visulization cannot be made.
#' @param main title of the plot
#' @param xlab label of x-axis
#' @param ylab label of y-axis
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth ggtitle xlab ylab
#' theme_bw theme element_text
#' @export
#' @return the plot which shows RNA age vs chronological age
#' @examples
#' data(fpkmExample)
#' fpkm_large = cbind(fpkm, fpkm, fpkm, fpkm)
#' fpkm_large = cbind(fpkm_large, fpkm_large, fpkm_large, fpkm_large)
#' colnames(fpkm_large) = paste0("sample",1:32)
#' chronage = data.frame(sampleid = colnames(fpkm_large), age = 1:32)
#' res = predict_age(exprdata = fpkm_large, exprtype = "FPKM",
#' chronage = chronage)
#' makeplot(res)

makeplot = function(res, main = "RNA age vs chronological age",
                    xlab = "chronological age", ylab = "RNA Age"){
    if(!"ChronAge"%in%colnames(res)){
        stop("Chronological age is not found in res data frame.")
    }
    if(nrow(na.omit(res))<30){
        warning("Less than 30 samples. The linear regression on RNA age vs
chronological age may not be reliable.")
    }
    ggplot(res, aes(x=ChronAge, y=RNAAge)) + geom_point(shape=19,size=1) +
        geom_smooth(method=lm) + ggtitle(main) + xlab(xlab) + ylab(ylab) +
        theme_bw() + theme(plot.title = element_text(hjust = 0.5))
}


