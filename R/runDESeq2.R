# run deseq2 from psudobulk data. Written by Ziyuan He
#' @title \code{runDEseq2}
#'
#' @description \code{runDEseq2} a wrapper function to run deseq2 from psudobulk sce object
#'
#' @param pseudo a pseudobulk SingleCellExperiment Object generated from MakeSoPseudoBulk function
#' @param design a string of formula of expermental design
#'
#' @return dds deseq2 object

runDEseq2 <- function(pseudo, design){
    # load library
    require('SummarizedExperiment')
    require("DESeq2")
    se <- SummarizedExperiment(assays=list(counts=counts(pseudo)),
                         colData= as.data.frame(colData(pseudo)))
    se$sizeFactor <- NULL
    dds <- DESeqDataSet(se, design = as.formula(design))
    dds <- DESeq(dds)
    return(dds)
}