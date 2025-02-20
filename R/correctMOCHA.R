#' @title \code{correctGenome}
#'
#' @description \code{correctGenome} Runs ChromVAR on a MOCHA Tile-Sample Accessibility Matrix object
#'
#'
#' @param TSAM RangedSummarizedExperiment object, containing the TSAM from getSampleTileMatrix, with a motif set added via addMotifSet()
#' @param genome genome to correct with. 
#'
#' @return chromVAR object
#'
#' @details This is a wrapper for basic SummarizedExperiment-based subsetting. 
#'
#'
#' @export
#' 
#' 


correctGenome<- function(TSAM, genome){

    S4Vectors::metadata(TSAM)$Genome <- genome
    return(TSAM)
}
