#' @title \code{getGeneMotifLinks}
#'
#' @description \code{getGeneMotifLinks} takes aset of genes, pseudobulk expression and accessibility, and finds region accessibility to gene expression correlations.
#'
#' @param ChromVARObj The SummarizedExperiment object output from getSampleTileMatrix containing your sample-tile matrices
#' @param SampleTileObj The SummarizedExperiment object output from getSampleTileMatrix containing your sample-tile matrices
#' @param cellPopulation A string denoting the cell population of interest, which must be present in SampleTileObj
#' @param DEGList A list of gene names to test. Must aligned with gene names from the SampleGeneObj. Should be differential.
#' @param DUMList a list of differential ChromVAR Z-scores to test. 
#' @param windowSize the size of the window, in basepairs, around each input region to search for co-accessible links
#' @param backgroundSize Integer. the size of background tileset to use. If NULL, the background tileset will be sized matched to the foreground set.
#' @param returnBackground Boolean. Set TRUE to return a list of both foreground and background correlations. 
#' @param numCores Optional, the number of cores to use with multiprocessing. Default is 1.
#' @param verbose Set TRUE to display additional messages. Default is FALSE.

#' @return foregroundDF A data.table correlation matrix
#'
#' @details The technical details of the zero-inflated correlation can be
#'          found here:
#'
#'               Pimentel, Ronald Silva, "Kendall's Tau and Spearman's Rho
#'               for Zero-Inflated Data" (2009). Dissertations.
#'
#'          while the implementation (scHOT R package), can be found here:
#'               http://www.bioconductor.org/packages/release/bioc/html/scHOT.html
#'
#'
#' @export
getGeneMotifLinks <- function(ChromVARObj,
                          SampleGeneObj,
                          DEGList = NULL,
                          DUMList = NULL,
                          returnBackground = FALSE,
                          numCores = 20,
                          verbose = TRUE) {

    
    . <- NULL

    if(length(SummarizedExperiment::assays(SampleGeneObj)) > 1){
        stop("SampleGeneObj needs to be subsetted to only one cell type.")
    }

    if(!all(colnames(ChromVARObj) %in% colnames(SampleGeneObj))){
        stop("Samples are not aligned between ChromVARObj and SampleGeneObj.")
    }

    if(!all(DEGList %in% rownames(SampleGeneObj))){
        stop("Not all gene names are found in SampleGeneObj.")
    }

    if(!all(DUMList %in% rownames(ChromVARObj))){
        stop("Motif names not found in ChromVARObj.")
    }

    if(length(DUMList)*2 > length(rownames(ChromVARObj))){
        stop("More than half of motifs are listed as Differentially Used. This is unlikely to be accurate. Provide a shorter DUMList.")
    }

    if(length(DEGList)*2 > length(rownames(SampleGeneObj))){
        stop("More than half of genes are listed as Differentially Expressed. This is unlikely to be accurate. Provide a shorter DEGList.")
    }

    ## Extract the matrices for both. 
    motifMat <- SummarizedExperiment::assays(ChromVARObj)$z
    exprMat <- SummarizedExperiment::assays(SampleGeneObj)[[1]]

    #Align the matrices.
    exprMat <- exprMat[,colnames(motifMat)]

    ## Find all tiles near each gene, and iterate over each one for ZI-spearman correlations. 
    iterList <- pbapply::pblapply(cl = NULL, X = seq_along(DUMList), function(x){

        #Extract a data.frame of all the tile intensities
        submotifMat <- motifMat[rownames(motifMat) %in% DUMList[x],]

        #Extract a data.frame of gene expression for the matching gene. 
        subExprMat <- exprMat[rownames(exprMat) %in% DEGList,]
        if(dim(subExprMat)[1] > 1){ stop('One gene name is duplicated in your gene expression matrix.')}

        list(submotifMat, subExprMat)
    })

    ## Now iterate over the list and run a Zero-inflated spearman for all combinations of tiles and genes. 
    cl <- parallel::makeCluster(numCores)
    foregroundDF <- pbapply::pblapply(cl = cl, X = iterList, motifGeneCorrelations)
    foregroundDF <- do.call('rbind', foregroundDF)
    parallel::stopCluster(cl)

    rm(iterList)

    ## Generate a background set of tiles and genes (non-DEGs) by finding a background set of tile-gene pairs.
    backGenes <- sample(rownames(exprMat)[! rownames(exprMat) %in% DEGList], length(DEGList), replace = FALSE)
    backMotifs <- sample(rownames(motifMat)[! rownames(motifMat) %in% DUMList], length(DUMList), replace = FALSE)
    
    iterList <- pbapply::pblapply(cl = NULL, X = seq_along(backMotifs), function(x){

        #Extract a data.frame of all the tile intensities
        submotifMat <- motifMat[rownames(motifMat) %in% backMotifs[x],]

        #Extract a data.frame of gene expression for the matching gene. 
        subExprMat <- exprMat[rownames(exprMat) %in% backGenes,]
        if(dim(subExprMat)[1] > 1){ stop('One gene name is duplicated in your gene expression matrix.')}

        list(submotifMat, subExprMat)
    })

    ## Now iterate over the list and run a Zero-inflated spearman for all combinations of tiles and genes. 
    cl <- parallel::makeCluster(numCores)
    backgroundDF <- pbapply::pblapply(cl = cl, X = iterList, motifGeneCorrelations)
    backgroundDF <- do.call('rbind', backgroundDF)
    parallel::stopCluster(cl)
  
    if (verbose) {
      message("Generating p-values.")
    }

    greatList <- unlist(pbapply::pblapply(foregroundDF$Correlation[which(foregroundDF$Correlation > 0)],
        function(x){
                  return(sum(x > backgroundDF$Correlation))
          }, cl = NULL))/length(backgroundDF$Correlation)

    lesserList <- unlist(pbapply::pblapply(foregroundDF$Correlation[which(foregroundDF$Correlation < 0)], 
          function(x){
                  return(sum(x < backgroundDF$Correlation))
          }, cl = NULL))/length(backgroundDF$Correlation)


    foregroundDF$pValues <- rep(NA, length(foregroundDF$Correlation))
    foregroundDF$pValues[which(foregroundDF$Correlation > 0)] = 1-greatList
    foregroundDF$pValues[which(foregroundDF$Correlation < 0)] = 1-lesserList

    foregroundDF$FDR <- p.adjust(foregroundDF$pValues, method = 'fdr')

    if(returnBackground){

        returnList <- list(foregroundDF, backgroundDF)
        names(returnList) <- c('Foreground', 'Background')
        return(returnList)
        
    }else{

        return(foregroundDF)
    }

}


motifGeneCorrelations <- function(iterList){

    motifMat <- iterList[[1]]
    exprMat <- iterList[[2]]

    spearmanCorr <- wCorr::weightedCorr(x = x, y = y, weights = w, method = "Spearman")
    te_corr <- unlist(lapply(1:nrows(exprMat), function(x){
        wCorr::weightedCorr(motifMat[1,],exprMat[x,],method = "Spearman")
    }))

    df <- data.frame(Motifs = rep(rownames(motifMat), length(te_corr)),
                     Genes = rownames(exprMat),
                     Correlations = te_corr)
    return(df)
}
