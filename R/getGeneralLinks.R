#' @title \code{getGeneralLinks}
#'
#' @description \code{getGeneralLinks} takes aset of genes, pseudobulk expression and accessibility, and finds region accessibility to gene expression correlations.
#'
#' @param Obj1 The SummarizedExperiment object output from getSampleTileMatrix containing your sample-tile matrices
#' @param SampleTileObj The SummarizedExperiment object output from getSampleTileMatrix containing your sample-tile matrices
#' @param cellPopulation A string denoting the cell population of interest, which must be present in SampleTileObj
#' @param Obj1List A list of gene names to test. Must aligned with gene names from the Obj2. Should be differential.
#' @param Obj2List a list of differential ChromVAR Z-scores to test. 
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
getGeneralLinks <- function(Obj1,
                          Obj2,
                          Obj1List = NULL,
                          Obj2List = NULL,
                          returnBackground = FALSE,
                          ZI_inflated = FALSE,
                          numCores = 20,
                          verbose = TRUE) {

    
    . <- NULL

    if(length(SummarizedExperiment::assays(Obj2)) > 1){
        stop("Obj2 needs to be subsetted to only one assay.")
    }

    if(length(SummarizedExperiment::assays(Obj1)) > 1){
        stop("Obj1 needs to be subsetted to only one assay.")
    }

    if(!all(colnames(Obj1) %in% colnames(Obj2))){
        stop("Samples are not aligned between Obj1 and Obj2.")
    }

    if(!all(Obj2List %in% rownames(Obj2))){
        stop("Obj2List not all found in Obj2.")
    }

    if(!all(Obj1List %in% rownames(Obj1))){
        stop("Obj1List not all found in Obj1.")
    }

    if(length(Obj2List)*2 > length(rownames(Obj2))){
        stop("Obj2List is too long. Provide an Obj2List that is less than half of the variables in Obj2.") 
    }

    if(length(Obj1List)*2 > length(rownames(Obj1))){
       stop("Obj1List is too long. Provide an Obj1List that is less than half of the variables in Obj1.") 
   }

    ## Extract the matrices for both. 
    mat1 <- SummarizedExperiment::assays(Obj1)[[1]]
    mat2 <- SummarizedExperiment::assays(Obj2)[[1]]

    #Align the matrices.
    mat2 <- mat2[,colnames(mat1)]

    ## Find all tiles near each gene, and iterate over each one for ZI-spearman correlations. 
    iterList <- pbapply::pblapply(cl = NULL, X = seq_along(Obj2List), function(x){

        #Extract a data.frame of all the tile intensities
        submat1 <- mat1[rownames(mat1) %in% Obj2List[x],]

        #Extract a data.frame of gene expression for the matching gene. 
        submat2 <- mat2[rownames(mat2) %in% Obj1List,]
        if(dim(submat2)[1] > 1){ stop('One gene name is duplicated in your gene expression matrix.')}

        list(submat1, submat2)
    })

    ## Now iterate over the list and run a Zero-inflated spearman for all combinations of tiles and genes. 
    cl <- parallel::makeCluster(numCores)
    if(ZI_inflated){
        foregroundDF <- pbapply::pblapply(cl = cl, X = iterList, tileGeneCorrelations)
    }else{
        foregroundDF <- pbapply::pblapply(cl = cl, X = iterList, generalCorrelations)
    }
    foregroundDF <- do.call('rbind', foregroundDF)
    parallel::stopCluster(cl)

    rm(iterList)

    ## Generate a background set of tiles and genes (non-DEGs) by finding a background set of tile-gene pairs.
    backMat2 <- sample(rownames(mat2)[! rownames(mat2) %in% Obj2List], length(Obj2List), replace = FALSE)
    backMat1 <- sample(rownames(mat1)[! rownames(mat1) %in% Obj1List], length(Obj1List), replace = FALSE)
    
    iterList <- pbapply::pblapply(cl = NULL, X = seq_along(backMat1), function(x){

        #Extract a data.frame of all the tile intensities
        submat1 <- mat1[rownames(mat1) %in% backMat1[x],]

        #Extract a data.frame of gene expression for the matching gene. 
        submat2 <- mat2[rownames(mat2) %in% backMat2,]
        if(dim(submat2)[1] > 1){ stop('One gene name is duplicated in your gene expression matrix.')}

        list(submat1, submat2)
    })

    ## Now iterate over the list and run a Zero-inflated spearman for all combinations of tiles and genes. 
    cl <- parallel::makeCluster(numCores)
    if(ZI_inflated){
        backgroundDF <- pbapply::pblapply(cl = cl, X = iterList, tileGeneCorrelations)

    }else{
        backgroundDF <- pbapply::pblapply(cl = cl, X = iterList, generalCorrelations)
    }
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


generalCorrelations <- function(iterList){

    mat1 <- iterList[[1]]
    mat2 <- iterList[[2]]

    spearmanCorr <- wCorr::weightedCorr(x = x, y = y, weights = w, method = "Spearman")
    te_corr <- unlist(lapply(1:nrows(mat2), function(x){
        wCorr::weightedCorr(mat1[1,],mat2[x,],method = "Spearman")
    }))

    df <- data.frame(Obj1 = rep(rownames(mat1), length(te_corr)),
                     Obj2 = rownames(mat2),
                     Correlations = te_corr)
    return(df)
}
