#' @title \code{getTileGeneLinks}
#'
#' @description \code{getTileGeneLinks} takes aset of genes, pseudobulk expression and accessibility, and finds region accessibility to gene expression correlations.
#'
#'
#' @param SampleGeneObj The SummarizedExperiment object output from getSampleGemeMatrix containing pseudobulked scRNA by cell type. 
#' @param cellPopulation A string denoting the cell population of interest, which must be present in SampleTileObj
#' @param gene_names A list of gene names to test. Must aligned with gene names from the SampleGeneObj. Should be differential.
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
#' 
#' Thoughts: Test interaction between a region and variable (or direct)

getTileGeneLinks <- function(SampleTileObj,
                          SampleGeneObj,
                          cellPopulation = "All",
                          gene_names = NULL,
                          windowSize = 1*10^6,
                          backgroundSize = NULL,
                          returnBackground = FALSE,
                          numCores = 20,
                          verbose = TRUE) {

    
    . <- NULL

    if(!all(gene_names %in% rownames(SampleGeneObj))){
        stop('Gene names are not all found in the Sample Gene Object.')
    }

    if(!all(cellPopulation %in% names(SummarizedExperiment::assays(SampleTileObj)))){
        stop("Cell populations not found within SampleTileObj.")
    }else{
        subSTM <- subsetMOCHAObject(SampleTileObj,  susbetBy = 'celltype', groupList = cellPopulation, subsetPeaks = TRUE)
        fullSTM <- combineSampleTileMatrix(subSTM)
    }

    if(!all(cellPopulation %in% names(SummarizedExperiment::assays(SampleGeneObj)))){
        stop("Cell populations not found within SampleGeneObj.")
    }else{
        subSGM <- subsetSampleGeneObject(SampleGeneObj,  susbetBy = 'celltype', groupList = cellPopulation)
        fullSGM <- combineSampleGeneObject(subSGM)
    }

    ## Extract the matrices for both. 
    accMat <- SummarizedExperiment::assays(fullSTM)[[1]]
    exprMat <- SummarizedExperiment::assays(fullSGM)[[1]]

    if(!all(colnames(accMat) %in% colnames(exprMat))){
        stop('Samples and/or cell type names are not the same between the SampleTileObject and SampleGeneObject. Please align sample and cell type names.')
    }

    #Align the matrices.
    exprMat <- exprMat[,colnames(accMat)]

    ## Extract all promoter tiles within the SampleTileObj. 
    allTiles <- SummarizedExperiment::rowRanges(SampleTileObj)
    if(!all(c('Gene', 'tileType') %in% colnames(SummarizedExperiment::mcols(allTiles)))){
        stop('Tiles within the SampleTileObj are not annotated. Please run annotateTiles() and try again.')
    }

    promotersOfInterest <- getDEGPromoters(allTiles, gene_names)

    ## Now take the promoters of interest, and expand it by the windowSize either direction.
    ## Reduce those ranges by gene so that promoters and TSS related to the same gene are compressed into one GRanges. 
    promoterWindow <- plyranges::stretch(plyranges::anchor_center(promotersOfInterest), extend = windowSize*2)
    promoterWindow <- plyranges::reduce_ranges(plyranges::group_by(promoterWindow, Genes))
    
    ## Find all tiles near each gene, and iterate over each one for ZI-spearman correlations. 
    iterList <- pbapply::pblapply(cl = NULL, X = seq_along(promoterWindow), function(x){

        #Find all tiles within the windowSize either direction of the promoter of interest. 
        overlapTiles <- plyranges::filter_by_overlaps(allTiles, promoterWindow[x])


        #Extract a data.frame of all the tile intensities
        subAccMat <- accMat[rownames(accMat) %in% GRangesToString(overlapTiles),]

        #Extract a data.frame of gene expression for the matching gene. 
        subExprMat <- exprMat[rownames(exprMat) %in% promoterWindow$Genes[x],]
        if(dim(subExprMat)[1] > 1){ stop('One gene name is duplicated in your gene expression matrix.')}

        list(subAccMat, subExprMat)
    })

    ## Now iterate over the list and run a Zero-inflated spearman for all combinations of tiles and genes. 
    cl <- parallel::makeCluster(numCores)
    foregroundDF <- pbapply::pblapply(cl = cl, X = iterList, tileGeneCorrelations)
    foregroundDF <- do.call('rbind', foregroundDF)
    parallel::stopCluster(cl)

    ## Generate a background set of tiles and genes (non-DEGs) by finding a background set of tile-gene pairs.
    backGenes <- sample(rownames(exprMat)[! rownames(exprMat) %in% gene_names], length(gene_names), replace = FALSE)
    
    if(length(backGenes) == 0){
        stop('No background geneset could be found.')
    }

    if(is.null(backgroundSize)){

        backTiles <- sample(rownames(accMat), length(unique(foreGroundDF$Tiles)), replace = FALSE)

    } else {

        backTiles <- sample(rownames(accMat), backgroundSize, replace = FALSE)

    }
    
    tilesPerGene <- length(backTiles) %/% length(backGenes)
    
    splitTileList <- split(backTiles, rep(1:tilesPerGene,length(backGenes)))

    iterList <- lapply(seq_along(splitTileList), function(x){

        subAccMat <- accMat[rownames(accMat) %in% splitTileList[[x]],]

        #Extract a data.frame of gene expression for the matching gene. 
        subExprMat <- exprMat[rownames(exprMat) %in% backGenes[x],]
        if(dim(subExprMat)[1] > 1){ stop('One gene name is duplicated in your gene expression matrix.')}

        list(subAccMat, subExprMat)

    })

    ## Now iterate over the list and run a Zero-inflated spearman for all combinations of tiles and genes. 
    cl <- parallel::makeCluster(numCores)
    backgroundDF <- pbapply::pblapply(cl = cl, X = iterList, tileGeneCorrelations)
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

    ## Annotate tile type. 
    foreground$tileType <- allTiles$tileType[match(StringsToGRanges(foreground$Tiles), GRangesToString(allTiles))]
    foreground$tileGenes <- allTiles$Gene[match(StringsToGRanges(foreground$Tiles), GRangesToString(allTiles))]

    if(returnBackground){

        returnList <- list(foregroundDF, backgroundDF)
        names(returnList) <- c('Foreground', 'Background')
        return(returnList)
        
    }else{

        return(foregroundDF)
    }

}


tileGeneCorrelations <- function(iterList){

    accMat <- iterList[[1]]
    exprMat <- iterList[[2]]

    te_corr <- unlist(lapply(1:nrows(accMat), function(x){
        weightedZISpearman(exprMat[1,], accMat[x,], ZI = TRUE)
    }))

    df <- data.frame(Tiles = rownames(accMat), 
                     Genes = rep(rownames(exprMat), length(te_corr)),
                     Correlations = te_corr)
    return(df)
}


getDEGPromoters <- function(tileGR, gene_names = NULL){
    
    if(!all(c('Gene', 'tileType') %in% colnames(GenomicRanges::mcols(tileGR)))){
        stop('Tiles within the SampleTileObj are not annotated. Please run annotateTiles() and try again.')
    }

    if(is.null(gene_names)){
        stop('No gene_names provided.')
    }

    ## Identify all promoter tiles related to the genes within the DEG List
    promoterTiles <- plyranges::filter(tileGR, tileType == 'Promoter')
    geneDF <-  do.call('rbind',lapply(seq_along(promoterTiles$Gene), function(x){
            geneList <- unlist(stringr::str_split(.,", "))
            data.frame(index = x, genes = geneList)
            }))
    percentOpen <- sum(gene_names %in% geneDF$genes)/length(gene_names)*100

    messages(stringr::str_interp("{percentOpen}% of the gene_names have an accessible promoter region."))

    subGeneDF <- dplyr::filter(geneDF, genes %in% gene_names) %>% distinct()
    promotersOfInterest <- promoterTiles[subGeneDF$index]
    mcols(promotersOfInterest) <- NULL
    promotersOfInterest$Genes = subGeneDF$genes

    return(promotersOfInterest)

}