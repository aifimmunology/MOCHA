#' @title \code{getCoAccessibleLinks}
#'
#' @description \code{getCoAccessibleLinks} takes an input set of regions (tiles) and finds co-accessible neighboring regions within a window. Co-accessibility is defined as the correlation between two region intensity (openness) across samples.
#'
#'
#' @param SampleTileObj The SummarizedExperiment object output from getSampleTileMatrix containing your sample-tile matrices
#' @param cellPopulation A string denoting the cell population of interest, which must be present in SampleTileObj
#' @param regions a GRanges object or vector or strings containing the regions on which to compute co-accessible links. Strings must be in the format "chr:start-end", e.g. "chr4:1300-2222".
#'   Can be the output from getDifferentialAccessibleTiles.
#' @param chrChunks This functions subsets by groups of chromosome, and then parallelizes within each group of chromosomes when running correlations. This method keeps memory
#'   low. To speed things up on high performing platforms, you can chunk out more than one chromosome at a time. Default is chrChunks = 1, so only one chromosome at a time.
#' @param windowSize the size of the window, in basepairs, around each input region to search for co-accessible links
#' @param numCores Optional, the number of cores to use with multiprocessing. Default is 1.
#' @param verbose Set TRUE to display additional messages. Default is FALSE.

#' @return TileGeneCorr A data.table correlation matrix
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
runTileToGeneLinks <- function(SampleTileObj,
                          SampleGeneObj,
                          cellPopulation = "All",
                          DEGList = NULL,
                          windowSize = 1*10^6,
                          numCores = 20) {

    
    . <- NULL

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


    ## Extract all promoter tiles within the SampleTileObj. 
    allTiles <- SummarizedExperiment::rowRanges(SampleTileObj)
    if(!all(c('Gene', 'tileType') %in% colnames(SummarizedExperiment::mcols(allTiles)))){
        stop('Tiles within the SampleTileObj are not annotated. Please run annotateTiles() and try again.')
    }

    promotersOfInterest <- getDEGPromoters(allTiles, DEGList)

    ## Now take the promoters of interest, and expand it by the windowSize either direction.
    ## Reduce those ranges by gene so that promoters and TSS related to the same gene are compressed into one GRanges. 
    promoterWindow <- plyranges::stretch(plyranges::anchor_center(promotersOfInterest), extend = windowSize*2)
    promoterWindow <- plyranges::reduce_ranges(plyranges::group_by(promoterWindow, Genes))
    
    ## Find all tiles near each gene, and iterate over each one for ZI-spearman correlations. 
    iterList <- pbapply::pblapply(cl = NULL, X = seq_along(promoterWindow), function(x){

        #Find all tiles within the windowSize either direction of the promoter of interest. 
        overlapTiles <- plyranges::filter_by_overlaps(allTiles, promoterWindow[x])


        #Extract a data.frame of all the tile intensities
        subAccMat <- accessibleMatrix[rownames(accessibleMatrix) %in% GRangesToString(overlapTiles),]

        #Extract a data.frame of gene expression for the matching gene. 
        subExprMat <- expressionMatrix[rownames(expressionMatrix) %in% promoterWindow$Genes[x],]
        if(dim(subExprMat)[1] > 1){ stop('One gene name is duplicated in your gene expression matrix.')}
        list(subAccMat, subExprMat)
    })

    ## Now iterate over the list and run a Zero-inflated spearman for all combinations of tiles and genes. 
    cl <- parallel::makeCluster(numCores)
    foregroundDF <- pbapply::pblapply(cl = cl, X = iterList, tileGeneCorrelations)
    foregroundDF <- do.call('rbind', foregroundDF)
    parallel::stopCluster(cl)

    ## Generate a background set of tiles and genes (non-DEGs) by finding a background set of tile-gene pairs.
    backGenes <- sample(rownames(exprMat)[! rownames(exprMat) %in% DEGList], length(DEGList), replace = FALSE)
    if(length(backGenes) == 0){
        stop('No background geneset could be found.')
    }
    backTiles <- sample(rownames(accMat), length(unique(foreGroundDF$Tiles)), replace = FALSE)
    
    tilesPerGene <- length(backTiles) %/% length(backGenes)
    
    back_df <- data.frame(Tiles = unlist(lapply(tilesPerGene, function(x){



                        })),
                             Genes = )

    iterList <- lapply(backGenes, function(x){

        

    })

        ## Now iterate over the list and run a Zero-inflated spearman for all combinations of tiles and genes. 
    cl <- parallel::makeCluster(numCores)
    backgroundDF <- pbapply::pblapply(cl = cl, X = iterList, ZISpearman)
    backgroundDF <- do.call('rbind', backgroundDF)
    parallel::stopCluster(cl)


    ## Test that set. 
    if (length(tile1) != length(tile2)) {
    stop("tile1 and tile2 must be the same length.")
    }

    fullObj <- combineSampleTileMatrix(SampleTileMatrix)

    backPeaks <- chromVAR::getBackgroundPeaks(fullObj)

    if (is.character(tile1) && is.character(tile2)) {
    nTile1 <- match(tile1, rownames(fullObj))
    nTile2 <- match(tile2, rownames(fullObj))
    } else if (is.numeric(tile1) && is.numeric(tile2)) {
    nTile1 <- tile1
    nTile2 <- tile2

    tile1 <- rownames(fullObj)[nTile1]
    tile2 <- rownames(fullObj)[nTile2]
    } else {
    stop("tile1 and tile 2 must both be either numbers (indices) or strings")
    }

    if (backNumber > 2450) {
    backNumber <- 2450
    if (verbose) {
        warning("backNumber too high, setting to maximum of 2450.")
    }
    } else if (backNumber <= 10) {
    stop("backNumber too low (<=10). We recommend at least 100.")
    }

    cl <- parallel::makeCluster(numCores)
    parallel::clusterExport(cl, varlist = c("nTile1", "nTile2", "backPeaks"), envir = environment())

    if (verbose) {
    message("Finding background peak pairs")
    }

    backgroundCombos <- pbapply::pblapply(seq_along(nTile1), function(x) {
    tmpMat <- expand.grid(backPeaks[nTile1[x], ], backPeaks[nTile2[x], ])
    tmpMat <- tmpMat[tmpMat[, 1] != tmpMat[, 2], ]
    tmpMat[sample.int(dim(tmpMat)[1], backNumber), ]
    }, cl = cl)





}

tileGeneCorrelations <- function(iterList){

    accMat <- iterList[[1]]
    exprMat <- iterList[[2]]

    te_corr <- unlist(lapply(1:nrows(accMat){
        weightedZISpearman(exprMat[1,], accMat[x,], ZI = TRUE)
    }))

    df <- data.frame(Tiles = rownames(accMat), 
                     Genes = rep(rownames(exprMat), length(te_corr)),
                     Correlations = te_corr)
    return(df)
}


getDEGPromoters <- function(tileGR, DEGList = NULL){
    
    if(!all(c('Gene', 'tileType') %in% colnames(GenomicRanges::mcols(tileGR)))){
        stop('Tiles within the SampleTileObj are not annotated. Please run annotateTiles() and try again.')
    }

    if(is.null(DEGList)){
        stop('No DEGList provided.')
    }

    ## Identify all promoter tiles related to the genes within the DEG List
    promoterTiles <- plyranges::filter(tileGR, tileType == 'Promoter')
    geneDF <-  do.call('rbind',lapply(seq_along(promoterTiles$Gene), function(x){
            geneList <- unlist(stringr::str_split(.,", "))
            data.frame(index = x, genes = geneList)
            }))
    percentOpen <- sum(DEGList %in% geneDF$genes)/length(DEGList)*100

    messages(stringr::str_interp("{percentOpen}% of the DEGList have an accessible promoter region."))

    subGeneDF <- dplyr::filter(geneDF, genes %in% DEGList) %>% distinct()
    promotersOfInterest <- promoterTiles[subGeneDF$index]
    mcols(promotersOfInterest) <- NULL
    promotersOfInterest$Genes = subGeneDF$genes

    return(promotersOfInterest)

}