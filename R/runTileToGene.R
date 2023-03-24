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
#' @param approximateTile If set to TRUE, it will use all tiles that overlap with the regions given, instead of finding an exact match to the regions variable. Default is FALSE.
#' @param ZI boolean flag that enables zero-inflated (ZI) Spearman correlations to be used. Default is TRUE. If FALSE, skip zero-inflation and calculate the normal Spearman.
#'
#' @return TileCorr A data.table correlation matrix
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
runTileToGene <- function(SampleTileObj,
                          SampleGeneObj,
                          cellPopulation = "All",
                          DEGList = NULL,
                          windowSize = 1*10^6) {

    
    . <- NULL

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