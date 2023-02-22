#' @title \code{getSampleTileMatrix}
#'
#' @description \code{getSampleTileMatrix} takes the output of peak calling with
#'   callOpenTiles and creates sample-tile matrices containing the signal
#'   intensity at each tile.
#'
#' @param tileResults a MultiAssayExperiment returned by callOpenTiles
#'   containing containing peak calling results.
#' @param cellPopulations vector of strings. Cell subsets in TileResults for
#'   which to generate sample-tile matrices. This list of group names must be
#'   identical to names that appear in the ArchRProject metadata.  If
#'   cellPopulations='ALL', then peak calling is done on all cell populations in
#'   the ArchR project metadata. Default is 'ALL'.
#' @param groupColumn Optional, the column containing sample group labels for
#'   determining consensus tiles within sample groups. Default is NULL, all
#'   samples will be used for determining consensus tiles.
#' @param threshold Threshold for consensus tiles, the minimum \% of samples
#'   (within a sample group, if groupColumn is set) that a peak must be called
#'   in to be retained. If set to 0, retain the union of all samples' peaks
#'   (this is equivalent to a threshold of 1/numSamples). It is recommended to
#'   tune this parameter to omit potentially spurious peaks.
#' @param numCores Optional, the number of cores to use with multiprocessing.
#'   Default is 1.
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#'
#' @return SampleTileMatrices a MultiAssayExperiment containing a sample-tile
#'   intensity matrix for each cell population
#'
#' @examples
#' \donttest{
#' # Starting from GRangesList
#' if (
#'   require(BSgenome.Hsapiens.UCSC.hg19) && 
#'   require(TxDb.Hsapiens.UCSC.hg38.refGene) && 
#'   require(org.Hs.eg.db)
#' ) {
#' tiles <- MOCHA::callOpenTiles(
#'   ATACFragments = MOCHA::exampleFragments,
#'   cellColData = MOCHA::exampleCellColData,
#'   blackList = MOCHA::exampleBlackList,
#'   genome = BSgenome.Hsapiens.UCSC.hg19,
#'   TxDb = TxDb.Hsapiens.UCSC.hg38.refGene,
#'   Org = org.Hs.eg.db,
#'   outDir = tempdir(),
#'   cellPopLabel = "Clusters",
#'   cellPopulations = c("C2", "C5"),
#'   numCores = 1
#' )
#' 
#' SampleTileMatrices <- MOCHA::getSampleTileMatrix(
#'   tiles,
#'   cellPopulations = c('C2', 'C5'),
#'   threshold = 0 # Take union of all samples' open tiles
#' )
#' }
#' }
#'
#' @export
getSampleTileMatrix <- function(tileResults,
                                cellPopulations = "ALL",
                                groupColumn = NULL,
                                threshold = 0.2,
                                numCores = 1,
                                verbose = FALSE) {
  if (class(tileResults)[1] != "MultiAssayExperiment") {
    stop("tileResults is not a MultiAssayExperiment")
  }

  # Any column can be used to group samples
  # Note that these are case-sensitive
  sampleData <- MultiAssayExperiment::colData(tileResults)
  validGroups <- colnames(sampleData)
  if (!is.null(groupColumn)) {
    if (!(groupColumn %in% validGroups)) {
      stop("`groupColumn` not found in the column data of tileResults.")
    }
  }

  if (length(cellPopulations) == 1 & any(tolower(cellPopulations) == "all")) {
    subTileResults <- tileResults
    cellPopulations <- names(tileResults)
  } else {
    if (all(cellPopulations %in% names(tileResults))) {
      subTileResults <- MultiAssayExperiment::subsetByAssay(tileResults, cellPopulations)
    } else {
      stop(paste(
        "All of `cellPopulations` must present in tileResults.",
        "Check `names(tileResults)` for possible cell populations."
      ))
    }
  }


  if (verbose) {
    message(stringr::str_interp("Extracting consensus tile set for each population:  ${paste(cellPopulations, collapse=', ')} "))
  }

  cl <- parallel::makeCluster(numCores)
 
  parallel::clusterExport(
      cl=cl,
      varlist=c("sampleData", "threshold", "groupColumn", "verbose"),
      envir=environment()
  )


  tilesByCellPop <- pbapply::pblapply(cl = cl, X = MultiAssayExperiment::experiments(subTileResults), 
    FUN = function(x) {

    # Get consensus tiles for this cell population for filtering
    singlePopulationConsensusTiles(
      x,
      sampleData,
      threshold,
      groupColumn,
      verbose = verbose
    )

  })

  errorMessages <- pbapply::pblapply(cl = cl, X = tilesByCellPop, 
    FUN = function(x) {
    if (any(grepl("Error", x))) {
      x[grep("Error", x)]
    } else {
      NA
    }
  })

  names(errorMessages) <- names(subTileResults)

  if (any(!is.na(errorMessages))) {
    stop(
      "Issues around thresholding and/or sample metadata. Please check user inputs, and attempt again",
      "If there are too few valid samples for a given cell type, use the variable cellPopulations to run this function on a subset of cell types, ",
      "Or, you can lower the threshold. ",
      "The following cell types were impacted:",
      paste(names(errorMessages)[!is.na(errorMessages)], collapse = ", ")
    )
  }

  allTiles <- sort(unique(do.call("c", tilesByCellPop)))

  if (verbose) {
    message(stringr::str_interp("Generating sample-tile matrix across all populations: ${paste(cellPopulations, collapse=', ')} "))
  }

  parallel::clusterExport(
      cl=cl,
      varlist=c("allTiles"),
      envir=environment()
  )

  # consensusTiles is used to  extract rows (tiles) from this matrix
  sampleTileIntensityMatList <- pbapply::pblapply(cl = cl,
   X = MultiAssayExperiment::experiments(tileResults), 
   FUN = function(x) {

    singlePopulationSampleTileMatrix(
      x,
      allTiles,
      NAtoZero = FALSE
    )

  })


  #Stop Clusters and run garbage collection. 
  parallel::stopCluster(cl)
  gc()

  # Order sampleData rows to match the same order as the columns
  maxMat <- which.max(lapply(sampleTileIntensityMatList, ncol))
  colOrder <- colnames(sampleTileIntensityMatList[[maxMat]])
  sampleData <- sampleData[match(colOrder, rownames(sampleData)), ]

  . <- NULL
  
  tilePresence <- lapply(tilesByCellPop, function(x) (allTiles %in% x)) %>%
    do.call("cbind", .) %>%
    as.data.frame()
  
  allTilesGR <- MOCHA::StringsToGRanges(allTiles)
  GenomicRanges::mcols(allTilesGR) <- tilePresence

  results <- SummarizedExperiment::SummarizedExperiment(
    sampleTileIntensityMatList,
    rowRanges = allTilesGR,
    colData = sampleData,
    metadata = MultiAssayExperiment::metadata(tileResults)
  )
  return(results)
}
