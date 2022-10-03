#' @title \code{getSampleTileMatrix}
#'
#' @description \code{getSampleTileMatrix} takes the output of peak calling with
#'   callOpenTiles and creates sample-tile matrices containing the signal intensity
#'   at each tile.
#'
#'
#' @param tileResults a MultiAssayExperiment returned by callOpenTiles containing
#'   containing peak calling results.
#' @param cellPopulations vector of strings. Cell subsets in TileResults for which to generate sample-tile matrices. This list of group names must be identical to names that appear in
#'   the ArchRProject metadata.  If cellPopulations='ALL', then peak
#'   calling is done on all cell populations in the ArchR project metadata.
#'   Default is 'ALL'.
#' @param groupColumn Optional, the column containing sample group labels for determining consensus tiles within sample groups. Default is NULL, all samples will be used for determining consensus tiles.
#' @param threshold Threshold for consensus tiles, the minimum % of samples (within a sample group) that a peak must be called in to be retained.
#' @param log2Intensity Boolean, set to TRUE to return the log2 of the sample-tile intensity matrix. Optional, default is FALSE.
#' @param numCores Optional, the number of cores to use with multiprocessing. Default is 1.
#'
#' @return SampleTileMatrices a MultiAssayExperiment containing a sample-tile intensity matrix
#'   for each cell population
#'
#' @details The technical details of the algorithm are found in XX.
#'
#' @references XX
#'
#' @export

getSampleTileMatrix <- function(tileResults,
                                cellPopulations = "ALL",
                                groupColumn = NULL,
                                threshold = 0.2,
                                log2Intensity = TRUE,
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

  if (length(cellPopulations) == 1 & any(cellPopulations == "ALL")) {
    subTileResults <- tileResults
    cellPopulations <- names(tileResults)
  } else {
    if (all(cellPopulations %in% names(tileResults))) {
      subTileResults <- tileResults[names(tileResults) %in% cellPopulations]
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

  tilesByCellPop <- parallel::mclapply(MultiAssayExperiment::experiments(subTileResults), function(x) {

    # Get consensus tiles for this cell population for filtering
    scMACS:::singlePopulationConsensusTiles(
      x,
      sampleData,
      threshold,
      groupColumn
    )
  }, mc.cores = numCores)

  allTiles <- sort(unique(do.call("c", tilesByCellPop)))

  if(any(grepl('Error', allTiles))){

    print(allTiles[grep('Error', allTiles)])
    stop('Issues around thresholding and/or sample metadata. Please check user inputs, and attempt again',
          'If there are too few valid samples for a given cell type, use the variable cellPopulations to run this function on a subset of cell types',
          'Or, you can lower the threshold. ')


  }

  if (verbose) {
    message(stringr::str_interp("Generating sample-tile matrix across all populations: ${paste(cellPopulations, collapse=', ')} "))
  }

  # consensusTiles is used to  extract rows (tiles) from this matrix
  sampleTileIntensityMatList <- parallel::mclapply(MultiAssayExperiment::experiments(tileResults), function(x) {
    scMACS:::singlePopulationSampleTileMatrix(
      x,
      allTiles,
      NAtoZero = FALSE,
      log2Intensity
    )
  }, mc.cores = numCores)
    
  # Order sampleData rows to match the same order as the columns
  maxMat <- which.max(lapply(sampleTileIntensityMatList, ncol))
  colOrder <- colnames(sampleTileIntensityMatList[[maxMat]])
  sampleData <- sampleData[match(colOrder, rownames(sampleData)),]
    
  tilePresence <- lapply(tilesByCellPop, function(x) (allTiles %in% x)) %>%
    do.call("cbind", .) %>%
    as.data.frame()
  allTilesGR <- scMACS::StringsToGRanges(allTiles)
  mcols(allTilesGR) <- tilePresence

  results <- SummarizedExperiment::SummarizedExperiment(
    sampleTileIntensityMatList,
    rowRanges = allTilesGR,
    colData = sampleData,
    metadata = append(
      list("Log2Intensity" = log2Intensity),
      MultiAssayExperiment::metadata(tileResults)
    )
  )
  return(results)
}
