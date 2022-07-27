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
#' @param join The join method to combine consensus tiles across sample groups. Can be "union" (default) or "intersect". 
#' @param NAtoZero Boolean, set to TRUE to convert NA intensities from peak calling (tiles with too low signal) to zeroes. Optional, default is FALSE.
#' @param log2Intensity Boolean, set to TRUE to return the log2 of the sample-tile intensity matrix. Optional, default is FALSE.
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
                                join = "union",
                                NAtoZero = FALSE,
                                log2Intensity = TRUE
                               ) {

  if (!(join %in% c("union", "intersect"))){
    stop("`join` must be either 'union' or 'intersect'")
  }
  
  if(class(tileResults)[1] != 'MultiAssayExperiment'){
    stop('tileResults is not a MultiAssayExperiment')  
  }
  
  # Any column can be used to group samples
  # Note that these are case-sensitive
  sampleData <- colData(tileResults)
  validGroups <- colnames(sampleData)
  if (!is.null(groupColumn)){
    if (!(groupColumn %in% validGroups)){
      stop("`groupColumn` not found in the column data of tileResults.")
    }
  }
  
  for (cellPop in cellPopulations){
      if (!(cellPop %in% names(tileResults))){
        stop(paste(
          "All of `cellPopulations` must present in tileResults.",
          "Check `names(tileResults)` for possible cell populations."
        ))
      }
  }
  
  if(length(cellPopulations) == 1){
    if (cellPopulations == 'ALL'){
      cellPopulations <- names(tileResults)
    }
  }
  
  experimentList <- list()
  idx <- 1
  for (cellPopulation in cellPopulations){
    
    message(str_interp("Extracting sampleTileMatrix for ${cellPopulation} population"))
    # Get the RaggedExperiment for this cellPopulation
    peaksExperiment <- tileResults[[cellPopulation]]
    
    # Get consensus tiles for this cell population for filtering
    consensusTiles <- scMACS:::singlePopulationConsensusTiles(
      peaksExperiment,
      cellPopulation,
      threshold,
      groupColumn,
      join
    )
    
    # consensusTiles is used to filter rows (tiles) from this matrix
    sampleTileIntensityMat <- scMACS:::singlePopulationSampleTileMatrix(
      peaksExperiment,
      consensusTiles,
      NAtoZero,
      log2Intensity
    )
    
    experimentList[[idx]] <- sampleTileIntensityMat
    
    idx <- idx + 1
  }
  
  names(experimentList) <- cellPopulations
  results <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = experimentList,
    colData = sampleData
  )

  return(results)
}
