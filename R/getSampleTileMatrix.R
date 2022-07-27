#' @title \code{create_peak_sampleMatrix}
#'
#' @description \code{create_peak_sampleMatrix} is a function that can transform
#'   a set of sample-specific peak calls into peak X sample matrix containing
#'   lambda1 measurements for each sample.
#'
#'
#' @param tileResults output of callOpenTiles
#' @return sampleTileIntensityMat a sample X peak matrix containing observed
#'   measurements for each sample at each peak.
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
