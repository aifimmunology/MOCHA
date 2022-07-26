#' @title \code{create_peak_sampleMatrix}
#'
#' @description \code{create_peak_sampleMatrix} is a function that can transform
#'   a set of tile intensities into peak X sample matrix for a custom set of tiles
#'
#'
#' @param tileResults custom input tile set
#' @param reproducibleTiles a vector containing the tileIDs to subset the
#'   sample-tile matrix
#' @return sampleTileIntensityMat a sample X peak matrix containing observed
#'   measurements for each sample at each peak.
#'
#' @details The technical details of the algorithm are found in XX.
#'
#' @references XX
#'
#' @export

getManualSampleTileMatrix <- function(tileResults,
                                cellPopulation,
                                reproducibleTiles,
                                NAtoZero = FALSE,
                                log2Intensity = TRUE
                               ) {

  # Get the RaggedExperiment for this cellPopulation
  peaksExperiment <- tileResults[[cellPopulation]]

  # Extract matrices of samples by peak tileIDs with TotalIntensity
  sampleTileIntensityMat <- RaggedExperiment::compactAssay(
    peaksExperiment,
    i = "TotalIntensity"
  )
  
  if (NAtoZero) {
    # Replace NAs with zeroes for zero-inflated hypothesis testing
    sampleTileIntensityMat[is.na(sampleTileIntensityMat)] <- 0
  }
  
  if (log2Intensity) {
    sampleTileIntensityMat <- log2(sampleTileIntensityMat+1)
  }

  # Filter to just peaks in the given reproducibleTiles
  sampleTileIntensityMat <- sampleTileIntensityMat[reproducibleTiles, ]

  return(sampleTileIntensityMat)

}

#' @title \code{create_peak_sampleMatrix}
#'
#' @description \code{create_peak_sampleMatrix} is a function that can transform
#'   a set of sample-specific peak calls into peak X sample matrix containing
#'   lambda1 measurements for each sample.
#'
#'
#' @param tileResults output of callOpenTiles
#' @param reproducibleTiles a vector containing the tileIDs to subset the
#'   sample-tile matrix
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
                                groupColumn,
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
  
  if(cellPopulations == 'ALL'){
      cellPopulations = names(tileResults)
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
  
  experimentList <- list()
  idx <- 1
  for (cellPop in cellPopulations){
    
    message(str_interp("Extracting sampleTileMatrix for ${cellPop} population"))
    # Get the RaggedExperiment for this cellPopulation
    peaksExperiment <- tileResults[[cellPopulation]]

    # Extract matrices of samples by peak tileIDs with TotalIntensity
    sampleTileIntensityMat <- RaggedExperiment::compactAssay(
      peaksExperiment,
      i = "TotalIntensity"
    )
    if (NAtoZero) {
      # Replace NAs with zeroes for zero-inflated hypothesis testing
      sampleTileIntensityMat[is.na(sampleTileIntensityMat)] <- 0
    }
    if (log2Intensity) {
      sampleTileIntensityMat <- log2(sampleTileIntensityMat+1)
    }

    # Extract matrices of samples by peak tileIDs with peak status
    # and set NA to FALSE (no peak here)
    samplePeakMat <- RaggedExperiment::compactAssay(
      peaksExperiment,
      i = "peak"
    )
    samplePeakMat[is.na(samplePeakMat)] <- FALSE
    
    # Get reproducible tiles for this cell population for filtering
    reproducibleTiles <- scMACS:::getReproducibleTiles(
      peaksExperiment,
      cellPopulation,
      threshold,
      groupColumn,
      join
    ) 

    # Filter to just peaks in the given reproducibleTiles
    sampleTileIntensityMat <- sampleTileIntensityMat[reproducibleTiles, ]
    # And add it to the experimentList for this cell population
    experimentList[[idx]] <- sampleTileIntensityMat
    
    idx <- idx + 1
  }
  
  browser()

  names(experimentList) <- cellPopulations
  results <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = experimentList,
    colData = sampleData
  )

  return(results)
}
