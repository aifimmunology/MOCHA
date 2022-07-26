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
