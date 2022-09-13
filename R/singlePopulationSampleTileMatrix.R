#' @title \code{singlePopulationSampleTileMatrix}
#'
#' @description \code{singlePopulationSampleTileMatrix} is a function that can transform
#'   a set of tile intensities into peak X sample matrix for a custom set of tiles
#'
#'
#' @param peaksExperiment peakset raggedExperiment, one celltype from the output of callOpenTiles
#' @param consensusTiles a vector containing the tileIDs to subset the
#'   sample-tile matrix
#' @return sampleTileIntensityMat a sample X peak matrix containing observed
#'   measurements for each sample at each peak.
#'
#' @details The technical details of the algorithm are found in XX.
#'
#' @references XX
#'
#' @keywords internal

singlePopulationSampleTileMatrix <- function(peaksExperiment,
                                             consensusTiles,
                                             NAtoZero = FALSE,
                                             log2Intensity = TRUE) {

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
    sampleTileIntensityMat <- log2(sampleTileIntensityMat + 1)
  }

  # Filter to just peaks in the given consensusTiles.
  prevDims <- dim(sampleTileIntensityMat)
  sampleTileIntensityMat <- sampleTileIntensityMat[rownames(sampleTileIntensityMat) %in% consensusTiles, , drop = FALSE]

  # If the tiles you are pulling don't exist in the peakExperiment, add a matrix of NAs or Zeros to represent those tiles
  if (any(!consensusTiles %in% rownames(sampleTileIntensityMat))) {
    filler <- ifelse(NAtoZero, 0, NA)
    rowList <- consensusTiles[!consensusTiles %in% rownames(sampleTileIntensityMat), drop = FALSE]
    emptyMat <- matrix(filler,
      nrow = length(rowList), ncol = ncol(sampleTileIntensityMat),
      dimnames = list(rowList, colnames(sampleTileIntensityMat))
    )
    sampleTileIntensityMat <- rbind(sampleTileIntensityMat, emptyMat)
  }
  currDims <- dim(sampleTileIntensityMat)
  message(
    "\tDimensions of sample-tile matrix\n",
    "\tbefore and after consensus tiles:\n",
    stringr::str_interp("\t${prevDims} before to ${currDims} after")
  )
  sampleTileIntensityMat <- sampleTileIntensityMat[
    order(rownames(sampleTileIntensityMat)),
    order(colnames(sampleTileIntensityMat)),
    drop = FALSE
  ]

  return(sampleTileIntensityMat)
}
