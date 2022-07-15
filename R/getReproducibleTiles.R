#' @title \code{get_reproducible_peaks}
#'
#' @description \code{get_reproducible_peaks} is an R helper function, part of
#'   the single-cell peak
#'
#' @export

getReproducibleTiles <- function(tileResults,
                                 cellPopulation,
                                 reproducibility = 0.2) {

  # Get the RaggedExperiment for this cellPopulation
  peaksExperiment <- tileResults[[cellPopulation]]
  nSamples <- length(colnames(peaksExperiment))

  # Extract matrices of samples by peak tileIDs with peak status
  # and set NA to FALSE (no peak here)
  samplePeakMat <- RaggedExperiment::compactAssay(
    peaksExperiment,
    i = "peak"
  )
  samplePeakMat[is.na(samplePeakMat)] <- FALSE

  # Count the number of TRUE peaks in each row and divide by nSamples
  percentTruePeaks <- rowSums(samplePeakMat) / nSamples

  # Keep only peaks with reproducibility above specified threshold
  reproduciblePeaks <- percentTruePeaks[percentTruePeaks >= reproducibility]
  reproduciblePeaks <- names(reproduciblePeaks)

  return(reproduciblePeaks)

}
