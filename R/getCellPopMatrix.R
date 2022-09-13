#' @title \code{getCellPopMatrix}
#'
#' @description \code{getCellPopMatrix} pulls out the SampleTileMatrix of tiles called in one
#'   given cell population.
#'
#' @param SampleTileObj The output from getSampleTileMatrix, a SummarizedExperiment of pseudobulk
#'   intensities across all tiles & cell types.
#' @param @cellPopulation - The cell population you want to pull out.
#'
#' @return sampleTileMatrix a matrix of samples by called tiles for a given cell population.
#'
#' @export
getCellPopMatrix <- function(SampleTileObj, cellPopulation) {
  tilesCalled <- GenomicRanges::mcols(SummarizedExperiment::rowRanges(SampleTileObj))[, cellPopulation]

  sampleTileMatrix <- SummarizedExperiment::assays(SampleTileObj)[[cellPopulation]][tilesCalled, , drop = FALSE]

  return(sampleTileMatrix)
}
