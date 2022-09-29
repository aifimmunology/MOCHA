#' @title \code{getCellPopMatrix}
#'
#' @description \code{getCellPopMatrix} pulls out the SampleTileMatrix of tiles called in one
#'   given cell population.
#'
#' @param SampleTileObj The output from getSampleTileMatrix, a SummarizedExperiment of pseudobulk 
#'   intensities across all tiles & cell types.

getCellPopMatrix <- function(SampleTileObj, cellPopulation, dropSamples = TRUE, NAtoZero = TRUE){
    
    tilesCalled <- GenomicRanges::mcols(SummarizedExperiment::rowRanges(SampleTileObj))[,cellPopulation]
        
    sampleTileMatrix <- SummarizedExperiment::assays(SampleTileObj)[[cellPopulation]][tilesCalled, , drop = FALSE]
		
		#Drop empty samples
		if(dropSamples){sampleTileMatrix <- sampleTileMatrix[, colSums(sampleTileMatrix, na.rm = TRUE) > 0] }
		
		if(NAtoZero){ sampleTileMatrix[is.na(sampleTileMatrix)] <- 0 }

    return(sampleTileMatrix)
}
