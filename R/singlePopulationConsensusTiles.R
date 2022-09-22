#' @title \code{singlePopulationConsensusTiles}
#'
#' @description \code{singlePopulationConsensusTiles} is an R helper function, part of
#'   the single-cell peak
#'
#' @keywords internal

singlePopulationConsensusTiles <- function(peaksExperiment,
                                           sampleData,
                                           threshold,
                                           groupColumn = NULL) {
    
  
  # Extract matrices of samples by peak tileIDs with peak status
  # and set NA to FALSE (no peak here)
  samplePeakMat <- RaggedExperiment::compactAssay(
    peaksExperiment,
    i = "peak"
  )
  
  # Identify any samples that had less than 5 cells, and therefore no peak calls
  emptySamples <- apply(samplePeakMat, 2, function(x) all(is.na(x) | !x))
  
  #Now we'll filter out those samples, and replace all NAs with zeros. 
  samplePeakMat <- samplePeakMat[,!emptySamples]
  samplePeakMat[is.na(samplePeakMat)] <- FALSE
  
  #The number of samples for reproducibility will be all samples with peak calls
  nSamples <- sum(!emptySamples)
  
  if (is.null(groupColumn)) {

    # Count the number of TRUE peaks in each row and divide by nSamples
    percentTruePeaks <- rowSums(samplePeakMat) / nSamples

    # Keep only peaks with reproducibility above specified threshold
    consensusPeaks <- percentTruePeaks[percentTruePeaks >= threshold]
    consensusPeaks <- names(consensusPeaks)
    
  } else {
    
    # Get consensus peaks for groupings of samples
    groups <- unique(sampleData[[groupColumn]])
    consensusPeaksByGroup <- list()
    for (group in groups){
      # Filter sample-peak matrix to samples in this group
      samplesInGroupDF <- sampleData[sampleData[[groupColumn]]==group,]
      samplesInGroup <- rownames(samplesInGroupDF)
      groupSamplePeakMat <- samplePeakMat[,samplesInGroup]

      # Keep only peaks with reproducibility above specified threshold
      percentTruePeaks <- rowSums(groupSamplePeakMat) / nSamples
      consensusPeaks <- percentTruePeaks[percentTruePeaks >= threshold]
      consensusPeaks <- names(consensusPeaks)
      
      consensusPeaksByGroup <- append(consensusPeaksByGroup, list(consensusPeaks))
      
    }
    names(consensusPeaksByGroup) <- groups
	
    consensusPeaks <- Reduce(union, consensusPeaksByGroup)
	
  }
  return(consensusPeaks)
}
