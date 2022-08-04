#' @title \code{singlePopulationConsensusTiles}
#'
#' @description \code{singlePopulationConsensusTiles} is an R helper function, part of
#'   the single-cell peak
#'
#' @keywords internal

singlePopulationConsensusTiles <- function(peaksExperiment,
                                 threshold,
                                 groupColumn = NULL,
                                 join = 'union') {
  
  if (!(join %in% c("union", "intersect"))){
    stop("`join` must be either 'union' or 'intersect'")
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

  # Extract matrices of samples by peak tileIDs with peak status
  # and set NA to FALSE (no peak here)
  samplePeakMat <- RaggedExperiment::compactAssay(
    peaksExperiment,
    i = "peak"
  )
  samplePeakMat[is.na(samplePeakMat)] <- FALSE
  nSamples <- length(colnames(peaksExperiment))
  
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
    
    if (join=="union"){
      consensusPeaks <- Reduce(union, consensusPeaksByGroup)
    } else if (join=="intersect"){
      consensusPeaks <- Reduce(intersect, consensusPeaksByGroup)
    }

  }
  return(consensusPeaks)
}
