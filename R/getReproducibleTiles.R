#' @title \code{get_reproducible_peaks}
#'
#' @description \code{get_reproducible_peaks} is an R helper function, part of
#'   the single-cell peak
#'
#' @export

getReproducibleTiles <- function(tileResults,
                                 cellPopulation = 'ALL',
                                 threshold = 0.2,
                                 groupColumn = NULL,
                                 join = "union") {
  
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
  
  if(class(tileResults)[1] != 'MultiAssayExperiment'){
    
      stop('tileResults is not a MultiAssayExperiment')  
    
  }else if(cellPopulation == 'ALL'){
    
      cellPopulation = names(tileResults)
    
  }
  
  
  if (!(cellPopulation %in% names(tileResults))){
    stop(paste(
      "`cellPopulation` must present in tileResults.",
      "Check `names(tileResults)` for possible cell populations."
    ))
  }

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
  

  if (is.null(groupColumn)) {

    # Count the number of TRUE peaks in each row and divide by nSamples
    percentTruePeaks <- rowSums(samplePeakMat) / nSamples

    # Keep only peaks with reproducibility above specified threshold
    reproduciblePeaks <- percentTruePeaks[percentTruePeaks >= threshold]
    reproduciblePeaks <- names(reproduciblePeaks)
    
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
      reproduciblePeaks <- percentTruePeaks[percentTruePeaks >= threshold]
      reproduciblePeaks <- names(reproduciblePeaks)
      
      consensusPeaksByGroup <- append(consensusPeaksByGroup, list(reproduciblePeaks))
      
    }
    names(consensusPeaksByGroup) <- groups
    
    if (join=="union"){
      reproduciblePeaks <- Reduce(union, consensusPeaksByGroup)
    } else if (join=="intersect"){
      reproduciblePeaks <- Reduce(intersect, consensusPeaksByGroup)
    }

  }
  return(reproduciblePeaks)
}
