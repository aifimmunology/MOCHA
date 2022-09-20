#' @title \code{plotReproducibility}
#'
#' @description \code{plotReproducibility} Extracts the peak reproducibility and generates a 
#'					heuristic plots that can be used to determine the reproducibility threshold
#'					used within getSampleTileMatrix. 
#' @param tileObject A MultiAssayExperiment object from callOpenTiles, 
#' @param cellPopulations the cell populations you want to visualize.
#' @param groupColumn Optional parameter, same as in getSampleTileMatrix, which defines whether you
#' 				  want to plot reproducibility within each 
#' defining the promoter region. The format is (upstream, downstream). 
#' Default is (2000, 100).
#'
#' @return SampleTileObj the input data structure with added gene annotations.
#'
#'
#' @export
plotReproducibility <- function(tileObject,
						  cellPopulations  = 'All',
                          groupColumn = NULL,
						  returnPlotList = FALSE,
						  returnDFs = FALSE,
						  numCores = 1
						  ) {
						  
						  
   if (length(cellPopulations) == 1 & all(tolower(cellPopulations) == "all")) {
    subTileResults <- tileObject
    cellPopulations <- names(tileObject)
  } else {
    if (all(cellPopulations %in% names(tileObject))) {
      subTileResults <- tileObject[names(tileObject) %in% cellPopulations]
    } else {
      stop(paste(
        "All of `cellPopulations` must present in tileResults.",
        "Check `names(tileResults)` for possible cell populations."
      ))
    }
  }
  
  sampleData = MultiAssayExperiment::colData(tileObject)
  
  alldf <- mclapply(names(subTileResults), function(x){
					cellTypeDF(peaksExperiment = subTileResults[[x]], sampleData = sampleData, 
									groupColumn = groupColumn, returnPlotList = returnPlotList)
					}, mc.cores = numCores)
					
  names(alldf) <- names(subTileResults)
  
  if(returnDFs){ return(alldf)}
  
  if(returnPlotList){
  
		if(!is.null(groupColumn)){
		
			allPlots <- lapply(seq_along(alldf), function(x) 
							ggplot(alldf[[x]], aes(x = Reproducibility, y = log(PeakNumber), group = GroupName, color = GroupName)) + 
								geom_point() + ggtitle(names(alldf)[x])
							
						)
		
		}else{
		
			allPlots <- lapply(seq_along(alldf), function(x) 
						ggplot(alldf[[x]], aes(x = Reproducibility, y = log(PeakNumber))) + 
								geom_point() + ggtitle(names(alldf)[x])
							
						)
		}
						
		return(allPlots)
  
  }else{
  
	combinedDF <- do.call('rbind', alldf)
	combinedDF$CellPop <- gsub("\\.","", gsub("[0-9]{1,3}","",rownames(combinedDF)))
	
	if(!is.null(groupColumn)){
		combinedDF$groups = paste(combinedDF$CellPop, combinedDF$GroupName, sep = "_")
	}else{
	
		combinedDF$groups = combinedDF$CellPop
	
	}
  
	p1 <- ggplot(combinedDF, aes(x = Reproducibility, y = log(PeakNumber), group = groups, color = groups)) + geom_line()
	return(p1)
 
  }
   
}


cellTypeDF <- function(peaksExperiment, sampleData,groupColumn = NULL, returnPlotList = FALSE){

  samplePeakMat <- RaggedExperiment::compactAssay(
							peaksExperiment,521
							i = "peak"
						)
  
  samplePeakMat[is.na(samplePeakMat)] <- FALSE
  nSamples <- length(colnames(peaksExperiment))
  
  if (is.null(groupColumn)) {

		# Count the number of TRUE peaks in each row and divide by nSamples
		peakTable <- table(rowSums(samplePeakMat) / nSamples)
		TruePeaks <- data.frame('Reproducibility' = as.numeric(as.character(names(peakTable))), 'PeakNumber' = unlist(as.list(peakTable)))
		rownames(TruePeaks) = NULL
  }else{
    
    # Get consensus peaks for groupings of samples 
    groups <- unique(sampleData[[groupColumn]])
    TruePeaks <- data.frame()
    for (group in groups){
      # Filter sample-peak matrix to samples in this group
      samplesInGroupDF <- sampleData[sampleData[[groupColumn]]==group,]
      samplesInGroup <- rownames(samplesInGroupDF)
      groupSamplePeakMat <- samplePeakMat[,samplesInGroup]
	  
	  peakTable <- table(rowSums(groupSamplePeakMat) / length(samplesInGroup))
      tmp <- data.frame('Reproducibility' = as.numeric(as.character(names(peakTable))), 'PeakNumber' = unlist(as.list(peakTable)))
	  tmp$GroupName = rep(group, length(peakTable))
	  rownames(tmp) = NULL
	  TruePeaks <- rbind(TruePeaks, tmp)
      
    }
	
  }
  
  return(TruePeaks)


}
