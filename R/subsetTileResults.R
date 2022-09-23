#' @title \code{subsetTileResults}
#'
#' @description \code{subsetTileResults} Extracts the peak reproducibility and generates a 
#'					heuristic plots that can be used to determine the reproducibility threshold
#'					used within getSampleTileMatrix. 
#' @param tileObject A MultiAssayExperiment object from callOpenTiles, 
#' @param cellPopulations the cell populations you want to visualize.
#' @param groupColumn Optional parameter, same as in getSampleTileMatrix, which defines whether you
#' 				  want to plot reproducibility within each 
#' @param returnPlotList Instead of one plot with all celltypes/conditions, it returns a list of plots for each cell types 
#' @param returnDFs Instead of a plot, returns a data.frame of the reproducibility across samples. 
#' 						If set to false, then it plots the data.frame instead of returning it. 
#' @param numCores Number of cores to multithread over. 
#'
#' @return SampleTileObj the input data structure with added gene annotations.
#'
#'
#' @export


subsetTileResults <- function(tileObject,
						  subsetBy,
                          groupList,
						  numCores = 1
						  ) {
						  
		
  sampleData <- MultiAssayExperiment::colData(tileObject)
   
  if (!(subsetBy %in% colnames(sampleData)) & tolower(subsetBy) != 'celltypes') {
    stop("Error: subsetBy must either be a column name within colData(tileObjects), or 'celltype'.")
  } 
  
  if (subsetBy %in% colnames(sampleData)) {
    if(all(groupList %in% unique(sampleData[[subsetBy]]))){ 
		stop("Error: groupList includes names not found within the object sample data. Please check groupList.")
	}
	
	
	
	
  } 
  
  if (tolower(subsetBy) == 'celltypes') {
    if(all(groupList %in% names(tileObject))){ 
		stop("Error: groupList includes celltypes not found within tileObject.")
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
							ggplot(alldf[[x]], aes(x = Reproducibility, y = PeakNumber, group = GroupName, color = GroupName)) + 
								geom_point() + ggtitle(names(alldf)[x]) + scale_y_continuous(trans='log2') + 
								ylab('Peak Number') + theme_bw()
							
						)
		
		}else{
		
			allPlots <- lapply(seq_along(alldf), function(x) 
						ggplot(alldf[[x]], aes(x = Reproducibility, y = PeakNumber)) + 
								geom_point() + ggtitle(names(alldf)[x]) + scale_y_continuous(trans='log2') + 
								ylab('Peak Number')  + theme_bw()
							
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
  
	p1 <- ggplot(combinedDF, aes(x = Reproducibility, y = PeakNumber, group = groups, color = groups)) + geom_line() + 
			scale_y_continuous(trans='log10') + ylab('Peak Number')  + theme_bw()
	
	 return
	 (p1)
  }
  
 
   
}