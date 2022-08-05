#' @title \code{create_peak_sampleMatrix}
#'
#' @description \code{create_peak_sampleMatrix} is a function that can transform
#'   a set of sample-specific peak calls into peak X sample matrix containing
#'   lambda1 measurements for each sample.
#'
#'
#' @param tileResults output of callOpenTiles
#' @return sampleTileIntensityMat a sample X peak matrix containing observed
#'   measurements for each sample at each peak.
#'
#' @details The technical details of the algorithm are found in XX.
#'
#' @references XX
#'
#' @export

getSampleTileMatrix <- function(tileResults,
                                cellPopulations = "ALL",
                                groupColumn = NULL,
                                threshold = 0.2,
                                join = "union",
                                NAtoZero = FALSE,
                                log2Intensity = TRUE,
				numCores = 1,
				verbose = FALSE
                               ) {

  if (!(join %in% c("union", "intersect"))){
    stop("`join` must be either 'union' or 'intersect'")
  }
  
  if(class(tileResults)[1] != 'MultiAssayExperiment'){
    stop('tileResults is not a MultiAssayExperiment')  
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
  
  if (!(all(cellPopulations %in% names(tileResults))) & tolower(cellPopulations) != 'all'){
       stop(paste(
          "All of `cellPopulations` must present in tileResults.",
          "Check `names(tileResults)` for possible cell populations."
        ))
  }
  
  if(length(cellPopulations) == 1){
    if (cellPopulations == 'ALL'){
      subTileResults <- tileResults
      cellPopulations = names(tileResults)
    }
  }else{

    subTileResults <- tileResults[names(tileResults) %in% cellPopulations]

  }
  

  if(verbose){message(str_interp("Extracting consensus tile set within each population:  ${cellPopulations} "))}

  peakList <- mclapply(experiments(subTileResults), function(x){
    
    # Get consensus tiles for this cell population for filtering
    scMACS:::singlePopulationConsensusTiles(
      x,
      threshold,
      groupColumn,
      join
    )

  }, mc.cores = numCores)

  allPeaks <- sort(unique(do.call('c',peakList)))

  if(verbose){message(str_interp("Generating Sample Tile Matrix Across All Populations: ${cellPopulations} "))}

  # consensusTiles is used to  extract rows (tiles) from this matrix
  sampleTileIntensityMatList <- mclapply(experiments(tileResults), function(x){

		scMACS:::singlePopulationSampleTileMatrix(
      				x,
      				allPeaks,
      				NAtoZero,
      				log2Intensity
    		)
   }, mc.cores = numCores)

  peakPresence <- lapply(peakList, function(x) (allPeaks %in% x)) %>% do.call('cbind', .) %>% as.data.frame()

  allPeakGR <- scMACS::StringsToGRanges(sort(allPeaks))
  mcols(allPeakGR) <- peakPresence


   sampleTileIntensityMatList2 <- lapply(sampleTileIntensityMatList, function(x){
					rownames(x) = NULL
					x

   })

  #names(experimentList) <- cellPopulations
  results <- SummarizedExperiment::SummarizedExperiment(
  		sampleTileIntensityMatList2,
    		rowRanges = allPeakGR,
		colData = sampleData
  )

  return(results)
}
