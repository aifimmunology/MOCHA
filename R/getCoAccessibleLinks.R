#' @title \code{getCoAccessibleLinks}
#'
#' @description \code{getCoAccessibleLinks} takes an input set of regions (tiles) and finds co-accessible neighboring regions within a window. Co-accessibility is defined as the correlation between two region intensity (openness) across samples.
#'
#'
#' @param SampleTileObj The SummarizedExperiment object output from getSampleTileMatrix containing your sample-tile matrices
#' @param cellPopulation A string denoting the cell population of interest, which must be present in SampleTileObj
#' @param regions a GRanges object or vector or strings containing the regions on which to compute co-accessible links. Strings must be in the format "chr:start-end", e.g. "chr4:1300-2222".
#'   Can be the output from getDifferentialAccessibleTiles.
#' @param windowSize the size of the window, in basepairs, around each input region to search for co-accessible links
#' @param numCores Optional, the number of cores to use with multiprocessing. Default is 1.
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @param ZI boolean flag that enables zero-inflated (ZI) Spearman correlations to be used. Default is TRUE. If FALSE, skip zero-inflation and calculate the normal Spearman.
#' 
#' @return TileCorr A data.table correlation matrix
#'
#' @details The technical details of the zero-inflated correlation can be
#'          found here:
#'
#'               Pimentel, Ronald Silva, "Kendall's Tau and Spearman's Rho
#'               for Zero-Inflated Data" (2009). Dissertations.
#'
#'          while the implementation (scHOT R package), can be found here:
#'               http://www.bioconductor.org/packages/release/bioc/html/scHOT.html
#'
#'
#' @export
getCoAccessibleLinks <- function(SampleTileObj, 
                                 cellPopulation = 'All', 
                                 regions, 
                                 windowSize = 1 * 10^6, 
                                 numCores = 1, 
                                 ZI = TRUE,
                                 verbose = FALSE) {

  if (methods::is(regions, "GRanges")) {
    regionDF <- as.data.frame(regions)
  } else if (methods::is(regions,"character")) {
    regionDF <- MOCHA::StringsToGRanges(regions) %>% as.data.frame()
  } else {
    stop('Invalid input type for "region": must be either "GRanges" or a character vector')
  }

  if(cellPopulation == 'All'){
  
	  tileDF <- do.call('cbind', as.list(SummarizedExperiment::assays(SampleTileObj)))
	
  }else if(length(cellPopulation) > 1 & all(cellPopulation %in% names(assays(SampleTileObj)))){
  
 	  tileDF <- do.call('cbind', assays(SampleTileObj)[names(assays(SampleTileObj)) %in% cellPopulation])
	
  }else if(length(cellPopulation) == 1 & all(cellPopulation %in% names(assays(SampleTileObj)))){
  
 	  tileDF <- MOCHA::getCellPopMatrix(SampleTileObj, cellPopulation, NAtoZero = TRUE)
	
  }else{
  
	  error('Cell type not found within SampleTileObj')
  
  }

  tileDF[is.na(tileDF)] = 0

  tileNames <- rownames(tileDF)

  start <- as.numeric(gsub("chr.*\\:|\\-.*", "", rownames(tileDF)))
  end <- as.numeric(gsub("chr.*\\:|.*\\-", "", rownames(tileDF)))
  chr <- gsub("\\:.*", "", rownames(tileDF))


  #Find all combinations to test

  allCombinations <- pbapply::pblapply(1:dim(regionDF)[1], function(y){

    keyTile <- which(start == regionDF$start[y] &
      end == regionDF$end[y]  &
      chr == regionDF$seqnames[y])

    windowIndexBool <- which(start > regionDF$start[y] - windowSize / 2 &
      end < regionDF$end[y] + windowSize / 2 &
      chr == regionDF$seqnames[y])

    windowIndexBool == windowIndexBool[windowIndexBool != keyTile]

    regionOfInterest <- tileNames[keyTile]
    allOtherRegions <- tileNames[windowIndexBool]

    # Var1 will always be our region of interest
    keyNeighborPairs <- data.frame(
      "Key" = regionOfInterest,
      "Neighbor" = allOtherRegions
    )

    keyNeighborPairs

  }, cl = numCores) %>% do.call('rbind', .)   %>% dplyr::distinct()

  chrNum <- paste(regionDF$chr,":",sep='')

  zi_spear_mat <- NULL

  return(allCombinations)
  
  for(i in chrNum){

    subTileDF <- tileDF[grepl(i, rownames(tileDF)),]
    subCombinations <- allCombinations[grepl(i, allCombinations$Key),]

    # General case for >1 pair
    zero_inflated_spearman <- unlist(pbapply::pblapply(1:nrow(subCombinations),
      function(x) {
        weightedZISpearman(
          x = subTileDF[subCombinations[x, "Key"], ],
          y = subTileDF[subCombinations[x, "Neighbor"], ],
          verbose = verbose,
          ZI = ZI
        )
      },
      cl = numCores
    ))

    # Create zero-inflated correlation matrix from correlation values
    zi_spear_mat_tmp <- data.table::data.table(
      Correlation = zero_inflated_spearman,
      Tile1 = subCombinations[, "Key"],
      Tile2 = subCombinations[, "Neighbor"]
    )

    zi_spear_mat <- rbind(zi_spear_mat,zi_spear_mat_tmp)

  }   
  return(zi_spear_mat)
}
