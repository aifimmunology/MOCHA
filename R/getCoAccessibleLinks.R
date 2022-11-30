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

  start <- as.numeric(gsub("chr.*\\:|\\-.*", "", rownames(tileDF)))
  end <- as.numeric(gsub("chr.*\\:|.*\\-", "", rownames(tileDF)))
  chr <- gsub("\\:.*", "", rownames(tileDF))

  # Initialize correlation datatable
  TileCorr <- NULL
  pb <- utils::txtProgressBar(min = 0, max = length(regions), initial = 0, style = 3)

  for (i in 1:length(regions)) {
    utils::setTxtProgressBar(pb, i)

    # Find all neighboring tiles in the window
    windowIndexBool <- which(start > regionDF$start[i] - windowSize / 2 &
      end < regionDF$end[i] + windowSize / 2 &
      chr == regionDF$seqnames[i])

    if (length(windowIndexBool) > 1) {

      # The region of interest should overlap with the tile at least partially.
      keyTile <- which(
        windowIndexBool == which(
          ((start <= regionDF$start[i] & end >= regionDF$start[i]) |
            (start >= regionDF$start[i] & end <= regionDF$end[i]) |
            (start <= regionDF$end[i] & end >= regionDF$end[i])) &
            chr %in% regionDF$seqnames[i]
        )
      )
      nextCorr <- MOCHA:::co_accessibility(tileDF[windowIndexBool, , drop = FALSE],
        filterPairs = TileCorr, numCores = numCores,
        index = keyTile, verbose = verbose, ZI = ZI
      )
      # For first iteration, TileCorr is NULL so it will be ignored
      TileCorr <- rbind(TileCorr, nextCorr)
    } else if (verbose) {
      warning(
        "No neighboring tiles found for given region: ",
        stringr::str_interp(
          "${regionDF$seqnames[i]}:${regionDF$start[i]}-${regionDF$end[i]}"
        ),
        stringr::str_interp(" with windowSize ${windowSize} basepairs")
      )
    }
  }

  close(pb)
  TileCorr
}
