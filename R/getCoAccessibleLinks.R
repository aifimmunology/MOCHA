#' @title \code{getCoAccessibleLinks}
#'
#' @description \code{getCoAccessibleLinks} does stuff 
#'
#'
#' @param SampleTileObj The SummarizedExperiment object output from getSampleTileMatrix
#' @param regions a GRanges containing the regions to compute co-accessible links for. 
#'   Can be the output from getDifferentialAccessibleTiles.
#' @param cellPopulation A string denoting the cell population of interest
#' @param numCores Optional, the number of cores to use with multiprocessing. Default is 1.
#'
#' @return 
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
getCoAccessibleLinks <- function(peakDT, regions, windowSize = 1 * 10^6, numCores = 1, verbose = FALSE) {
  peakDF <- as.data.frame(peakDT[, -"tileID"])
  rownames(peakDF) <- peakDT$tileID
  
  start <- as.numeric(gsub("chr.*\\:|\\-.*", "", peakDT$tileID))
  end <- as.numeric(gsub("chr.*\\:|.*\\-", "", peakDT$tileID))
  chr <- gsub("\\:.*", "", peakDT$tileID)
  
  regionDF <- as.data.frame(regions)
  
  # Initialize Correlation Datatable.
  PeakCorr <- NULL
  pb <- txtProgressBar(min = 0, max = length(regions), initial = 0)
  
  for (i in 1:length(regions)) {
    # print(i)
    setTxtProgressBar(pb, i)
    
    ## Find all peaks in the next window
    window_tmp <- which(start > regionDF$start[i] - windowSize / 2 &
                          end < regionDF$end[i] + windowSize / 2 &
                          chr == regionDF$seqnames[i])
    
    if (length(window_tmp) > 1) {
      
      # The region of interest should overlap with the peak at least partially.
      keyPeak <- which(
        window_tmp == which(
          ((start <= regionDF$start[i] & end >= regionDF$start[i]) |
             (start >= regionDF$start[i] & end <= regionDF$end[i]) |
             (start <= regionDF$end[i] & end >= regionDF$end[i])) &
            chr %in% regionDF$seqnames[i]
        )
      )
      
      nextCorr <- scMACS:::co_accessibility(peakDF[window_tmp, ],
                                   filterPairs = PeakCorr, numCores = numCores,
                                   index = keyPeak, verbose = verbose
      )
      PeakCorr <- rbind(PeakCorr, nextCorr)
    } else if (verbose) {
      print("Warning: No correlations found for given region.")
    }
  }
  
  close(pb)
  gc()
  return(PeakCorr)
}