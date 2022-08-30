#' @title \code{FilterCoAccessibleLinks}
#'
#' @description \code{FilterCoAccessibleLinks} does stuff 
#'
#'
#' @param PeakCorr 
#' @param threshold 
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
FilterCoAccessibleLinks <- function(PeakCorr, threshold = 0.5) {
  if (!any(abs(PeakCorr$Correlation) > threshold)) {
    stop("Error: There are no values above the threshold.")
  }
  
  PeakCorr1 <- PeakCorr[abs(PeakCorr$Correlation) > threshold, ]
  start1 <- as.numeric(gsub("chr.*\\:|\\-.*", "", PeakCorr1$Peak1))
  end1 <- as.numeric(gsub("chr.*\\:|.*\\-", "", PeakCorr1$Peak1))
  
  start2 <- as.numeric(gsub("chr.*\\:|\\-.*", "", PeakCorr1$Peak2))
  end2 <- as.numeric(gsub("chr.*\\:|.*\\-", "", PeakCorr1$Peak2))
  
  PeakCorr1$chr <- gsub("\\:.*", "", PeakCorr1$Peak2)
  PeakCorr1$start <- apply(data.table(start1, start2), 1, min)
  PeakCorr1$end <- apply(data.table(end1, end2), 1, max)
  
  return(PeakCorr1)
}


#############################################