#' @title \code{co_accessibility}
#'
#' @description \code{co_accessibility} allows you to determine whether 2 peaks
#'              are co-accessible using a zero-inflated spearman correlation.
#'
#'
#' @param mat sample-peak matrix with regions to analyze
#' @param numCores integer to determine # of paralell cores
#'
#' @return a 3-column data.frame containing
#'         - Correlation= Zero-inflated Spearman Correlation
#'         - Peak1= location of co-accessible region 1
#'         - Peak2= location of co-accessible region 2
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
#' @references XX
#' @examples
#' Generate
#' mat1 <- matrix(pmax(0, rnorm(1000)), ncol = 100)
#' row.names(mat1) <- paste("A", 1:10, sep = "_")
#' ziSpear_mat <- co_accessibility(mat1, numCores = 5)
#' head(ziSpear_mat)
#'
#' @internal
#' @noRd

co_accessibility <- function(subMat, filterPairs = NULL, index = NULL, numCores = 40, verbose = FALSE) {

  ## Prep data.table of peaks for correlations.
  ## mat1 <- as.data.table(mcols(GR))
  ## rownames(mat1) <- GRangesToString(GR)

  # mat1 <- subMat[,-c('tileID','chr', 'start', 'end')]
  # rownames(mat1) = subMat$tileID
  pairwise_combos <- RcppAlgos::comboGrid(rownames(subMat), rownames(subMat), repetition = FALSE)

  ## Pull out only combinations that involve the peak of interest
  if (!is.null(index) & class(index) == "integer") {
    pairwise_combos <- pairwise_combos[pairwise_combos[, "Var1"] %in% rownames(subMat)[index] |
      pairwise_combos[, "Var2"] %in% rownames(subMat)[index], ]
  } else if (verbose) {
    print("No valid index given. Testing all peaks within range.")
  }


  ## Filter out any pairs that have already been tested, and see if anything remains.
  if (!is.null(filterPairs) & length(pairwise_combos) > 2) {

    # If there's more than one pair of peaks, filter out all the previously tested ones.
    pairwise_combos <- pairwise_combos[!(pairwise_combos[, "Var1"] %in% filterPairs$Peak1 &
      pairwise_combos[, "Var2"] %in% filterPairs$Peak2), ]
  } else if (!is.null(filterPairs) & length(pairwise_combos) > 0) {

    # If only one pair of peaks is left, and that pair has already been tested, then return NULL
    if (pairwise_combos[1] %in% filterPairs$Peak1 & pairwise_combos[2] %in% filterPairs$Peak2) {
      return(NULL)
    }
  }


  ### Loop through all pairwise
  ### combinations of peaks


  ## if only one pair of peaks left, then nrows will hit an error. length works whether there is one pair or many left.
  # If no peaks
  if (length(pairwise_combos) == 0 & verbose) {
    print("Warning: No peaks found in neighborhood of region of interest")
    return(NULL)
  } else if (length(pairwise_combos) == 0) {
    return(NULL)
  } else if (length(pairwise_combos) == 2) {

    # If only one pair of peaks to test, then it's no longer a data.frame, but a vector.
    zero_inflated_spearman <- scMACS:::weightedZISpearman(
      x = subMat[pairwise_combos[1], ],
      y = subMat[pairwise_combos[2], ]
    )

    zi_spear_mat <- data.table(
      Correlation = zero_inflated_spearman,
      Peak1 = pairwise_combos[1],
      Peak2 = pairwise_combos[2]
    )
  } else {
    zero_inflated_spearman <- unlist(mclapply(1:nrow(pairwise_combos),
      function(x) {
        scMACS:::weightedZISpearman(
          x = subMat[pairwise_combos[x, "Var1"], ],
          y = subMat[pairwise_combos[x, "Var2"], ]
        )
      },
      mc.cores = numCores
    ))


    # return(zero_inflated_spearman)
    ### Create zero-inflated correlation matrix
    ### from correlation values,
    zi_spear_mat <- data.table(
      Correlation = zero_inflated_spearman,
      Peak1 = pairwise_combos[, "Var1"],
      Peak2 = pairwise_combos[, "Var2"]
    )
  }

  return(zi_spear_mat)
}

