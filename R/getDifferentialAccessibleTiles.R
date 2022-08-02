#' @title \code{get_differential_accessible_Regions}
#'
#' @description \code{get_differential_accessible_Regions} allows you to
#'   determine whether regions of chromatin are differentially accessible
#'   between groups by conducting a test
#'
#'
#'
#' @export
getDifferentialAccessibleTiles <- function(sampleTileMatrix,
                                           tileResults,
                                           fdr_control = 0.2,
                                           numCores = 2) {
  metaFile <- MultiAssayExperiment::colData(tileResults)

  ## Get group labels
  positive_samples <- metaFile$SampleCellType[metaFile$Class == "Positive"]
  group <- ifelse(
    names(sampleTileMatrix)[2:ncol(sampleTileMatrix)] %in% positive_samples, 1, 0
  )

  ## Estimate differential accessibility
  res_pvals <- mclapply(sampleTileMatrix$tileID,
    function(x) {
      estimate_differential_accessibility(sampleTileMatrix, x, group, F)
    },
    mc.cores = nCores
  )

  ## combine results into single objects
  res_pvals <- rbindlist(res_pvals)

  ##
  thresh <- find_alpha_threshold(res_pvals, fdr_control)

  res_pvals$FDR_Controlled <- res_pvals$P_value < thresh

  discoveries <- sum(res_pvals$FDR_Controlled)
  print(paste(
    discoveries, " differential regions found at FDR",
    fdr_control
  ))
  return(res_pvals)
}
