#' @title Conduct a differential test between open regions of two sample groups
#'
#' @description \code{getDifferentialAccessibleTiles} allows you to
#'   determine whether regions of chromatin are differentially accessible
#'   between groups by conducting a test
#'
#' @param SampleTileObj The SummarizedExperiment object output from
#'  getSampleTileMatrix
#' @param cellPopulation A string denoting the cell population of interest
#' @param groupColumn The column containing sample group labels
#' @param foreground The foreground group of samples for differential comparison
#' @param background The background group of samples for differential comparison
#' @param signalThreshold Optional value. 
#'  This is Minimum median intensity required to keep tiles for
#'  differential testing to increase statistical power in small sample cohorts.
#'  Default is NULL, at which point the optimal threshold will be found.
#' @param minZeroDiff Minimum difference in average dropout rates across groups
#'  require to keep tiles for differential testing. Default is 0.5 (50\%).
#' @param fdrToDisplay False-discovery rate used only for standard
#'  output messaging. Default is 0.2.
#' @param qValueMethod String describing qvalue method. Can be 'standard', or 'experimental'. 
#'           See methods in the MOCHA manuscript for description of the experimental. 
#'           Otherwise, 'standard' applies standard q value. 
#' @param outputGRanges Outputs a GRanges if TRUE and a data.frame if
#'  FALSE. Default is TRUE.
#' @param numCores The number of cores to use with multiprocessing.
#'  Default is 1.
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#'
#' @return full_results The differential accessibility results as a GRanges or
#'   matrix data.frame depending on the flag `outputGRanges`.
#'
#' @examples
#' \dontrun{
#' cellPopulation <- "MAIT"
#' foreground <- "Positive"
#' background <- "Negative"
#' # Standard output will display the number of tiles found below a false-discovery rate threshold.
#' # This parameter does not filter results and only affects the aforementioned message.
#' fdrToDisplay <- 0.2
#' # Choose to output a GRanges or data.frame.
#' # Default is TRUE
#' outputGRanges <- TRUE
#' # SampleTileMatrices is the output of MOCHA::getSampleTileMatrix
#' differentials <- MOCHA::getDifferentialAccessibleTiles(
#'   SampleTileObj = SampleTileMatrices,
#'   cellPopulation = cellPopulation,
#'   groupColumn = groupColumn,
#'   foreground = foreground,
#'   background = background,
#'   fdrToDisplay = fdrToDisplay,
#'   outputGRanges = outputGRanges,
#'   numCores = numCores
#' )
#' }
#' @export
#' @keywords downstream

getDifferentialAccessibleTiles <- function(SampleTileObj,
                                           cellPopulation,
                                           groupColumn,
                                           foreground,
                                           background,
                                           minZeroDiff = 0.5,
                                           qValueMethod = 'experimental',
                                           signalThreshold = NULL,
                                           qValueThreshold = 0.2,
                                           outputGRanges = TRUE,
                                           numCores = 1,
                                           verbose = FALSE) {
  if (!all(cellPopulation %in% names(SummarizedExperiment::assays(SampleTileObj)))) {
    stop("cellPopulation was not found within SampleTileObj. Check available cell populations with `colData(SampleTileObj)`.")
  }

  if (!is.null(signalThreshold)) {
    if (verbose & signalThreshold < 1) {
      warning("Setting the signalThreshold too low will reduce statistical power. You may inspect the distribution of intensities to set a threshold that will remove highly sparse regions.")
    }
  }

  if(!qValueMethod %in% c('standard', 'experimental')){
  
      stop("qValueMethod must either be set to 'standard' or 'experimental'")     
  }
    
  metaFile <- SummarizedExperiment::colData(SampleTileObj)

  if (!(groupColumn %in% colnames(metaFile))) {
    stop(stringr::str_interp("Provided groupCol '{groupColumn}' not found in the provided SampleTileObj"))
  }
  if (!(foreground %in% metaFile[[groupColumn]])) {
    stop(stringr::str_interp("Provided foreground value is not present in the column {groupColumn} in the provided SampleTileObj"))
  }
  if (!(background %in% metaFile[[groupColumn]])) {
    stop(stringr::str_interp("Provided background value is not present in the column {groupColumn} in the provided SampleTileObj"))
  }


  # Get group labels
  foreground_samples <- metaFile[metaFile[, groupColumn] == foreground, "Sample"]
  background_samples <- metaFile[metaFile[, groupColumn] == background, "Sample"]


  #Run a for-loop over all cell populations
  DAT_list = list()
  for(cellPop in cellPopulation){

      # This will only include called tiles
      sampleTileMatrix <- MOCHA::getCellPopMatrix(SampleTileObj, cellPop, NAtoZero = FALSE)

      # Enforce that the samples included are in foreground and background groups -
      # this can onl be an A vs B comparison, i.e. this ignores other groups in groupCol

      sampleTileMatrix <- sampleTileMatrix[, colnames(sampleTileMatrix) %in%
                                           c(foreground_samples, background_samples), drop = FALSE]


      group <- as.numeric(colnames(sampleTileMatrix) %in% foreground_samples)

      #############################################################################
      # Prioritize high-signal tiles

      # Log2 transform the matrix
      # (input must not be log2 transformed prior to this)
      sampleTileMatrix <- log2(sampleTileMatrix + 1)

      medians_a <- matrixStats::rowMedians(sampleTileMatrix[, which(group == 1), drop = FALSE], 
                                           na.rm = TRUE)
      medians_b <- matrixStats::rowMedians(sampleTileMatrix[, which(group == 0), drop = FALSE], 
                                           na.rm = TRUE)
      #medians_a[is.na(medians_a)] = 0
      #medians_b[is.na(medians_b)] = 0
      # Set NAs to zero
      sampleTileMatrix[is.na(sampleTileMatrix)] <- 0

      zero_A <- rowMeans(sampleTileMatrix[, which(group == 1), drop = FALSE] == 0)
      zero_B <- rowMeans(sampleTileMatrix[, which(group == 0), drop = FALSE] == 0)

      diff0s <- abs(zero_A - zero_B)

      if(is.null(signalThreshold)){
        log2FC_filter = quantile(c(medians_a,medians_b), na.rm = TRUE, probs = 0.5)
      }else{
        log2FC_filter <- signalThreshold
      }
      idx <- which(medians_a > log2FC_filter | medians_b > log2FC_filter | diff0s >= minZeroDiff)    

      ############################################################################
      # Estimate differential accessibility

      ## Let's create a matrix to iterate over. 
      cl <- parallel::makeCluster(numCores)
      res_pvals <- pbapply::pbapply(cl = cl, sampleTileMatrix[idx,], 
                                MARGIN = 1, estimate_differential_accessibility,
                                    group = group)
      parallel::stopCluster(cl)
      res_pvals <- do.call(rbind, res_pvals)               
      ##Add back in other untested tiles for visibility.
      removedTiles = which(!rownames(sampleTileMatrix) %in% names(idx))
      filteredTiles = data.frame( P_value = rep(NA, length(removedTiles)),
              TestStatistic =  rep(NA,length(removedTiles)),
              Log2FC_C = rep(NA, length(removedTiles)),
              MeanDiff = rep(NA, length(removedTiles)),
              Case_mu =  rep(NA, length(removedTiles)),
              Case_rho =  rep(NA, length(removedTiles)),
              Control_mu =  rep(NA,length(removedTiles)),
              Control_rho =  rep(NA, length(removedTiles))
            )
      rownames(filteredTiles) = rownames(sampleTileMatrix)[removedTiles]

      ## Bind both together and reorder to match original order.
      res_pvals = rbind(res_pvals, filteredTiles)
      res_pvals = res_pvals[rownames(sampleTileMatrix),] 
      res_pvals = cbind(data.frame(Tile = rownames(res_pvals)), 
                        res_pvals)


      #############################################################################
      # Apply FDR on filtered regions

      filtered_res <- res_pvals[idx, ]

      if (!all(is.na(filtered_res$P_value[filtered_res$P_value <= 0.95]))) {

        if(qValueMethod == 'experimental'){
            pi0_reduced <- qvalue::pi0est(filtered_res$P_value[filtered_res$P_value <= 0.95],
              pi0.method = "bootstrap",
              lambda = seq(0, 0.6, .05)
            )

            filtered_res$FDR <- qvalue::qvalue(filtered_res$P_value, pi0 = pi0_reduced$pi0)$qvalues
        }else{

            filtered_res$FDR <- qvalue::qvalue(filtered_res$P_value)$qvalues
        }
      } else {
        filtered_res$FDR <- NA # TODO Handle appropriately
      }
      #############################################################################

      # Join with original results
      full_results <- dplyr::left_join(res_pvals, filtered_res[, c("Tile", "FDR")], by = "Tile")

      # Set results to NA where FDR is NA
      na.idx <- which(is.na(full_results$FDR))
      full_results$P_value[na.idx] <- NA
      full_results$TestStatistic[na.idx] <- NA

      sampleTileMatrix[is.na(sampleTileMatrix)] <- 0
      meansA <- rowMeans(sampleTileMatrix[, group == 1, drop = FALSE], na.rm = T)
      meansB <- rowMeans(sampleTileMatrix[, group == 0, drop = FALSE])
      full_results$MeanDiff <- meansA - meansB
      full_results$MeanDiff[is.na(full_results$FDR)] <- NA
      full_results$CellPopulation <- rep(cellPop, length(meansA))
      full_results$Foreground <- rep(foreground, length(meansA))
      full_results$Background <- rep(background, length(meansA))


      colnames(full_results) <- c(
        "Tile", "P_value", "Test_Statistic", "Log2FC_C",
        "MeanDiff", "Avg_Intensity_Case", "Pct0_Case",
        "Avg_Intensity_Control", "Pct0_Control", "FDR",
        "CellPopulation", "Foreground", "Background"
      )

      full_results <- full_results[, c(
        "Tile", "CellPopulation", "Foreground", "Background", "P_value", "Test_Statistic", "FDR", "Log2FC_C",
        "MeanDiff", "Avg_Intensity_Case", "Pct0_Case",
        "Avg_Intensity_Control", "Pct0_Control"
      )]

      discoveries <- sum(full_results$FDR <= qValueThreshold, na.rm = TRUE)
      if (verbose & !is.null(signalThreshold)) {
        message(
          discoveries, " differential regions found at FDR ",
          qValueThreshold
        )
      }
      DAT_list = append(DAT_list, list(full_results))
     
  }
  
  if(is.null(signalThreshold) & !all(is.na(unlist(lapply(DAT_list, function(XX) XX$FDR))))){
      message('Optimizing thresholds')
      cl <- parallel::makeCluster(numCores)
      DAT_list <- pbapply::pblapply(cl = cl, DAT_list, optimizeThreshold,
                                    qValue = qValueThreshold, minZeroDiff = minZeroDiff)
      parallel::stopCluster(cl)
      
  }else{
  
      warning('No signal threshold optimization was done, because a signalThreshold was given.')
      
  }
      
  full_results <- do.call('rbind', DAT_list)
  if (outputGRanges) {
    full_results <- MOCHA::differentialsToGRanges(full_results)
  }

  return(full_results)
}



parallelDifferential <- function(rowVals, group){
        cbind(Tile = x, estimate_differential_accessibility(sampleTileMatrix[x, ], group))
}


optimizeThreshold <- function(differentials, qValue, minZeroDiff){
        diffs = dplyr::filter(as.data.frame(differentials), !is.na(Test_Statistic))
        diffs$Avg_Intensity_Case[is.na(diffs$Avg_Intensity_Case)] = 0
        diffs$Avg_Intensity_Control[is.na(diffs$Avg_Intensity_Control)] = 0
        step1 = 0.05
        quantiles1 = seq(0.5, 0.95, by = step1)                                       
        while(step1 > 0.01){
            subDiffs <- calculateSignificantTiles(diffs, quantileNum = quantiles1, minZeroDiff = minZeroDiff)
            numSet <- unlist(lapply(subDiffs, function(XX)sum(XX$FDR < qValue)))
            names(numSet) = quantiles1
            ## Identify top 3 quantiles
            topNum = as.numeric(names(numSet)[order(numSet, decreasing = TRUE)][c(1:3)])
            ## refine quantiles for next round, and step size
            step1 = step1/2
            if(step1 > 0.01){
                quantiles1 <- seq(min(topNum), max(topNum), by = step1)    
            }
            ## repeat
        }
        finalDiffs <- subDiffs[[which(names(numSet) == topNum[1])]]
        newDiffs <- rbind(finalDiffs, diffs[!rownames(diffs) %in% rownames(finalDiffs),])
        return(newDiffs[rownames(differentials),])
}
                                    

calculateSignificantTiles <- function(diffDF, 
                                quantileNum, minZeroDiff = 0.5){

    diff0s = diffDF$Pct0_Case - diffDF$Pct0_Control
    newThresh <- stats::quantile(c(diffDF$Avg_Intensity_Case, na.rm = TRUE,
                            diffDF$Avg_Intensity_Control), probs = quantileNum)
    
    filtList = list()
    for(i in newThresh){
        filtered_res <- diffDF[diffDF$Avg_Intensity_Case > i | 
                                diffDF$Avg_Intensity_Control > i | 
                               diff0s > minZeroDiff,]

        if(length(filtered_res) > 0){
             pi0_reduced <- qvalue::pi0est(filtered_res$P_value[
                 filtered_res$P_value <= 0.95],
                  pi0.method = "bootstrap",
                  lambda = seq(0, 0.6, .05)
                )

             filtered_res$FDR <- qvalue::qvalue(filtered_res$P_value,
                                    pi0 = pi0_reduced$pi0)$qvalues
        }
        filtList <- append(filtList, list(filtered_res))
    }
    return(filtList)
}
                                    
Test_Statistic <- sampleTileMatrix <- NULL