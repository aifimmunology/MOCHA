#' @title \code{getDifferentialAccessibleTiles}
#'
#' @description \code{getDifferentialAccessibleTiles} allows you to
#'   determine whether regions of chromatin are differentially accessible
#'   between groups by conducting a test
#' 
#' @param SampleTileObj The SummarizedExperiment object output from getSampleTileMatrix
#' @param cellPopulation A string denoting the cell population of interest
#' @param groupColumn The column containing sample group labels
#' @param foreground The foreground group of samples for differential comparison
#' @param background The background group of samples for differential comparison
#' @param fdrToDisplay Optional false-discovery rate used only for standard output messaging 
#' @param numCores Optional, the number of cores to use with multiprocessing. Default is 1.
#' 
#' @return 
#'
#' @export


getDifferentialAccessibleTiles <- function(SampleTileObj,
					   cellPopulation,
					   groupColumn,
					   foreground,
					   background,
                                           fdrToDisplay = 0.2,
                                           numCores = 2) {
  
  if(!any(names(assays(SampleTileObj)) %in% cellPopulation)){

	stop('cellPopulation was not found within SampleTileObj. Check available cell populations with `colData(SampleTileObj)`.')
  }

  metaFile <- SummarizedExperiment::colData(SampleTileObj)
  
  if (!(groupColumn %in% colnames(metaFile))){
    stop(str_interp("Provided groupCol '{groupColumn}' not found in the provided SampleTileObj"))
  }
  if (!(foreground %in% metaFile[[groupColumn]])) {
    stop(str_interp("Provided foreground value is not present in the column {groupColumn} in the provided SampleTileObj"))
  }
  if (!(background %in% metaFile[[groupColumn]])) {
    stop(str_interp("Provided background value is not present in the column {groupColumn} in the provided SampleTileObj"))
  }

  # Get group labels
  positive_samples <- metaFile[metaFile[,groupColumn] == foreground, 'Sample']
  negative_samples <- metaFile[metaFile[,groupColumn] == background, 'Sample']


  tilesCalled <- mcols(SummarizedExperiment::rowRanges(SampleTileObj))[,cellPopulation]
  sampleTileMatrix <- scMACS:::getCellTypeMatrix(SampleTileObj, cellPopulation)
  sampleTileMatrix <- sampleTileMatrix[, which(colnames(SampleTileObj) %in% c(positive_samples, negative_samples))]
    
  # We need to enforce that NAs were set to zeros in getSampleTileMatrix
  miscMetadata <- S4Vectors::metadata(SampleTileObj)
  if (!miscMetadata$NAtoZero){
    #stop("The provided SampleTileObj was generated with NAtoZero=FALSE. Please run getSampleTileMatrix with NAtoZero=TRUE For differential analysis.")
    sampleTileMatrix[is.na(sampleTileMatrix)] = 0
  }
  group <- as.numeric(colnames(sampleTileMatrix) %in% positive_samples)

  #############################################################################
  # Prioritize high-signal tiles

  medians_a <- matrixStats::rowMedians(sampleTileMatrix[,which(group==1)], na.rm=T)
  medians_b <- matrixStats::rowMedians(sampleTileMatrix[,which(group==0)], na.rm=T)
    
  zero_A <- rowMeans(is.na(sampleTileMatrix[,which(group==1)]))
  zero_B <- rowMeans(is.na(sampleTileMatrix[,which(group==0)]))

  diff0s = abs(zero_A-zero_B)

  Log2Intensity = metadata(SampleTileObj)$Log2Intensity
  log2FC_filter = ifelse(Log2Intensity, 12, 2^12)
  idx = which(medians_a > log2FC_filter| medians_b > log2FC_filter | diff0s >= 0.5)

  ############################################################################    
  # Estimate differential accessibility


  res_pvals <- mclapply(rownames(sampleTileMatrix),
    function(x) {
      cbind(Tile= x, estimate_differential_accessibility(sampleTileMatrix[x,], group, F))
    },
    mc.cores = numCores
  )

  # Combine results into single objects
  res_pvals <- rbindlist(res_pvals)

  #############################################################################
  # Apply FDR on filtered regions           

  filtered_res = res_pvals[idx,]
  
  
  pi0_reduced = qvalue::pi0est(filtered_res$P_value[filtered_res$P_value <= 0.95],
                                     pi0.method = 'bootstrap',
					lambda = seq(0,0.6,.05))

  filtered_res$FDR = qvalue::qvalue(filtered_res$P_value, pi0 = pi0_reduced$pi0)$qvalues
 
  #############################################################################
    
  # Join with original results 
  full_results = dplyr::left_join(res_pvals, filtered_res[,c('Tile','FDR')], by='Tile')
  na.idx = which(is.na(full_results$FDR))
  full_results$P_value[na.idx] <- NA
  full_results$TestStatistic[na.idx] <- NA

  sampleTileMatrix[is.na(sampleTileMatrix)] <- 0
  meansA = rowMeans(sampleTileMatrix[,group==1], na.rm=T)
  meansB = rowMeans(sampleTileMatrix[,group==0])
  full_results$MeanDiff = meansA-meansB
  full_results$CellPopulation = rep(cellPopulation, length(meansA))
  full_results$Foreground = rep(foreground, length(meansA))
  full_results$Background = rep(background, length(meansA))


  colnames(full_results) <- c('Tile','P_value','Test-Statistic','Log2FC_C',
                                    'MeanDiff','Avg_Intensity_Case','Pct0_Case',
                                    'Avg_Intensity_Control','Pct0_Control','FDR',
				     'CellPopulation','Foreground','Background')
    
  full_results = full_results[, c('Tile','CellPopulation','Foreground','Background', 'P_value','Test-Statistic','FDR','Log2FC_C',
                                    'MeanDiff','Avg_Intensity_Case','Pct0_Case',
                                    'Avg_Intensity_Control','Pct0_Control')]
  

  discoveries <- sum(full_results$FDR <= 0.2, na.rm = TRUE)
  message(
    discoveries, " differential regions found at FDR ",
    fdrToDisplay
  )
  
  full_results

}
