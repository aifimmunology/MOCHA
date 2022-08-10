#' @title \code{get_differential_accessible_Regions}
#'
#' @description \code{get_differential_accessible_Regions} allows you to
#'   determine whether regions of chromatin are differentially accessible
#'   between groups by conducting a test
#'
#'
#'
#' @export
getDifferentialAccessibleTiles <- function(SampleTileObj,
					   cellType,
					   groupColumn,
					   ForeGround,
					   BackGround,
                                           fdr_control = 0.2,
                                           numCores = 2) {
  if(!any(names(assays(SampleTileObj)) %in% cellType)){

	stop('Error: Cell type not found within Summarized Experiment object. Check cell type name.')
  }

  metaFile <- SummarizedExperiment::colData(SampleTileObj)



  ## Get group labels
  positive_samples <- metaFile[metaFile[,groupColumn] == ForeGround, 'Sample']
  negative_samples <- metaFile[metaFile[,groupColumn] == BackGround, 'Sample']

  ## Estimate differential accessibility
  peaksCalled <- mcols(rowRanges(SampleTileObj))[,cellType]
  sampleTileMatrix <- SummarizedExperiment::assays(SampleTileObj)[[cellType]][peaksCalled,
				colnames(SampleTileObj) %in% c(positive_samples, negative_samples)]
  
  sampleTileMatrix[is.na(sampleTileMatrix)] <- 0

  group <- as.numeric(colnames(sampleTileMatrix) %in% positive_samples)

  #############################################################################
  ## prioritize high-signal tiles

  medians_a <- matrixStats::rowMedians(sampleTileMatrix[,which(group==1)], na.rm=T)
  medians_b <- matrixStats::rowMedians(sampleTileMatrix[,which(group==0)], na.rm=T)
    
  zero_A <- rowMeans(is.na(sampleTileMatrix[,which(group==1)]))
  zero_B <- rowMeans(is.na(sampleTileMatrix[,which(group==0)]))

  diff0s = abs(zero_A-zero_B)

  Log2Intensity = metadata(SampleTileMatrices)$Log2Intensity
  log2FC_filter = ifelse(Log2Intensity, 12, 2^12)
  idx = which(medians_a > log2FC_filter| medians_b > log2FC_filter | diff0s >= 0.5)

 #############################################################################    
  ## Estimate differential accessibility


  res_pvals <- mclapply(rownames(sampleTileMatrix),
    function(x) {
      cbind(Peak= x, estimate_differential_accessibility(sampleTileMatrix[x,], group, F))
    },
    mc.cores = numCores
  )

  ## combine results into single objects
  res_pvals <- rbindlist(res_pvals)

  #############################################################################
  ## Apply FDR on filtered regions           

  filtered_res = res_pvals[idx,]
  
  
  pi0_reduced = qvalue::pi0est(filtered_res$P_value[filtered_res$P_value <= 0.95],
                                     pi0.method = 'bootstrap',
					lambda = seq(0,0.6,.05))

  filtered_res$FDR = qvalue::qvalue(filtered_res$P_value, pi0 = pi0_reduced$pi0)$qvalues
 
  #############################################################################
    
  ## Join with original results 
  full_results = dplyr::left_join(res_pvals, filtered_res[,c('Peak','FDR')], by='Peak')
  na.idx = which(is.na(full_results$FDR))
  full_results$P_value[na.idx] <- NA
  full_results$TestStatistic[na.idx] <- NA

  sampleTileMatrix[is.na(sampleTileMatrix)] <- 0
  meansA = rowMeans(sampleTileMatrix[,group==1], na.rm=T)
  meansB = rowMeans(sampleTileMatrix[,group==0])
  full_results$MeanDiff = meansA-meansB
  full_results$CellType = rep(cellType, length(meansA))
  full_results$ForeGround = rep(ForeGround, length(meansA))
  full_results$BackGround = rep(BackGround, length(meansA))


  colnames(full_results) <- c('Peak','P_value','Test-Statistic','Log2FC_C',
                                    'MeanDiff','Avg_Intensity_Case','Pct0_Case',
                                    'Avg_Intensity_Control','Pct0_Control','FDR',
				     'CellType','ForeGround','Background')
    
  full_results = full_results[, c('Peak','CellType','ForeGround','Background', 'P_value','Test-Statistic','FDR','Log2FC_C',
                                    'MeanDiff','Avg_Intensity_Case','Pct0_Case',
                                    'Avg_Intensity_Control','Pct0_Control')]

                                    
  #############################################################################

  discoveries <- sum(full_results$FDR <= 0.2, na.rm = TRUE)
  print(paste(
    discoveries, " differential regions found at FDR",
    fdr_control
  ))
  return(full_results)

}
