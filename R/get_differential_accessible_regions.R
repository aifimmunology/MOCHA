#' @title \code{get_differential_accessible_Regions}
#'
#' @description \code{get_differential_accessible_Regions} allows you to determine whether regions of chromatin are 
#'              differentially accessible betwen groups by conducting a test 
#'
#'

#'
#' @export
get_differential_accessible_regions <- function(sample_peak_matrix,
                                               metaFile, nCores, seed=1){


        ## Get group labels 
        positive_samples <- metaFile$SampleCellType[metaFile$Class=='Positive']
        group <- ifelse(names(sample_peak_matrix)[2:ncol(sample_peak_matrix)] %in% positive_samples,1,0)
    
    
        #############################################################################
        ## prioritize high-signal tiles
        mat = as.matrix(sample_peak_matrix[, 2:ncol(sample_peak_matrix)])
        mat = log2(mat+1)
        mat[mat==0] <- NA
        
     
        medians_a = matrixStats::rowMedians(mat[,which(group==1)], na.rm=T)
        medians_b = matrixStats::rowMedians(mat[,which(group==0)], na.rm=T)
    
        zero_A <- rowMeans(is.na(mat[,which(group==1)]))
        zero_B <- rowMeans(is.na(mat[,which(group==0)]))
    
        diff0s = abs(zero_A-zero_B)
    
        ## Prioritize Regions 
        log2FC_filter =12
        idx = which(medians_a > log2FC_filter| medians_b > log2FC_filter | diff0s >= 0.5)
        
        #############################################################################    
        ## Estimate differential accessibility
        res_pvals <- mclapply(sample_peak_matrix$tileID, 
           function(x) estimate_differential_accessibility(sample_peak_matrix,x,group, F),
                 mc.cores=nCores
           )

        ## combine results into single objects
        res_pvals = rbindlist(res_pvals)
       
        #############################################################################
        ## Apply FDR on filtered regions             
        filtered_res = res_pvals[idx,]
            
        pi0_reduced = qvalue::pi0est(filtered_res$P_value[filtered_res$P_value <= 0.95],
                                     pi0.method = 'bootstrap', 
             lambda = seq(0,0.6,.05))
        filtered_res$FDR = qvalue::qvalue(filtered_res$P_value, 
                                          pi0 = pi0_reduced$pi0)$qvalues
        #############################################################################
    
        ## Join with original results 
        full_results = left_join(res_pvals, filtered_res[,c('Peak','FDR')], by='Peak')
        na.idx = which(is.na(full_results$FDR))
        full_results$P_value[na.idx] <- NA
        full_results$TestStatistic[na.idx] <- NA
    
        ## 
        mat[is.na(mat)] <- 0
        meansA = rowMeans(mat[,group==1], na.rm=T)
        meansB = rowMeans(mat[,group==0])
        full_results$MeanDiff = meansA-meansB
    
        colnames(full_results) <- c('Peak','P_value','Test-Statistic','Log2FC_C',
                                    'MeanDiff','Avg_Intensity_Case','Pct0_Case',
                                    'Avg_Intensity_Control','Pct0_Control','FDR')
    
        full_results = full_results[, c('Peak','P_value','Test-Statistic','FDR','Log2FC_C',
                                    'MeanDiff','Avg_Intensity_Case','Pct0_Case',
                                    'Avg_Intensity_Control','Pct0_Control')]
                                    
        #############################################################################
        return(full_results)

}