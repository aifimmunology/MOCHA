#' @title \code{get_differential_accessible_Regions}
#'
#' @description \code{get_differential_accessible_Regions} allows you to determine whether regions of chromatin are 
#'              differentially accessible betwen groups by conducting a test 
#'
#'

#'
#' @export
get_differential_accessible_regions <- function(sample_peak_matrix,
                                               metaFile, fdr_control=0.2, nCores=2){


        ## Get group labels 
        positive_samples <- metaFile$SampleCellType[metaFile$Class=='Positive']
        group <- ifelse(names(sample_peak_matrix)[2:ncol(sample_peak_matrix)] %in% positive_samples,1,0)
        
        ## Estimate differential accessibility
        res_pvals <- mclapply(sample_peak_matrix$tileID, 
           function(x) estimate_differential_accessibility(sample_peak_matrix,x,group, F),
                 mc.cores=nCores
           )

        ## combine results into single objects
        res_pvals = rbindlist(res_pvals)

        ## 
        thresh <- find_alpha_threshold(res_pvals, fdr_control)

        res_pvals$FDR_Controlled <- res_pvals$P_value < thresh 

        discoveries <- sum(res_pvals$FDR_Controlled)
        print(paste(discoveries, ' differential regions found at FDR',
                    fdr_control))
        return(res_pvals)

}