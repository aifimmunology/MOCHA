#' @title \code{get_differential_accessible_Regions}
#'
#' @description \code{get_differential_accessible_Regions} allows you to determine whether regions of chromatin are 
#'              differentially accessible betwen groups by conducting a test 
#'
#'

#'
#' @export
get_differential_accessible_regions <- function(Union_peaks,sample_specific_peaks,
                                               metaFile, fdr_control=0.2, nCores=2){

        ## Create sample specific matrix
        sample_peak_matrix <- create_peak_sampleMatrix(sample_specific_peaks,
                                                       Union_peaks)

        ## Get group labels 
        positive_samples <- metaFile$Sample[metaFile$Class=='Positive']
        group <- ifelse(names(sample_peak_matrix)[2:ncol(sample_peak_matrix)] %in% positive_samples,1,0)

        
        ## Estimate differential accessibility
        res_pvals <- mclapply(Union_peaks, 
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