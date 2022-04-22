#' @title \code{plot_differential_region}
#'
#' @description \code{plot_differential_region} allows you to graph a particular
#'              region across two different conditions
#'
#' @param res_pvals: output from differential accessibility
#' @param threshold: fdr threshold to control for 
#'
#' 
#' @return nominal_threshold: the largest nominal threshold 
#'         that controls for false discovery at the user-defined
#'         threshold.
#'
#'
#' @export   

find_alpha_threshold <- function(res_pvals, threshold){
    
    ## create grid search 
    alpha_search <- seq(from=1e-6, to=5e-2, length=1000)
    
    ## calculate rate of discoveries
    ## in permuted vs. true p-values 
    estimates <- sapply(alpha_search, 
           function(x)

        sum(res_pvals$Permuted_Pvalue<x) / sum(res_pvals$P_value < x)
           )
    
    ## identify maximum p-value that
    ## controls for false discoveries
    max_threshold <- max(which(estimates < threshold))
    
    ## return nominal threshold
    nominal_threshold <- alpha_search[max_threshold]
    return(nominal_threshold)
    
    
}

