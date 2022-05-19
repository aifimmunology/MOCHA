#' @title \code{co_accessibility}
#'
#' @description \code{co_accessibility} allows you to determine whether 2 peaks
#'              are co-accessible using a zero-inflated spearman correlation.
#'
#'
#' @param mat: sample-peak matrix with regions to analyze
#' @param numCores: integer to determine # of paralell cores
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
#' @example 
#' Generate 
#' mat1 = matrix(pmax(0, rnorm(1000)), ncol=100); row.names(mat1) <- paste('A',1:10,sep='_')
#' ziSpear_mat <- co_accessibility(mat1, numCores=5)
#' head(ziSpear_mat)
#' @export

co_accessibility <- function(mat1, numCores=40){
    
    ### generate all pairwise combinations
    N = nrow(mat1)
    pairwise_combos = expand.grid(1:N, 1:N)
    
    ### Loop through all pairwise 
    ### combinations of peaks 
    zero_inflated_spearman <- unlist(mclapply(1:nrow(pairwise_combos),
             function(x)
                 weightedZISpearman(x=mat1[pairwise_combos$Var1[x],],
                                 y=mat1[pairwise_combos$Var2[x],]
                                 ),
             mc.cores=numCores
             ))
    
    ### Create zero-inflated correlation matrix
    ### from correlation values, 
    zi_spear_mat <- data.frame(Correlation=zero_inflated_spearman,
                               Peak1= row.names(mat1)[pairwise_combos$Var1],
                               Peak2= row.names(mat1)[pairwise_combos$Var2]
                               )
    
    return(zi_spear_mat)
    
}