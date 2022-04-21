#' @title \code{differential_accessibility}
#'
#' @description \code{differential_accessibility} allows you to determine whether regions of chromatin are 
#'              differentially accessible betwen groups by conducting a test 
#'
#'
#' @param tmp_mat: sample-peak matrix 
#' @param tileID: chromosome region to analyze 
#' @param group: vector of group indices (0,1)
#' @param savePlotToFile: boolean: do you want to save this file. Defaults to current directory
#'
#' @return a table indicating the results for differential accessibility
#'         which includes the following pieces of informations
#'         - Peak= Peak ID
#'         - P_value= P-value from two-part wilcoxon test
#'         - TestStatistic= statistic from two-part wilcoxon test
#'         - HL_EffectSize=hodges_lehmann effect size,
#'         - Case_mu= Median of nonzero values in case condition,
#'         - Case_rho=proportion of zeros in case condition ,
#'         - Control_mu= Median of nonzero values in control condition
#'         - Control_rho==proportion of zeros in control condition ,
#' 
#'
#' @details The technical details of the algorithm are found in XX.
#'
#' @references XX
#'
#' @export

estimate_differential_accessibility <- function(tmp_mat,tileID,group, savePlotToFile=FALSE){
    ## Filter out 
    ## 'tileID' column
    idx = which(tmp_mat$tileID == tileID)
    test_vec = as.numeric(tmp_mat[idx,2:ncol(tmp_mat)])

    ## log-transform the values 
    data=log2(test_vec+1)
   
    ## conduct two part test
    two_part_results <- TwoPart(test_vec, group=group, test='wilcoxon', point.mass=0)

    ## permutation test
    n_a = sum(group)
    n_b = length(group) - n_a

    sampled_test_vec <- sample(test_vec, size=n_a+n_b, replace=F)
    permuted_pvalue <- TwoPart(sampled_test_vec, group=group, test='wilcoxon', 
                               point.mass=0)$pvalue

    
    ## filter non-zero values for group1
    nonzero_dx <- data[group==1]
    nonzero_dx=nonzero_dx[nonzero_dx!=0]
    
    ## filter non-zero values for group0
    nonzero_control <- data[group==0]
    nonzero_control=nonzero_control[nonzero_control!=0]
    
    df = data.frame(
        Vals = data,
        Group= group,
        Group_label=ifelse(group==1, 'case','Control')
    )
    
    ## create all pairwise combinations 
    pairwise_matrix <- data.table::as.data.table(expand.grid(nonzero_dx, nonzero_control))
    pairwise_matrix$diff <- pairwise_matrix[,1] - pairwise_matrix[,2]
    
    hodges_lehmann <- median(pairwise_matrix$diff)
    res = data.frame(
        Peak=tileID,
        P_value=two_part_results$pvalue,
        TestStatistic=two_part_results$statistic,
        HL_EffectSize=hodges_lehmann,
        Case_mu=median(nonzero_dx),
        Case_rho=mean(data[group==1]==0),
        Control_mu=median(nonzero_control),
        Control_rho=mean(data[group==0]==0),
        Permuted_Pvalue = permuted_pvalue
        
        )
    

    if(savePlotToFile){
        
        fname =paste(res$Peak,'.png',sep='')
        png(fname)
        p <- ggplot(df,
                   aes(x=Vals, col=Group_label, fill=Group_label))+
               geom_histogram(aes(y=0.1*..density..),
                 alpha=0.5,position='identity',binwidth=0.1)+facet_wrap(~df$Group_label, ncol=1, scales='free')+ThemeMain+
                ggtitle(paste("Differential Accessibility",res$Peak,sep='\n'))+
        xlab('Log2 Intensity')+ylab('Frequency')+
            xlim(c(min(df$Vals)-1, max(df$Vals)+1))
        print(p)
        dev.off()
        print(paste('file saved to ',fname))
        
    }
    
    return(res)

}


find_optimal_threshold <- function(res_pvals, threshold, numCores=10){
    

    search <- seq(from=1e-5, to=1e-2, length=500)
    estimates <- sapply(search, 
           function(x)

        sum(res_pvals$Permuted_Pvalue<x) / sum(res_pvals$P_value < x)
           )
    
    max_threshold <- max(which(estimates < threshold))
    nominal_threshold <- search[max_threshold]
    return(nominal_threshold)
    
    
}

