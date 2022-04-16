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

estimate_differential_accessibility <- function(tmp_mat,tileID,group, savePlotToFile=FALSE, doPermutationTest = FALSE){
    ## Filter out 
    ## 'tileID' column
    idx = which(tmp_mat$Tile == tileID)
    test_vec = as.numeric(tmp_mat[idx,2:ncol(tmp_mat)])

    ## log-transform the values 
    data=log2(test_vec+1)
   
    ## conduct two part test
    two_part_results <- TwoPart(test_vec, group=group, test='wilcoxon', point.mass=0)

    if(doPermutationTest){
        ## permutation test
        n_a = sum(group)
        n_b = length(group) - n_a

        permutations <- mclapply(1:100000, function(x){
            sampled_test_vec <- sample(test_vec, size=n_a+n_b, replace=F)
            TwoPart(sampled_test_vec, group=group, test='wilcoxon', point.mass=0)$statistic

                },
                                 mc.cores=10
                       )

        p_value <- mean(two_part_results$statistic > permutations)
        
    }

    
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
    pairwise_matrix <- as.data.table(expand_grid(nonzero_dx, nonzero_control))
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
        Control_rho=mean(data[group==0]==0)
        
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


#' @title \code{calculate_dropout_rate}
#'
#' @description \code{differential_accessibility} allows you to determine whether regions of chromatin are 
#'              differentially accessible betwen groups by conducting a test 
#'
#'
#' @param case_frags: list of fragment files in case group
#' @param control_frags: list of fragment files in control group
#' @param case_peaks_gr: list of peaks in case group
#' @param control_peaks_gr: list of peaks in control group
#' @param numCores = number of cores to parallelize with
#' @param numReps = number of times to replicate dropout estimation
#'
#' @return median dropout rate
#' 
#'
#' @details The technical details of the algorithm are found in XX.
#'
#' @references XX
#'
#' @export


calculate_dropout_rate <- function(case_frags, control_frags,
                                   case_peaks_gr, control_peaks_gr, numCores=10, numReps=50){

    print('calculating dropout rate')
    ## calculate nFrags per sample
    nfrags_case = sapply(case_frags, length)
    nfrags_control=sapply(control_frags, length)

    median_frags_case <- median(nfrags_case)
    median_frags_control <- median(nfrags_control)

        if(median_frags_case > median_frags_control){
            
            print('case group has more fragments than control. Calculate dropout relative to control')

            

            dropout_rate <- unlist(mclapply(1:numReps, function(x){
                
                # select sample
                selected_sample <- names(sample(which(nfrags_case > median_frags_case),1))
                selected_sample_frags <- case_frags[selected_sample][[1]]
                
                # calculate dropout
                calculate_dropout_per_sample(selected_sample_frags, 
                                             control_peaks_gr,
                                             median_frags_control)
                },
                     mc.cores=numCores
                

                    ))
            return(median(dropout_rate))


        } else  if(median_frags_case < median_frags_control){
            print('control group has more fragments than case. Calculate dropout relative to case')                      
            
            dropout_rate <- unlist(mclapply(1:numReps, function(x){
                
                # select sample
                selected_sample <- names(sample(which(nfrags_control > median_frags_control),1))
                selected_sample_frags <- control_frags[selected_sample][[1]]
                
                calculate_dropout_per_sample(selected_sample_frags, 
                                             case_peaks_gr,
                                             median_frags_case)
                
                }
                                            
                                            
                                            ,
                     mc.cores=numCores

                    ))
            return(dropout_rate)
        }
        
}

#' @title \code{calculate_dropout_per_sample}
#'
#' @description \code{differential_accessibility} allows you to determine whether regions of chromatin are 
#'              differentially accessible betwen groups by conducting a test 
#'
#'
#' @param selected_sample_frags: a fragment file used to estimate dropout
#' @param peaks_gr: peakset from same condition as sample fragment file
#' @param downsample_target: # of fragments to downsample to estimate dropout
#'
#' @return dropout rate for a given iteration
#' 
#'
#' @details The technical details of the algorithm are found in XX.
#'
#' @references XX
#'
#' @export


calculate_dropout_per_sample <- function(selected_sample_frags, peaks_gr, downsample_target){

                downsample <- selected_sample_frags[sample(length(selected_sample_frags),  downsample_target)]
                overlaps <- count_overlaps(peaks_gr, downsample)
                dropout_rate <- mean(overlaps == 0)
        return(dropout_rate)

}