#' @title \code{differential_accessibility}
#'
#' @description \code{differential_accessibility} allows you to determine whether regions of chromatin are 
#'              differentially accessible betwen groups by conducting a test 
#'
#'
#' @param groupA: list of accessible peaks for each sample in group A
#' @param groupB: list of accessible peaks for each sample in group B
#' @param candidatePeaks: vector of strings containing PeakIDs (unique peak identifiers)
#'                        used to conduct differential accessibility
#' 
#' @return a table indicating the results for differential accessibility
#'         which includes the following pieces of informations
#'         - Peak= Peak ID
#'         - ES_wilc = effect size of the wilcoxon rank sum test
#'         - Wilcoxon = the statistical significance of the Rank sum test
#'         - ES_Chisq = effect size of the chi square test
#'         - Chisquare = tthe statistical significance of the chi-square test
#'         - MinPval = minimum p-value across both tests
#'         - L1A_Avg= Avg lambda1 value across samples in group A
#'         - L1B_Avg= Avg lambda1 value across samples in group B

#'
#' @details The technical details of the algorithm are found in XX.
#'
#' @references XX
#'
#' @export


differential_accessibility <- function(groupA, groupB, numCores ){
    
    ### Get the sample IDs 
    ### and remove the Union 
    ### from the groups 
    sample_names_grA <- names(groupA[['scMACS_peaks']])
    sample_names_grA <- sample_names_grA[sample_names_grA != 'Union']
    
    sample_names_grB <- names(groupB[['scMACS_peaks']])
    sample_names_grB <- sample_names_grB[sample_names_grB != 'Union']
    
    ### Error message: 
    ###   - if either group has < 2 samples without at least 5 cells
    ###     we cannot determine differential accessibility. 
    ###     as both wilcoxon and X^2 test break
    
    if(length(sample_names_grB) < 5 | length(sample_names_grA) <5){
        
        print('not enough samples with at least 5 cells to determine differential accessibility. Returning the set-difference.')
        
        ## Find the union 
        groupA_peaks <- groupA[['scMACS_peaks']]$Union
        groupB_peaks <- groupB[['scMACS_peaks']]$Union
        
        ### choose peaks in at least 2 samples
        groupA_peaks <- groupA_peaks[SamplesWithPeak>1]
        groupB_peaks <- groupB_peaks[SamplesWithPeak>1]   
        
        ### Common peaks
        commonPeaks <- intersect(groupA_peaks$PeakID, groupB_peaks$PeakID)
        
        ### Unique in A 
        ### Unique in B 
        uniqueA <- groupA_peaks[!PeakID %in% commonPeaks]
        uniqueB <- groupB_peaks[!PeakID %in% commonPeaks]
        
        
        res= list(UniquePeaksA=uniqueA, 
                  UniquePeaksB=uniqueB,
                  CommonPeaks = commonPeaks)
                                
        return(res)
        
    } else {

        ### If there are enough samples 
        ### proceed to do differential
        ### accessibility analyses 

        ### calculate # of cells 
        ### per group 

        numCellsA <- sapply(groupA[['scMACS_peaks']][sample_names_grA],
                            function(x) x$numCells[1]
                            )

        numCellsB <- sapply(groupB[['scMACS_peaks']][sample_names_grB],
                            function(x) x$numCells[1]
                            )
        
        ## Find the union 
        groupA_peaks <- groupA[['scMACS_peaks']]$Union
        groupB_peaks <- groupB[['scMACS_peaks']]$Union
        
        ### choose peaks in at least 2 samples
        groupA_peaks <- groupA_peaks[SamplesWithPeak>1]
        groupB_peaks <- groupB_peaks[SamplesWithPeak>1]  
        
        candidatePeaks <- unique(c(groupA_peaks$PeakID,
                                   groupB_peaks$PeakID)
                                 )
        
        diffPeaks <- mclapply(candidatePeaks, 
                                  function(x)
                                      findDifferentialPeak(x),
                                  mc.cores = numCores
                                  )

        diffPeaks <- do.call(rbind, diffPeaks)
       
            
        }
                                      
        return(diffPeaks)                      
        
    
}


### 
findDifferentialPeak  <- function(candidatePeak){

                            vals_groupA <- as.numeric(sapply(groupA[['scMACS_peaks']][sample_names_grA],
                                                function(x) x$lambda1[x$PeakID ==candidatePeak]
                                                )
                                                      )

                            vals_groupB <- as.numeric(sapply(groupB[['scMACS_peaks']][sample_names_grB],
                                                function(x)  x$lambda1[x$PeakID ==candidatePeak]
                                                ) 
                                                      )

                            vals_groupB[is.na(vals_groupB)] <- 0 
                            vals_groupA[is.na(vals_groupA)] <- 0 

                            A = sum(vals_groupA >0)
                            B = length(vals_groupA)-A
                            C = sum(vals_groupB >0)
                            D = length(vals_groupB) - C
                            contingency_mat <- rbind(c(A,B),
                                                     c(C,D))

                            Mat <- as.table(contingency_mat)
                            dimnames(Mat) <- list(Group = c("A", "B"),
                                                  Peak = c("Positive","Total"))
                      
                            wilcoxon_p <- wilcox.test(vals_groupA, vals_groupB)
                            fishers_test <- fisher.test(x=Mat)

                      
                            res = data.frame(
                                Peak=candidatePeak,
                                ES_wilc=wilcoxon_p$statistic,
                                Wilcoxon=wilcoxon_p$p.value,
                                ES_Fisher=fishers_test$estimate,
                                Fisher = fishers_test$p.value,
                                MinPval = min(wilcoxon_p$p.value, fishers_test$p.value),
                                L1A_avg=round(mean(vals_groupA),4),
                                L1B_avg=round(mean(vals_groupB),4)
                                )
                      
                            return(res)

                
                  
}

