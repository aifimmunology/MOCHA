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


differential_accessibility <- function(groupA, groupB, candidatePeaks='chr1:817000-817499', doPDRAnalysis=FALSE , numCores ){
    
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
    
    if(length(sample_names_grB) < 3 | length(sample_names_grA) <3){
        
        print('not enough samples with at least 5 cells to determine differential accessibility. At least one group has less < 2 samples with > 5 cells, therefore no differential accessibility can be made.')
        
        res=data.frame(Peak = candidatePeaks,
                       ES_wilc=rep(NA, length(candidatePeaks)),
                       Wilcoxon=rep(NA, length(candidatePeaks)),
                       ES_Chisq=rep(NA, length(candidatePeaks)),
                       Chisquare = rep(NA, length(candidatePeaks)),
                       MinPval = rep(NA, length(candidatePeaks)),
                       L1A_avg=rep(NA, length(candidatePeaks)),
                       L1B_avg=rep(NA, length(candidatePeaks))
                                )
     
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


        ################################################################################################
        ################################################################################################
        ### The following two functions are 
        ### internally defined functions for 
        ### PDR Analyses which can optionally
        ### be set as "TRUE" if selected 
        ### in the parameter settings 
        
        ### lambda1-subsample determines the 
        ### value of lambda1 for a given subsample 
        ### of cells selected for a given sample
        
        get_lambda1_subsample <- function(sampleSize, ind_sample, candidatePeak, FinalBins){

                cellList <- unique(ind_sample$RG)
                subsample <- ind_sample[ind_sample$RG %in% sample(cellList,sampleSize, replace=F)]

                tmp <- scMACS::calculate_intensities(subsample,
                                              FinalBins,
                                              NBDistribution=NULL,
                                              normalizeBins=FALSE,
                                              theta=0.001, 
                                              width=20
                                             )
                tmpL1 <- round( tmp$lambda1[tmp$PeakID ==candidatePeak],4)
                print(tmpL1)
                return(tmpL1)

        }


        ### the PDR: probability detection rate 
        ### is a technical measure defined to capture
        ### dropout rate in single-cell ATAC data
        ### to better qualify whether differential
        ### accessibility is a product of technical noise 
        ### or true biological signal. 

        calculate_PDR <- function(groupB, numCellsB, vals_groupB, candidatePeak, FinalBins){

                min_num_sample <- min(numCellsB)
                most_open_sample <-  which.max(vals_groupB)

                print(paste("Smallest sample has cells =", min_num_sample))

                ind_sample <- groupB$SampleFragments[[most_open_sample]]

                lambda1_subsamples <- unlist(mclapply(1:25, function(x) 
                    try(get_lambda1_subsample(min_num_sample, ind_sample,candidatePeak,FinalBins )),
                         mc.cores=numCores,
                         mc.preschedule=TRUE,
                         mc.allow.recursive=FALSE
                         ))

                pdr = mean(lambda1_subsamples > 0)
                return( pdr)
        }

        
        ### Find differential Peaks is a 
        ### wrapper that calculates  
        
        findDifferentialPeak_pdr  <- function(candidatePeak){
        
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

        
            pdrA = calculate_PDR(groupA, numCellsA, vals_groupA ,candidatePeak, FinalBins_groupA)
            pdrB = calculate_PDR(groupB, numCellsB, vals_groupB,candidatePeak, FinalBins_groupB)    

            chisq_mat <- rbind(c(sum(vals_groupA >0), length(vals_groupA)),
                                               c(sum(vals_groupB >0), length(vals_groupB)))

            Mat <- as.table(chisq_mat)
            dimnames(Mat) <- list(Group = c("A", "B"),
                                                  Peak = c("Positive","Total"))
                      
            wilcoxon_p <- wilcox.test(vals_groupA, vals_groupB)
            chisq_p <- chisq.test(chisq_mat)

                   
            res = data.frame(

                        Peak=candidatePeak,
                        ES_wilc=wilcoxon_p$statistic,
                        Wilcoxon=wilcoxon_p$p.value,
                        ES_Chisq=chisq_p$statistic,
                        Chisquare = chisq_p$p.value,
                        MinPval = min(wilcoxon_p$p.value, chisq_p$p.value),
                        PDR_A=pdrA,
                        L1A_avg=mean(vals_groupA),
                        PDR_B= pdrB,
                        L1B_avg=mean(vals_groupB)
                )

            return(res)
        
        }
        
              
        if(doPDRAnalysis){
            
            
            FinalBins_groupA <- determine_dynamic_range(SimpleList(groupA$SampleFragments),
                                                 ArchRProj, 
                                                 binSize=500, 
                                                 doBin=FALSE)

            FinalBins_groupB <- determine_dynamic_range(SimpleList(groupB$SampleFragments),
                                                 ArchRProj, 
                                                 binSize=500, 
                                                 doBin=FALSE)
            
            
        ## If a sample had > 5 cells but no reads,
        ## lambda1 is designated as 0. So here 
        ## we change the NAs (R-artifact) to 0s 
                
            diffPeaks <- mclapply(candidatePeaks, 
                                  function(x)
                                      findDifferentialPeak_pdr(x),
                                  mc.cores = 3
                                  )

            diffPeaks <- do.call(rbind, diffPeaks)

        } else { 
            
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

                            contingency_mat <- rbind(c(sum(vals_groupA >0), length(vals_groupA)),
                                               c(sum(vals_groupB >0), length(vals_groupB)))

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

                
            diffPeaks <- mclapply(candidatePeaks, 
                              function(x)
                                  findDifferentialPeak(x),
                              mc.cores = numCores
                              )
                              
            diffPeaks <- do.call(rbind, diffPeaks)
    
            
            
            
            
        }



                              
                                      
        return(diffPeaks)                      
        
    }

}

