#' @title \code{differential_accessibility}
#'
#' @description \code{differential_accessibility} allows you to determine whether regions of chromatin are 
#'              differentially accessible betwen groups by conducting a test 
#'
#'
#' @param groupA: list of accessible peaks for each sample in group A
#' @param groupB: list of accessible peaks for each sample in group B
#' 
#' @return a genomic ranges object that 
#' to calculate the probability of a (+) peak
#'
#' @details The technical details of the algorithm are found in XX.
#'
#' @references XX
#'

differential_accessibility <- function(groupA, groupB, candidatePeak='chr1:817000-817499' ){
    
    
    sample_names_grA <- names(groupA[['scMACS_peaks']])
    sample_names_grA <- sample_names_grA[sample_names_grA != 'Union']
    
    sample_names_grB <- names(groupB[['scMACS_peaks']])
    sample_names_grB <- sample_names_grB[sample_names_grB != 'Union']
    
    if(length(sample_names_grB) < 2 | length(sample_names_grA) <2){
        
        print('not enough samples with at least 5 cells to determine differential accessibility')
         
        res = data.frame(
            Peak=candidatePeak,
            ES=NA,
            P.value=NA,
            PDR_A=NA,
            L1A_avg=NA,
            PDR_B= NA,
            L1B_avg=NA
        )
     
        return(res)
        
    } else {
    
        numCellsA <- sapply(groupA[['scMACS_peaks']][sample_names_grA],
                            function(x) x$numCells[1]
                            )

        numCellsB <- sapply(groupB[['scMACS_peaks']][sample_names_grB],
                            function(x) x$numCells[1]
                            )


        vals_groupA <- as.numeric(sapply(groupA[['scMACS_peaks']][sample_names_grA],
                            function(x) x$lambda1[x$PeakID ==candidatePeak]
                            )
                                  )

        vals_groupB <- as.numeric(sapply(groupB[['scMACS_peaks']][sample_names_grB],
                            function(x)  x$lambda1[x$PeakID ==candidatePeak]
                            ) 
                                  )

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
        vals_groupB[is.na(vals_groupB)] <- 0 
        vals_groupA[is.na(vals_groupA)] <- 0 


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


        pdrA = calculate_PDR(groupA, numCellsA, vals_groupA ,candidatePeak, FinalBins_groupA)
        pdrB = calculate_PDR(groupB, numCellsB, vals_groupB,candidatePeak, FinalBins_groupB)    


        res = data.frame(
            Peak=candidatePeak,
            ES=wilcox.test(vals_groupA, vals_groupB)$statistic,
            P.value=wilcox.test(vals_groupA, vals_groupB)$p.value,
            PDR_A=pdrA,
            L1A_avg=mean(vals_groupA),
            PDR_B= pdrB,
            L1B_avg=mean(vals_groupB)
            )

        return(res)
        
        
    }

}

