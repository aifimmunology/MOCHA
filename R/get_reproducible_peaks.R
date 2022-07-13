#' @title \code{get_reproducible_peaks}
#'
#' @description \code{get_reproducible_peaks} is an R helper function, part of the single-cell peak 
#'

get_reproducible_peaks <- function(sample_specific_peaks, metaFile, 
                                   reproducibility=0.2, 
                                   fname='cd16_peaks.png'){
    
    
    ## Extract Peak list from Sample Specific Peaks Object 
    sample_names <-  sapply(sample_specific_peaks, function(x) names(x))


    ## Reformat peak information to 
    ## obtain the # of peaks reproducible
    ## across 0-50% of the samples
    lambda1 <- lapply(1:length(sample_specific_peaks),
                      function(x) 
                          data.table(Sample=sample_names[x],
                                     Intensity=sample_specific_peaks[[x]][[1]]$TotalIntensity,
                                     tileID=sample_specific_peaks[[x]][[1]]$tileID,
                                     Peak=sample_specific_peaks[[x]][[1]]$peak,
                                     Index=ifelse(metaFile$Class[x]=='Positive',1,0)
                                    )
                      )
    lambda_mat <- rbindlist(lambda1)

    #######################################################################
    #######################################################################
    # Plot sample-specific reproducible peaks

    counts_by_group <- plot_sample_specific_reproducible_peaks(lambda_mat,
                                                                 save_plot=TRUE,
                                                                 fname=fname)
    
    ### plot Venn diagram 
    Repro_Limit = reproducibility
    reproducible_casePeaks = counts_by_group$CasePeaks[V1>= Repro_Limit]
    reproducible_controlPeaks = counts_by_group$ControlPeaks[V1>= Repro_Limit]

    # Generate the two sets for 
    # the venn diagramß
    set1 <- reproducible_casePeaks$tileID
    set2 <- reproducible_controlPeaks$tileID
    common_peaks <- intersect(set1,set2)
    union_peaks <- union(set1, set2)
    
    results <- list(CasePeaks=set1,
                    ControlPeaks=set2,
                    Intersection=common_peaks,
                    Union=union_peaks)
                    
    return(results)
}