#' @title \code{plot_sample_specific_reproducible_peaks.R}
#'
#' @description \code{plot_sample_specific_reproducible_peaks.R} allows you to
#'              visualize common peaks by setting a proportion-based threshold
#' 
#'
#' @param peak_mat a peak-sample matrix 
#' @param index a vector of 0, 1, denoting the two groups 
#' @param fname a string containing where the plots wills be saved
#'        the samples
#'
#' @return a plot indicating # of common peaks by group using 
#'         proportion-based thresholding 
#'
#'
#' @export

plot_sample_specific_reproducible_peaks.R <- function(peak_mat,
                                                      save_plot=FALSE,
                                                      fname='tmp.png'
                                                     )
    {
   
    unique_case <- unique(peak_mat$Sample[peak_mat$Index==1])
    unique_control <- unique(peak_mat$Sample[peak_mat$Index==0])

    ## calculate reproducible pekas
    CaseSamples_with_Peak <- peak_mat[Index==1, 
                                        sum(Peak)/length(unique_case), 
                                        by=Tile]
    CaseSamples_with_Peak = CaseSamples_with_Peak[V1 >0]
    
    ControlSamples_with_Peak <- peak_mat[Index==0,
                                           sum(Peak)/length(unique_control), 
                                           by=Tile]
    ControlSamples_with_Peak = ControlSamples_with_Peak[V1 >0]    

    ## Reproducible peaks 
    reproducibility_perc <- seq(0, 1, by=0.05)
    
    ## 
    Case_peak_nums<- sapply(reproducibility_perc, function(x) 
        nrow(CaseSamples_with_Peak[V1>=x])
                                         )
    Control_peak_nums <- sapply(reproducibility_perc, function(x) 
        nrow(ControlSamples_with_Peak[V1>=x])
                                           )

    ## create matrix for plotting 
    rep_df <- data.frame(nPeaks = c(Case_peak_nums, 
                                    Control_peak_nums),
                         Reproducibility=rep(reproducibility_perc,2),
                         Group = c(rep('Case',length(reproducibility_perc)),
                                   rep('Control',length(reproducibility_perc)))
               )

    if(save_plot){
        png(fname)
        p <- ggplot(rep_df,
               aes(x=Reproducibility,
                   y=nPeaks,
                   col=Group,
                   fill=Group))+geom_point() + geom_smooth(span=0.95)+ThemeMain+
            ggtitle('Peaks Common Across Samples:\n')+ xlab('Proportion of Samples')+
            ylab('# Peaks')+xlim(0,0.5)
        print(p)
        dev.off()
    }
    
    return(list(CasePeaks = CaseSamples_with_Peak, 
                ControlPeaks =ControlSamples_with_Peak)
           )
    
}