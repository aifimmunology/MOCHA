#' @title \code{calculate_local_noise}
#'
#' @description \code{calculate_local_noise} is an R helper function, part of the single-cell peak calling algorithm scMACS by (XXX et al, 2021) that determines the "background noise" located 
#' around each candidate peak region.  
#'
#'
#' @param lambda1 a vector, intensity (reads per cell) from which to calculate the 
#'        background noise in a neighborhood around a given region
#' @param width integer, an positive integer indicating how many bins around region j
#'        to calculate background noise. For example, for 500-bp regions, width=20 will 
#'.       calculate noise-level in a 10kbp neighborhood. 
#'
#' @return local_noise_estimates a vector length(lambda1) containing noise estimates 
#'         around each value of lambda1. For windows at the boundaries with < Width/2 
#'         neighbors either to the right or left, a 1-sided background noise estimate
#'         is provided using a one-sided search. 
#'
#'
#' @references XX
#'

calculate_local_noise <- function(lambda1, width=20){
    
    if(class(lambda1)!='numeric'){
        stop('lambda1 must be a numeric vector containing ')
    }

    if(class(width)!='numeric' | width %%1 !=0 | width < 0){
        stop('width must be a positive integer indicating the number of windows to search to calculate local noise estimates of lambda1')
    }
    
    ### set internal parameters 
    N = length(lambda1)
    half_width = round(width/2)
        
    ### use the roll mean function to calculate
    ### a centered local estimate of noise 
    local_noise_estimates <- zoo::rollmean(lambda1, width,  na.rm=TRUE, fill=NA, align='center')
    
    ### boundary cases that are < width/2 to the ends of the
    ### genome will have NA values. This provides a 1-sided 
    ### estimate of noise searching width/2 tiles left and right 
    ### of the lambda1 estimates 
    
    left_boundary <- rollmean(lambda1, half_width,  na.rm=TRUE, fill=NA, align='left')
    right_boundary <- rollmean(lambda1, half_width,  na.rm=TRUE, fill=NA, align='right')

    ### replace boundary cases with half-width estimates 
    ### of noise 
    local_noise_estimates[1:half_width] <- left_boundary[1:half_width]
    local_noise_estimates[(N-half_width):N] <- right_boundary[(N-half_width):N]
    
    return(local_noise_estimates)
    
}
