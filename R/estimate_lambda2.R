#' @title \code{Helper Functions}
#'
#' @description \code{Helper Functions} is a set of helper functions used to determine maximum 
#'              intensity transformations into an estimate of lambda. These helper functions 
#'              are called internally and not meant to be called directly outside of the 
#'              calculate_intensities function call. 
#'
#'
#' @details The technical details of the algorithm are found in XX.
#'
#' @references XX
#'

max_nbinom_mass_function <- function(mean.param, dispersion.param , n,k){
  
  f_x <- dnbinom(k, size=dispersion.param, mu = mean.param)
  F_x <- pnbinom(k, size=dispersion.param, mu = mean.param)
  nMinus1 =n-1
  
  prob = n*(F_x^(nMinus1)) * f_x
  
  return(round(prob,5))
  
}


calculateLambda <- function(mat, observedMaxInt,fixedTheta){
  
  lowerLimit <- 0
  upperLimit <- 1   
  n = mat$numCells[1]
  
  lambdaSearch <- seq(lowerLimit, upperLimit, by=.0001)
  
  
  probsNB <- sapply(lambdaSearch, 
                    function(x) max_nbinom_mass_function(x, 
                                                         dispersion.param = fixedTheta,
                                                         n = n, 
                                                         k = observedMaxInt)
  )
  
  return(lambdaSearch[which.max(probsNB)])                  
  
  
}

calculateNBDistribution <- function(mat, theta=0.001){
  maxObserved <- max(mat$maxIntensity)
  
  n = mat$numCells[1]  
  NBDistribution <- sapply(0:maxObserved, 
                           function(x) calculateLambda(mat,x, fixedTheta=theta)
  )
  
  
  return(NBDistribution)
}


