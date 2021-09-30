# Estimating Maximum Helper functions 
#  setwd('/home/jupyter/variousCellTypes/')
# load(dir()[1])
# mat = countsMatrixList[[180]]
estimateEmpiricalDistribution <- function(mat){
  
  df = data.frame(round(table(mat$maxIntensity) / nrow(mat),6))
  colnames(df) <- c('MaxInt','Probability')
  return(df)
  
}

max_poisson_mass_function <- function(lambda, n, k){
  
  f_x = dpois(x=k, lambda=lambda)
  F_x = ppois(q=k, lambda=lambda)
  nMinus1 =n-1
  prob = n* (F_x^(nMinus1)) * f_x
  
  return(round(prob,6))
  
}      

max_nbinom_mass_function <- function(mean.param, dispersion.param , n,k){
  
  f_x <- dnbinom(k, size=dispersion.param, mu = mean.param)
  F_x <- pnbinom(k, size=dispersion.param, mu = mean.param)
  nMinus1 =n-1
  
  prob = n*(F_x^(nMinus1)) * f_x
  
  return(round(prob,5))
  
}

obtainTheoreticalEstimates_NB <- function(mean.param, dispersion.param, n,maxObserved){
  
  values = c(0:maxObserved)
  probDistn <- sapply(values, function(K) round(max_nbinom_mass_function(mean.param = mean.param,
                                                                         dispersion.param=dispersion.param, 
                                                                         n=n,
                                                                         k=K),5)
  )
  # Divide by normalizing constant 
  probDistn <- probDistn / sum(probDistn)
  
  theoreticalProbs <- data.frame(Probability =probDistn,
                                 MaxInt = values
  )
  
  return(theoreticalProbs)
}

calcKL <- function(empiricalDist, theoreticalDistList, i){
  probMat <- rbind(empiricalDist$Probability,
                   theoreticalDistList[[i]]$Probability)
  
  philentropy::KL(probMat)
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

calculateNBDistribution <- function(mat, theta=0.1){
  maxObserved <- max(mat$maxIntensity)
  
  n = mat$numCells[1]
  empiricalDist <- estimateEmpiricalDistribution(mat)
  
  NBDistribution <- sapply(0:maxObserved, 
                           function(x) calculateLambda(mat,x, fixedTheta=theta)
  )
  
  
  return(NBDistribution)
}


