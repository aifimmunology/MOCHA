#' @title \code{co_accessibility}
#'
#' @description \code{co_accessibility} allows you to determine whether 2 peaks
#'              are co-accessible using a zero-inflated spearman correlation.
#'
#'
#' @param mat: sample-peak matrix with regions to analyze
#' @param numCores: integer to determine # of paralell cores
#'
#' @return a 3-column data.frame containing 
#'         - Correlation= Zero-inflated Spearman Correlation
#'         - Peak1= location of co-accessible region 1
#'         - Peak2= location of co-accessible region 2
#' 
#' @details The technical details of the zero-inflated correlation can be
#'          found here: 
#'
#'               Pimentel, Ronald Silva, "Kendall's Tau and Spearman's Rho 
#'               for Zero-Inflated Data" (2009). Dissertations. 
#'
#'          while the implementation (scHOT R package), can be found here:
#'               http://www.bioconductor.org/packages/release/bioc/html/scHOT.html
#'
#' @references XX
#' @example 
#' Generate 
#' mat1 = matrix(pmax(0, rnorm(1000)), ncol=100); row.names(mat1) <- paste('A',1:10,sep='_')
#' ziSpear_mat <- co_accessibility(mat1, numCores=5)
#' head(ziSpear_mat)
#' @export

co_accessibility <- function(subMat, filterPairs = NULL, index = NULL, numCores=40, verbose = FALSE){
  
  ## Prep data.table of peaks for correlations. 
  ##mat1 <- as.data.table(mcols(GR))
  ##rownames(mat1) <- GRangesToString(GR)
  
  #mat1 <- subMat[,-c('tileID','chr', 'start', 'end')]
  #rownames(mat1) = subMat$tileID
  pairwise_combos = RcppAlgos::comboGrid(rownames(subMat),rownames(subMat), repetition = FALSE)
  
  ## Filter out any pairs that have already been tested
  if(!is.null(filterPairs)){
    
    pairwise_combos <- pairwise_combos[!(pairwise_combos[,'Var1'] %in% filterPairs$Peak1 &
                                           pairwise_combos[,'Var2'] %in% filterPairs$Peak2),]
    
  }
  
  if(!is.null(index) & class(index) =='integer'){
    
    pairwise_combos <- pairwise_combos[pairwise_combos[,'Var1'] %in% rownames(subMat)[index] | 
                                         pairwise_combos[,'Var2'] %in% rownames(subMat)[index] ,]
    
  }else if(verbose){ print('No valid index given. Testing all peaks within range.')}
  
  
  ### Loop through all pairwise 
  ### combinations of peaks 
  
  
  ## nrows hits an error when there's only one pair of peaks that needs to be tested. 
  if(length(pairwise_combos) == 0 & verbose){
    print('Warning: No peaks found in neighborhood of region of interest')
    return(NULL)
  }else if(length(pairwise_combos) == 0) { return(NULL)}
  else if(length(pairwise_combos) == 2) { 
    
    #If only one pair of peaks to test, then it's no longer a data.frame, but a vector.
    zero_inflated_spearman <- weightedZISpearman(x= subMat[pairwise_combos[1],],
                                                 y= subMat[pairwise_combos[2],])
    
    zi_spear_mat <- data.table(Correlation=zero_inflated_spearman,
                               Peak1= pairwise_combos[1],
                               Peak2= pairwise_combos[2])
    
    
  }else{
  
    zero_inflated_spearman <- unlist(mclapply(1:nrow(pairwise_combos),
                                            function(x)
                                              weightedZISpearman(x= subMat[pairwise_combos[x,'Var1'],],
                                                                 y= subMat[pairwise_combos[x,'Var2'],]
                                              ),
                                            mc.cores=numCores
                                    ))
    
    
    #return(zero_inflated_spearman)
    ### Create zero-inflated correlation matrix
    ### from correlation values, 
    zi_spear_mat <- data.table(Correlation=zero_inflated_spearman,
                               Peak1= pairwise_combos[,'Var1'],
                               Peak2= pairwise_combos[,'Var2']
    )
    
  }
  
  return(zi_spear_mat)
  
}         




FindCoAccessibleLinks <- function(peakDT,regions, windowSize = 2*10^6, numCores = 1,  verbose = FALSE){
  
  peakDF <- as.data.frame(peakDT[,-'tileID'])
  rownames(peakDF) <- peakDT$tileID
  
  start <- as.numeric(gsub('chr.*\\:|\\-.*','',peakDT$tileID))
  end <- as.numeric(gsub('chr.*\\:|.*\\-','',peakDT$tileID))
  chr <- gsub('\\:.*','',peakDT$tileID)
  
  regionDF <- as.data.frame(regions)
 
  #Initialize Correlation Datatable. 
  PeakCorr <- NULL
  pb = txtProgressBar(min = 0, max = length(wideList), initial = 0)  
  
  for(i in 1:length(regions)){
      print(i)
      setTxtProgressBar(pb,i)
      
      ## Find all peaks in the next window
      window_tmp <- which(start > regionDF$start[i] - 2*10^6 &
                           end < regionDF$end[i] + 2*10^6 & 
                           chr == regionDF$seqnames[i])
      
      if(length(window_tmp) > 1){
        
        #The region of interest should overlap with the peak at least partially.
        keyPeak <- which(
                          window_tmp == which(
                            ((start <= regionDF$start[i] & end >= regionDF$start[i]) |
                             (start >= regionDF$start[i] & end<= regionDF$end[i])  |
                              (start <= regionDF$end[i] & end >= regionDF$end[i])) & 
                           chr %in% regionDF$seqnames[i]))
        
        nextCorr <- co_accessibility(peakDF[window_tmp,],  filterPairs = PeakCorr, numCores=numCores,
                                     index = keyPeak, verbose = verbose)
        PeakCorr <- rbind(PeakCorr, nextCorr)
        
      }else if(verbose){
        print('Warning: No correlations found for given region.')
      }
      
  }

  close(pb)
  gc()
  return(PeakCorr)
  
}  

FilterCoAccessibleLinks <- function(PeakCorr, threshold = 0.5){
  
  if(!any(abs(PeakCorr$Correlation) > threshold)){
    stop('Error: There are no values above the threshold.')  
  
  }
  
  PeakCorr1 <- PeakCorr[abs(PeakCorr$Correlation) > threshold,]
  start1 <- as.numeric(gsub('chr.*\\:|\\-.*','',  PeakCorr1$Peak1))
  end1 <- as.numeric(gsub('chr.*\\:|.*\\-','', PeakCorr1$Peak1))
  
  start2 <- as.numeric(gsub('chr.*\\:|\\-.*','', PeakCorr1$Peak2))
  end2 <- as.numeric(gsub('chr.*\\:|.*\\-','', PeakCorr1$Peak2))
  
  PeakCorr1$chr <- gsub('\\:.*','',PeakCorr1$Peak2)
  PeakCorr1$start <-  apply(data.table(start1, start2), 1, min)
  PeakCorr1$end <-  apply(data.table(end1, end2), 1, max)

  return(PeakCorr1)
}


#############################################












