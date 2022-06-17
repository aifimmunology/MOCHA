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

co_accessibility <- function(GR, filterPairs = NULL, index = NULL, numCores=40, verbose = FALSE){
  
  ## Prep data.table of peaks for correlations. 
  mat1 <- as.data.frame(mcols(GR))
  rownames(mat1) <- GRangesToString(GR)
  
  
  ### generate all pairwise combinations
  N = nrow(mat1)
  pairwise_combos = expand.grid(1:N, 1:N)
  pairwise_combos = pairwise_combos[!duplicated(t(apply(pairwise_combos[1:2], 1, sort))), ]
  pairwise_combos$Peak1 <- row.names(mat1)[pairwise_combos$Var1]
  pairwise_combos$Peak2 <- row.names(mat1)[pairwise_combos$Var2]
  pairwise_combos <- pairwise_combos[pairwise_combos$Peak1 != pairwise_combos$Peak2,]
  
  ## Filter out any pairs that have already been tested
  if(!is.null(filterPairs)){
    
    pairwise_combos <- pairwise_combos[!(pairwise_combos$Peak1 %in% filterPairs$Peak1 &
                                           pairwise_combos$Peak2 %in% filterPairs$Peak2),]
    
  }
  
  if(!is.null(index) & class(index) == 'numeric'){
    
    pairwise_combos <- pairwise_combos[pairwise_combos$Peak1 %in% rownames(mat1)[index] | 
                                         pairwise_combos$Peak2 %in% rownames(mat1)[index] ,]
    
  }else if(verbose){ print('No valid index given. Testing all peaks within range.')}
  
  ### Loop through all pairwise 
  ### combinations of peaks 
  
  if(nrow(pairwise_combos) == 0){return(NULL)}
  #return(mat1)
  zero_inflated_spearman <- unlist(mclapply(1:nrow(pairwise_combos),
                                            function(x)
                                              weightedZISpearman(x=mat1[pairwise_combos$Var1[x],],
                                                                 y=mat1[pairwise_combos$Var2[x],]
                                              ),
                                            mc.cores=numCores
  ))
  
  #return(zero_inflated_spearman)
  ### Create zero-inflated correlation matrix
  ### from correlation values, 
  zi_spear_mat <- data.frame(Correlation=zero_inflated_spearman,
                             Peak1= row.names(mat1)[pairwise_combos$Var1],
                             Peak2= row.names(mat1)[pairwise_combos$Var2],
                             seqnames = as.character(seqnames(GR))[pairwise_combos$Var1],
                             start = min(start(GR)[pairwise_combos$Var1],start(GR)[pairwise_combos$Var2] ),
                             end = max(end(GR)[pairwise_combos$Var1],end(GR)[pairwise_combos$Var2] )
  )
  
  return(zi_spear_mat)
  
}         




FindCoAccessibleLinks <- function(peakDT,regions, windowSize = 2*10^6, numCores = 1,  verbose = FALSE){
  
  
  tileList <- StringsToGRanges(peakDT$tileID)
  peakMat <- as.data.frame(peakDT)[,-1]
  mcols(tileList) <- peakMat
  
  wideRange <- stretch(regions, windowSize)
  
  if(length(regions) > 1){
    
    wideList <- split(wideRange,seq(1,length(wideRange)))
    window_tmp <- plyranges::filter_by_overlaps(tileList,wideList[[1]])
  }else{
    
    window_tmp <- plyranges::filter_by_overlaps(tileList,wideRange)
  }
  
  
  
  ##Run correlations for the first window
  PeakCorr <- co_accessibility(window_tmp, numCores=numCores)
  
  if(length(regions) > 1){
    for(i in 2:length(wideList)){
      #print(i)
      ## Find all peaks in the next window
      window_tmp <- plyranges::filter_by_overlaps(tileList,wideList[[i]])
      keyPeak <- queryHits(findOverlaps(window_tmp, regions[i]))
      nextCorr <- co_accessibility(window_tmp, PeakCorr, numCores=numCores,
                                   index = keyPeak, verbose = verbose)
      PeakCorr <- rbind(PeakCorr, nextCorr)
      
    }
  }
  
  return(PeakCorr)
  
}  
















