#' @title \code{determine_dynamic_range}
#'
#' @description \code{determine_dynamic_range} is an R helper function, part of the single-cell peak calling
#' algorithm scMACS by (XXX et al, 2021) that determines which genomic regions, or bins,
#' will be used for de-novo peak calling. The function "determine_dynamic_range" is a wrapper function
#' around the dynamic_bins function
#'
#'
#' @param AllFragmentsList: List of fragments by arrow file
#' @param GeneralWindowSize: Window size for sliding window generated over longer fragments.
#' @param WindowSizeRange: The sliding window function will generate a smaller window at the end of a longer fragment,
#if the fragment is not evenly divisible by GeneralWindowSize. If that smaller window is less than or equal to
#WindowSizeRange, then it will be merged with the preceding window to great a larger window. This means that
#longer fragments will be broken up into bins that are between WindowSizeRange and WindowSizeRange+GeneralWindowSize in length.
#' @param coreNum an integer indicating the number of cores to use
#' @param doBins is a boolean variable. When true, then it will into windows according to GeneralWindowSize and WindowSizeRange.

#' @return a data.table, countsByBin, that returns the two intensity parameters required
#' to calculate the probability of a (+) peak
#'
#' @details The technical details of the algorithm are found in XX.
#'
#' @references XX
#'


determine_dynamic_range <- function(AllFragmentsList, ArchRProject, binSize, doBin){
  TotalRange <- DynamicBins(AllFragmentsList = AllFragmentsList,
                            doBin = FALSE,
                            coreNum = 30)

  if(class(AllFragmentsList)!='list'){
    stop('AllFragmentsList must be a list of arrow files')
  }

  if(class(ArchRProject)!='ArchRProject'){
    stop('ArchRProject must be an ArchR Project')
  }

  if(class(binSize)!='numeric' | binSize <0){
    stop(paste('invalid binSize!: binSize must be an integer value > 0 indicating the width of the genomic region check "binSize=',
               binSize,'" input', sep=''))

  }

  if(binSize > 5000){
    warning('binSize is > 5,000bp. Do you want to set the tile size smaller for peakcalling?')
  }

  if(class(doBin) != 'logical'){
    stop('doBin user-input must be a TRUE/FALSE boolean input')
  }

  blackList <- getBlacklist(ArchRProject)

  #Let's subtract out areas of fragments that overlap with blacklist regions.
  #Let's not remove the region entirely, because the regions may be long or only

  TotalRangesFilt<- setdiff_ranges(TotalRange, blackList)

  RangeBins <- plyranges::stretch(anchor_end(TotalRangesFilt),
                                  extend = start(TotalRangesFilt)%% binSize)

  FinalBins <- stretch(anchor_start(RangeBins),extend = (binSize - end(RangeBins)%% binSize)) %>%
    reduce_ranges() %>% slide_ranges(width = binSize, step = binSize)%>% filter(width(.) ==binSize)

  return(FinalBins)

  print('fragments loaded')
}
