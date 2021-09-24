#' @title \code{dynamic_bins}
#'
#' @description \code{dynamic_bins} is an R helper function, part of the single-cell peak calling
#' algorithm scMACS by (XXX et al, 2021) that determines which genomic regions, or bins,
#' will be used for de-novo peak calling. The function "dynamic_bins" generates bins for peak calling based on the actual fragments present within a sample.
#AllFragmentsList: List of fragments by arrow file
#'
#'
#' @param AllFragmentsList List of fragments by arrow file
#' @param GeneralWindowSize Window size for sliding window generated over longer fragments.
#' @param WindowSizeRange The sliding window function will generate a smaller window at the end of a longer fragment,
#if the fragment is not evenly divisible by GeneralWindowSize. If that smaller window is less than or equal to
#WindowSizeRange, then it will be merged with the preceding window to generate a larger window. This means that
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
#' @export
dynamic_bins <- function(
    AllFragmentsList, 
    GeneralWindowSize = 500, 
    WindowSizeRange = 100, 
    coreNum = 1, 
    doBins = TRUE){

  all_chr <- unique(unlist(lapply(AllFragmentsList, seqlevels)))

  AllFrags <- lapply(
      frags,
      function(frag) {
          split(frag, seqnames(frag))
      }
  )
    
  AllFrags <- lapply(
    all_chr, 
    function(chr) { 
        chr_frag <- lapply(AllFrags, function(frag) {
            frag[[chr]]
        })
        plyranges::bind_ranges(chr_frag)
    })
    
  AllFrags <- mclapply(
      AllFrags, 
      function(x) {
        plyranges::reduce_ranges(x, counts = plyranges::n())
  }, mc.cores = coreNum)

  AllFrags <- plyranges::bind_ranges(AllFrags)

  if(doBins){
    print("01 - Fragments bound into one set of non-overlapping ranges")

    SmallFrags <- AllFrags[width(AllFrags) <= GeneralWindowSize + WindowSizeRange]

    BigFrags <- AllFrags[width(AllFrags) > GeneralWindowSize + WindowSizeRange] 
    BigFrags_counts <- BigFrags$counts

    BigFrags <- BigFrags %>%
      GenomicRanges::slidingWindows(width = GeneralWindowSize, step = GeneralWindowSize)

    BigFrags_splits <- diff(c(0,BigFrags@partitioning@end))
    BigFrags <- unlist(BigFrags)

    BigFrags$counts <- base::rep(BigFrags_counts, BigFrags_splits)

    print("02 - Longer fragments broken into sliding windows")
    FinalBins <- plyranges::bind_ranges(BigFrags, SmallFrags)

  } else {
    FinalBins <- AllFrags
  }

  FinalBins
}
