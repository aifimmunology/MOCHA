#' @title \code{dynamic_bins}
#'
#' @description \code{dynamic_bins} is an R helper function, part of the single-cell peak calling
#' algorithm MOCHA by (Zaim, Pebworth, et. al. 2022) that determines which genomic regions, or bins,
#' will be used for de-novo peak calling. The function "dynamic_bins" generates bins for peak calling based on the actual fragments present within a sample.
# AllFragmentsList: List of fragments by arrow file
#'
#'
#' @param AllFragmentsList: List of fragments by arrow file
#' @param GeneralWindowSize: Window size for sliding window generated over longer fragments.
#' @param WindowSizeRange: The sliding window function will generate a smaller window at the end of a longer fragment,
# if the fragment is not evenly divisible by GeneralWindowSize. If that smaller window is less than or equal to
# WindowSizeRange, then it will be merged with the preceding window to great a larger window. This means that
# longer fragments will be broken up into bins that are between WindowSizeRange and WindowSizeRange+GeneralWindowSize in length.
#' @param coreNum an integer indicating the number of cores to use
#' @param doBins is a boolean variable. When true, then it will into windows according to GeneralWindowSize and WindowSizeRange.

#' @return a data.table, countsByBin, that returns the two intensity parameters required
#' to calculate the probability of a (+) peak
#'
#'
#' @noRd
#'


dynamic_bins <- function(AllFragmentsList, GeneralWindowSize, WindowSizeRange, coreNum, doBin) {
  AllFragsList <- plyranges::reduce_ranges(
    AllFragmentsList,
    counts = plyranges::n()
  )
  AllFrags <- plyranges::bind_ranges(AllFragsList) %>% plyranges::reduce_ranges(counts = sum(counts))

  if (doBin) {
    message("01 - Fragments bound into one set of non-overlapping ranges")

    SmallFrags <- AllFrags[width(AllFrags) <= GeneralWindowSize + WindowSizeRange]

    BigFrags <- AllFrags[width(AllFrags) > GeneralWindowSize + WindowSizeRange] %>%
      plyranges::slide_ranges(width = GeneralWindowSize, step = GeneralWindowSize) %>%
      anchor_end() %>%
      mutate(partition = (as.character(partition))) %>%
      mutate(width = ifelse(width <= WindowSizeRange, width + 2, width - 1)) %>%
      group_by(partition) %>%
      plyranges::reduce_ranges(minoverlap = 2) %>%
      mutate(width = width + 1)
    message("02 - Longer fragments broken into sliding windows")
    FinalBins <- bind_ranges(BigFrags, SmallFrags)
  } else {
    FinalBins <- AllFrags
  }

  return(FinalBins)
}
