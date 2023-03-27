
### getCoverage turns sample/celltype-specific fragments lists into
### sample-specific coverage files for each sample-celltype.
# @popFrags - GRangesList of fragments for all sample/celltypes
# @normFactor - Normalizaiton factor. Can be either be one, in which case all coverage files will be normalized by the same value, or the same length as the GRangesList
# @filterEmpty - True/False flag on whether or not to carry forward regions without coverage.
# @numCores - number of cores to parallelize over
# @verbose - Boolean variable to determine verbosity of output.

getCoverage <- function(popFrags, normFactor, TxDb, cl, filterEmpty = FALSE, verbose = FALSE) {
  score <- NULL
  if (length(normFactor) == 1) {
    normFactor <- rep(normFactor, length(popFrags))
  } else if (length(normFactor) != length(popFrags)) {
    stop("Length of normFactor is equal to length of popFrags. Please either give 1 value, or a vector of equal length to popFrags.")
  }

  if (verbose) {
      message(paste("Counting", paste0(names(popFrags), collapse = ", "), sep = " "))
  }

  popFragList <- lapply(seq_along(popFrags), function(x){

    list(popFrags[[x]], normFactor[[x]], filterEmpty)

  })

  # Summarize the coverage over the region window at a single basepair resolution
  popCounts <- pbapply::pblapply(popFragList, calculateCoverage, cl = cl)

  popCounts <- lapply(popCounts, function(x){
          GenomeInfoDb::seqinfo(x) <- GenomeInfoDb::seqinfo(TxDb)[GenomicRanges::seqnames(GenomeInfoDb::seqinfo(x))]
          x
  })
    
  names(popCounts) <- names(popFrags)

  return(popCounts)
}

######## calculateCoverage: Function that takes in a GRanges fragment object and generates coverage GRanges.
## @param ref
calculateCoverage <- function(ref){

  popFrags <- ref[[1]]
  Num <- ref[[2]]
  filterEmpty <- ref[[3]]

  counts_gr <- plyranges::compute_coverage(popFrags) %>% plyranges::mutate(score = score / Num)

  if (filterEmpty) {
    plyranges::filter(counts_gr, score > 0)
  } else {
    counts_gr
  }

}



######## getSpecificCoverage: Function that takes in a GRangesList of coverage per sample/celltype and finds
## the average coverage intensity for just specific regions.
## @covFiles: GRangesList of coverage for each sample
## @regions: Regions to count intensities over, must be non-overlapping and non-adjacent ( > 1 bp apart).
## @numCores: number of cores to parallelize over.
getSpecificCoverage <- function(covFiles, regions, numCores = 1) {
  score <- NewScore <- WeightedScore <- . <- NULL
  counts <- parallel::mclapply(covFiles, function(x) {
    x %>%
      plyranges::mutate(NewScore = score) %>%
      plyranges::join_overlap_intersect(regions) %>%
      plyranges::mutate(WeightedScore = NewScore * GenomicRanges::width(.)) %>%
      plyranges::reduce_ranges(score = mean(WeightedScore))
  }, mc.cores = numCores)

  return(counts)
}
