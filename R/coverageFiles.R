
### getCoverage turns sample/celltype-specific fragments lists into
### sample-specific coverage files for each sample-celltype.
# @popFrags - GRangesList of fragments for all sample/celltypes
# @normFactor - Normalizaiton factor. Can be either be one, in which case all coverage files will be normalized by the same value, or the same length as the GRangesList
# @filterEmpty - True/False flag on whether or not to carry forward regions without coverage.
# @numCores - number of cores to parallelize over
# @verbose - Boolean variable to determine verbosity of output.

getCoverage <- function(popFrags, normFactor, TxDb, filterEmpty = FALSE, numCores = 1, verbose = FALSE) {
  score <- NULL
  if (length(normFactor) == 1) {
    normFactor <- rep(normFactor, length(popFrags))
  } else if (length(normFactor) != length(popFrags)) {
    stop("Length of normFactor is equal to length of popFrags. Please either give 1 value, or a vector of equal length to popFrags.")
  }

  cl <- parallel::makeCluster(numCores)
  parallel::clusterExport(
      cl=cl,
      varlist=c("popFrags", "normFactor", "frags", "TxDb", "filterEmpty", "verbose"),
      envir=environment()
  )

  # Summarize the coverage over the region window at a single basepair resolution
  popCounts <- pbapply::pblapply(cl = cl, 
              X = seq_along(popFrags), 
              FUN = function(x) {
    if (verbose) {
      message(paste("Generate coverage files for", names(popFrags)[x], sep = " "))
    }

    Num <- normFactor[[x]]
    counts_gr <- plyranges::compute_coverage(popFrags[[x]]) %>% plyranges::mutate(score = score / Num)

    GenomeInfoDb::seqinfo(counts_gr) <- GenomeInfoDb::seqinfo(TxDb)[GenomicRanges::seqnames(GenomeInfoDb::seqinfo(counts_gr))]

    if (filterEmpty) {
      plyranges::filter(counts_gr, score > 0)
    } else {
      counts_gr
    }
  })

  parallel::stopCluster(cl)
  gc()

  names(popCounts) <- names(popFrags)

  return(popCounts)
}



######## getSpecificCoverage: Function that takes in a GRangesList of coverage per sample/celltype and finds
## the average coverage intensity for just specific regions.
## @covFiles: GRangesList of coverage for each sample
## @regions: Regions to count intensities over, must be non-overlapping and non-adjacent ( > 1 bp apart).
## @numCores: number of cores to parallelize over.
getSpecificCoverage <- function(covFiles, regions, numCores = 1) {
  
  score <- NewScore <- WeightedScore <- . <- NULL
  
  cl <- parallel::makeCluster(numCores)
    parallel::clusterExport(
        cl=cl,
        varlist=c("regions"),
        envir=environment()
    )

  counts <- pbapply::pblapply(
    cl = cl,
    X = covFiles, 
    FUN = function(x) {
    x %>%
      plyranges::mutate(NewScore = score) %>%
      plyranges::join_overlap_intersect(regions) %>%
      plyranges::mutate(WeightedScore = NewScore * GenomicRanges::width(.)) %>%
      plyranges::reduce_ranges(score = mean(WeightedScore))
  })

  parallel::stopCluster(cl)

  return(counts)
}
