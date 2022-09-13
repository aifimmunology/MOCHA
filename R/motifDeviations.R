

### ChromVar-inspired analysis for how specific peaksets use less or more
### Designed to be used with scMACs normalization protocol. Use with other normalization schemes may not control for technical bias.



######## getMotifCoverage: Function that takes in Sample-specific coverage data and finds
## the average coverage for each motif location within a list of motif locations
## @covFiles: GRangesList of coverage for each sample
## @MotifLocations: GRangesList of Motif locations
## @numCores: number of cores to parallelize over.

getMotifCoverages <- function(covFiles, MotifLocations, numCores = 1) {
  allMotifs <- stack(MotifLocations)
  counts <- mclapply(covFiles, function(x) {
    x %>%
      plyranges::mutate(NewScore = score) %>%
      plyranges::join_overlap_intersect(allMotifs) %>%
      plyranges::mutate(WeightedScore = NewScore * width(.)) %>%
      plyranges::group_by(name) %>%
      plyranges::reduce_ranges(score = mean(WeightedScore))
  }, mc.cores = numCores)
  return(counts)
}


######## matchPeaks: Function that takes in Sample-Tile matrix and a foreground set of peaks
## and finds peakset that matches the intensity within the rest of the peaks (or a predefined background set).
## @covFiles: GRangesList of coverage for each sample
## @MotifLocations: GRangesList of Motif locations
## @numCores: number of cores to parallelize over.


matchPeaks <- function(SampleTileMatrix, ForegroundPeakset, BackgroundSet = NULL, method = "Mean", N = 50, numCores = 1) {
  if (class(ForegroundPeakset)[1] == "Character") {
    ForegroundPeakSet <- MOCHA::StringsToGRanges(ForegroundPeakset)
  } else if (class(ForegroundPeakset)[1] != "GRanges") {
    stop("Error: Foreground Peakset must either be a GRanges object or a vector of string in the format Chr1:100-200")
  }

  if (is.null(BackgroundSet)) {
    BackgroundSet <- plyranges::filter_by_non_overlaps(MOCHA::StringsToGRanges(SampleTileMatrix$tileID), ForegroundPeakset)
  } else if (!is.null(BackgroundSet) & class(BackgroundSet)[1] == "Character") {
    BackgroundSet <- MOCHA::StringsToGranges(BackgrounfSet)
  } else if (!is.null(BackgroundSet) & class(BackgroundSet)[1] != "GRanges") {
    stop("Error: BackgroundSet must either be a GRanges object, a vector of string in the format Chr1:100-200, or set to NULL")
  }

  Count1 <- findOverlaps(MOCHA::StringsToGRanges(SampleTileMatrix$tileID), ForegroundPeakset)
  Count2 <- findOverlaps(MOCHA::StringsToGRanges(SampleTileMatrix$tileID), BackgroundSet)

  # Check if the background and foreground set overlap with any tiles in common.
  if (any(queryHits(Count1) %in% queryHits(Count2))) {
    stop("Error: Background & Foreground set overlap for Sample-Tile matrix")
  }

  ForeMat <- SampleTileMatrix[queryHits(Count1), -"tileID"]

  BackMat <- SampleTileMatrix[queryHits(Count2), -"tileID"]

  if (toLower(Method) == "mean") {
    avgForeMat <- rowMeans(as.matrix(ForeMat))
    avgBackMat <- rowMeans(as.matrix(BackMat))
  } else if (toLower(Method) == "median") {
    avgForeMat <- rowMedians(as.matrix(ForeMat))
    avgBackMat <- rowMedians(as.matrix(BackMat))
  } else if (toLower(Method) == "sd") {
    avgForeMat <- rowSds(as.matrix(ForeMat))
    avgBackMat <- rowMSds(as.matrix(BackMat))
  } else {
    stop("Error: Method not recognize")
  }


  matchedPeak <- mclapply(avgForeMat, function(x) {
    order(abs(x - avgBackMat))[c(1:N)]
  }, mc.cores = numCores)

  if (returnGR) {
    matchedPeakReg <- BackgroundSet[sort(unlist(unique(matchedPeak)))]
    return(matchedPeakReg)
  } else {
    subMat <- BackMat[sort(unlist(unique(matchedPeak))), ]
    return(subMat)
  }
}



normalized_wilcoxon <- function(data, group, point.mass = 0, test = "wilcoxon") {
  # Function for calculating two-part statistics
  Index1 <- c(group == 1)
  Group1 <- data[Index1]
  Group0 <- data[!Index1]
  n1 <- length(Group1)
  n2 <- length(Group0)
  obs <- c(n1, n2)
  success <- c(sum(Group1 != point.mass), sum(Group0 != point.mass))
  pointmass <- obs - success
  if (sum(success) == 0) {
    T2 <- 0
    B2 <- 0
  } else if ((success[1] == 0) | (success[2] == 0)) {
    T2 <- 0
    B2 <- prop.test(pointmass, obs)$statistic
  } else if ((success[1] == 1) | (success[2] == 1)) {
    T2 <- 0
    B2 <- prop.test(pointmass, obs)$statistic
  } else {
    uniq1 <- length(unique(Group1[Group1 != point.mass]))
    uniq2 <- length(unique(Group0[Group0 != point.mass]))
    if ((uniq1 < 2) & (uniq2 < 2)) {
      T2 <- 0
      if (sum(pointmass) == 0) {
        B2 <- 0
      } else {
        B2 <- prop.test(pointmass, obs)$statistic
      }
    } else if (sum(pointmass) == 0) {
      B2 <- 0
      if (test == "t.test") {
        T2 <- t.test(data ~ group)$statistic
      }
      if (test == "wilcoxon") {
        W <- wilcox.test(data ~ group, exact = FALSE)$statistic
        mu <- (n1 * n2) / 2
        sigma <- sqrt((n1 * n2 * (n1 + n2 + 1)) / 12)
        T2 <- ((W - mu) - 0.5) / sigma
      }
    } else {
      B2 <- prop.test(pointmass, obs)$statistic
      contIndex <- data != point.mass
      cont <- data[contIndex]
      cGroup <- group[contIndex]
      n1c <- sum(cGroup == 1)
      n2c <- sum(cGroup == 0)
      if (test == "t.test") {
        T2 <- t.test(cont ~ cGroup)$statistic
      }
      if (test == "wilcoxon") {
        W <- wilcox.test(cont ~ cGroup, exact = FALSE)$statistic
        mu <- (n1c * n2c) / 2
        sigma <- sqrt((n1c * n2c * (n1c + n2c + 1)) / 12)
        T2 <- ((W - mu) - 0.5) / sigma
      }
    }
  }
  return(T2)
}



##### getMotifDev

getMotifDev <- function(MotifCoverage, ForegroundPeaks, BackgroundPeaks,
                        numCores = 1, testType = "Wilcoxon", verbose = FALSE) {
  ForeSubset <- mclapply(MotifCoverage, function(x) plyranges::filter_by_overlaps(x, ForegroundPeaks), mc.cores = numCores)
  BackSubset <- mclapply(MotifCoverage, function(x) plyranges::filter_by_overlaps(x, BackgroundPeaks), mc.cores = numCores)


  ForeDF <- mclapply(seq_along(MotifCoverage), function(x) {
    tmp1 <- data.frame(GRangesToString(ForeDF[[x]]), ForeDF[[x]]$name, ForeDF[[x]]$score)
    colnames(tmp1) <- c("Peak", "Motif", names(ForeSubset)[x])
    tmp1
  }, mc.cores = 25) %>% purrr::reduce(dplyr::full_join, by = c("Peak", "Motif"))

  ForeDF[is.na(ResDFList)] <- 0

  ResDF <- mclapply(seq_along(MotifCoverage), function(x) {
    tmp1 <- data.frame(GRangesToString(BackSubset[[x]]), BackSubset[[x]]$name, BackSubset[[x]]$score)
    colnames(tmp1) <- c("Peak", "Motif", names(BackSubset)[x])
    tmp1
  }, mc.cores = 25) %>% purrr::reduce(dplyr::full_join, by = c("Peak", "Motif"))


  ResDF[is.na(ResDF)] <- 0

  motifList <- names(MotifLocations)

  motifDev <- mclapply(seq_along(motifList), function(x) {
    timeSub <- TimeDFList %>% dplyr::filter(Motif == motifList[x])
    resSub <- ResDFList %>% dplyr::filter(Motif == motifList[x])
    ## Iterate over all samples, comparing each sample to the average.
    motif1Dev <- lapply(seq_along(dim(timeSub)[2] - 2), function(y) {
      comb -> rbind(timeSub[, y + 2], resSub[, y + 2])
      group1 <- c(rep(1, dim(timeSub)[1]), rep(0, dim(resSub)[1]))

      normalized_wilcoxon(comb, group1, test = testType)
    }) %>% cbind()

    colnames(motif1Dev) <- colnames(timeSub)[, -c(1:2)]

    motif1Dev
  }, mc.cores = numCores) %>% do.call("rbind", .)

  rownames(motifDev) <- motifList

  return(motifDev)
}
