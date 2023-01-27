#' @title \code{testCoAccessibilityChromVar}
#'
#' @description \code{testCoAccessibilityChromVar} takes an input set of tile pairs and tests whether they are significantly different compared to a background set found via ChromVAR
#'
#' @param STObj The SummarizedExperiment object output from getSampleTileMatrix containing your sample-tile matrices
#' @param tile1 vector of indices or tile names (chrX:100-2000) for tile pairs to test (first tile in each pair)
#' @param tile2 vector of indices or tile names (chrX:100-2000) for tile pairs to test (second tile in each pair)
#' @param backNumber number of ChromVAR-matched background pairs. Default is 1000.
#' @param highMem Boolean to control memory usage. Default is FALSE. Only set highMem to TRUE if you have plenty of memory and want to run this function faster.s
#' @param numCores Optional, the number of cores to use with multiprocessing. Default is 1.
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @param ZI boolean flag that enables zero-inflated (ZI) Spearman correlations to be used. Default is TRUE. If FALSE, skip zero-inflation and calculate the normal Spearman.
#'
#' @return foreGround A data.frame with Tile1, Tile2, Correlation, and p-value for that correlation compared to the background
#'
#'
#' @export
testCoAccessibilityChromVar <- function(STObj,
                                        tile1,
                                        tile2,
                                        numCores = 1,
                                        ZI = TRUE,
                                        backNumber = 1000,
                                        highMem = FALSE,
                                        verbose = TRUE) {
  . <- NULL

  if (length(tile1) != length(tile2)) {
    stop("tile1 and tile2 must be the same length.")
  }

  fullObj <- getBackGroundObj(STObj)

  backPeaks <- chromVAR::getBackgroundPeaks(fullObj)

  if (is.character(tile1) && is.character(tile2)) {
    nTile1 <- match(tile1, rownames(fullObj))
    nTile2 <- match(tile1, rownames(fullObj))
  } else if (is.numeric(tile1) && is.numeric(tile2)) {
    nTile1 <- tile1
    nTile2 <- tile2

    tile1 <- rownames(fullObj)[nTile1]
    tile2 <- rownames(fullObj)[nTile2]
  } else {
    stop("tile1 and tile 2 must both be either numbers (indices) or strings")
  }

  if (backNumber > 2450) {
    backNumber <- 2450
    if (verbose) {
      warning("backNumber too high, setting to maximum of 2450.")
    }
  } else if (backNumber <= 10) {
    stop("backNumber too low (<=10). We recommend at least 100.")
  }

  cl <- parallel::makeCluster(numCores)
  parallel::clusterExport(cl, varlist = c("nTile1", "nTile2", "backPeaks"), envir = environment())

  if (verbose) {
    message("Finding background peak pairs")
  }

  backgroundCombos <- pbapply::pblapply(seq_along(nTile1), function(x) {
    tmpMat <- expand.grid(backPeaks[nTile1[x], ], backPeaks[nTile2[x], ])
    tmpMat <- tmpMat[tmpMat[, 1] != tmpMat[, 2], ]
    tmpMat[sample.int(dim(tmpMat)[1], backNumber), ]
  }, cl = cl)

  parallel::stopCluster(cl)

  gc()

  cl <- parallel::makeCluster(numCores)

  accMat <- SummarizedExperiment::assays(fullObj)[[1]]

  ## Test original pairs of locations
  combPairs <- data.frame(tile1, tile2)
  subAccMat <- accMat[unique(c(nTile1, nTile2)), ]

  if (verbose) {
    message("Identifying foreground")
  }

  parallel::clusterExport(cl, varlist = c("subAccMat", "combPairs"), envir = environment())
  foreGround <- runCoAccessibility(subAccMat, combPairs, ZI, verbose, cl)
  if (any(is.na(foreGround$Correlation))) {
    if (verbose) {
      warning("All foreground correlations are undefined")
    }
  }
  parallel::stopCluster(cl)

  gc()

  ## Now we need to test the background set

  if (verbose) {
    message("Identifying background correlations.")
  }

  if (highMem) {
    allBackCombos <- do.call("rbind", backgroundCombos)

    uniqueBackCombos <- unique(allBackCombos)

    combPairs <- data.frame(
      Tile1 = rownames(accMat)[uniqueBackCombos[, 1]],
      Tile2 = rownames(accMat)[uniqueBackCombos[, 2]]
    )

    subAccMat <- accMat[unique(c(uniqueBackCombos[, 1], uniqueBackCombos[, 2])), ]

    cl <- parallel::makeCluster(numCores)

    if (verbose) {
      message(
        paste("Generating Background correlations for all",
          length(tile1), "pairs",
          sep = " "
        )
      )
    }

    parallel::clusterExport(cl, varlist = c("subAccMat", "combPairs"), envir = environment())

    uniquebackGround <- runCoAccessibility(subAccMat, combPairs, ZI, verbose, cl)

    backGround <- dplyr::left_join(allBackCombos, uniquebackGround) %>%
      group_by(row_number(.) %/% backNumber) %>%
      group_map(~.x)

    parallel::stopCluster(cl)
  } else {
    ### low memory implementation

    backGround <- list()

    for (i in seq_along(backgroundCombos)) {
      combPairs <- data.frame(
        Tile1 = rownames(accMat)[backgroundCombos[[i]][, 1]],
        Tile2 = rownames(accMat)[backgroundCombos[[i]][, 2]]
      )

      subAccMat <- accMat[unique(c(backgroundCombos[[i]][, 1], backgroundCombos[[i]][, 2])), ]

      cl <- parallel::makeCluster(numCores)
      if (verbose) {
        message(paste("Generating Background correlations for", i, "of", length(tile1), "pairs", sep = " "))
      }
      parallel::clusterExport(cl, varlist = c("subAccMat", "combPairs"), envir = environment())

      tmp_background <- runCoAccessibility(subAccMat, combPairs, ZI, verbose, cl)

      parallel::stopCluster(cl)

      gc()

      # backGround <- append(backGround, tmp_background)
      backGround <- append(backGround, list(tmp_background))
    }
  }


  cl <- parallel::makeCluster(numCores)
  parallel::clusterExport(cl, varlist = c("foreGround", "backGround"), envir = environment())

  pValues <- pbapply::pblapply(seq_along(backGround), function(x) {
    cor1 <- foreGround$Correlation[x]
    
    if (is.na(cor1)){
      NA
    } else if (cor1 >= 0) {
      sum(cor1 > backGround[[x]]$Correlation) / length(backGround[[x]]$Correlation)
    } else if (cor1 < 0) {
      sum(cor1 < backGround[[x]]$Correlation) / length(backGround[[x]]$Correlation)
    }
  }, cl = cl) %>% unlist()
  parallel::stopCluster(cl)

  gc()

  foreGround$pValues <- pValues

  return(foreGround)
}

#' @title \code{testCoAccessibilityRandom}
#'
#' @description \code{testCoAccessibilityRandom} takes an input set of tile pairs and tests whether they are significantly different compared to random, non-overlapping background set.
#' @param STObj The SummarizedExperiment object output from getSampleTileMatrix containing your sample-tile matrices
#' @param tile1 vector of indices or tile names (chrX:100-2000) for tile pairs to test (first tile in each pair)
#' @param tile2 vector of indices or tile names (chrX:100-2000) for tile pairs to test (second tile in each pair)
#' @param backNumber number of ChromVAR-matched background pairs. Default is 1000.
#' @param numCores Optional, the number of cores to use with multiprocessing. Default is 1.
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @param ZI boolean flag that enables zero-inflated (ZI) Spearman correlations to be used. Default is TRUE. If FALSE, skip zero-inflation and calculate the normal Spearman.
#'
#' @return foreGround A data.frame with Tile1, Tile2, Correlation, and p-value for that correlation compared to the background
#'
#'
#' @export

testCoAccessibilityRandom <- function(STObj,
                                      tile1,
                                      tile2,
                                      numCores = 1,
                                      ZI = TRUE,
                                      backNumber = 1000,
                                      verbose = TRUE) {
  . <- NULL
  
  if (length(tile1) != length(tile2)) {
    stop("tile1 and tile2 must be the same length.")
  }

  fullObj <- getBackGroundObj(STObj)

  if (is.character(tile1) && is.character(tile2)) {
    nTile1 <- match(tile1, rownames(fullObj))
    nTile2 <- match(tile1, rownames(fullObj))
  } else if (is.numeric(tile1) && is.numeric(tile2)) {
    nTile1 <- tile1
    nTile2 <- tile2

    tile1 <- rownames(fullObj)[nTile1]
    tile2 <- rownames(fullObj)[nTile2]
  } else {
    stop("tile1 and tile 2 must both be either numbers (indices) or strings")
  }

  if (backNumber >= length(rownames(fullObj)) - c(length(tile1) + length(tile2))) {
    backNumber <- length(rownames(fullObj)) - c(length(tile1) + length(tile2))
    if (verbose) {
      warning("backNumber too high. Reset to all background combinations.")
    }
  } else if (backNumber <= 10) {
    stop("backNumber too low (<=10). We recommend 1000.")
  }

  cl <- parallel::makeCluster(numCores)

  accMat <- SummarizedExperiment::assays(fullObj)[[1]]

  ## Test original pairs of locations
  combPairs <- data.frame(tile1, tile2)
  subAccMat <- accMat[unique(c(nTile1, nTile2)), ]

  if (verbose) {
    message("Identifying foreground")
  }
  parallel::clusterExport(cl, varlist = c("subAccMat", "combPairs"), envir = environment())
  foreGround <- runCoAccessibility(subAccMat, combPairs, ZI, verbose, cl)
  
  if (any(is.na(foreGround$Correlation))) {
    if (verbose) {
      warning("All foreground correlations are undefined")
    }
  }
  
  parallel::stopCluster(cl)
  gc()

  if (verbose) {
    message("Finding background peak pairs")
  }

  backGroundTiles <- rownames(fullObj)[!rownames(fullObj) %in% c(tile1, tile2)]

  backgroundCombos <- data.frame(
    Tile1 = sample(backGroundTiles, backNumber),
    Tile2 = sample(backGroundTiles, backNumber)
  )
  backgroundCombos <- backgroundCombos[backgroundCombos[, 1] != backgroundCombos[, 2], ]

  ## Now we need to test the background set

  if (verbose) {
    message("Identifying background correlations.")
  }
  cl <- parallel::makeCluster(numCores)
  subAccMatB <- accMat[c(backgroundCombos$Tile1, backgroundCombos$nTile2), ]
  parallel::clusterExport(cl, varlist = c("backgroundCombos", "subAccMatB"), envir = environment())
  backGround <- runCoAccessibility(subAccMatB, backgroundCombos, ZI, verbose, cl)
  parallel::stopCluster(cl)

  cl <- parallel::makeCluster(numCores)
  parallel::clusterExport(cl, varlist = c("foreGround", "backGround"), envir = environment())

  pValues <- pbapply::pblapply(seq_along(tile1), function(x) {
    cor1 <- foreGround$Correlation[x]
    
    if (is.na(cor1)){
      NA
    } else if (cor1 >= 0) {
      sum(cor1 > backGround$Correlation) / length(backGround$Correlation)
    } else if (cor1 < 0) {
      sum(cor1 < backGround$Correlation) / length(backGround$Correlation)
    }
  }, cl = cl) %>% unlist()
  parallel::stopCluster(cl)

  gc()

  foreGround$pValues <- pValues

  return(foreGround)
}

#' @title \code{runCoAccessibility}
#'
#' @description \code{runCoAccessibility}
#'
#' @param STObj The SummarizedExperiment object output from getSampleTileMatrix containing your sample-tile matrices
#' @param accMat accessibility matrix to use for correlations
#' @param pairs data.frame for pairs of tiles to test. Must be tileNames (chrX:100-200).
#' @param numCores Optional, the number of cores to use with multiprocessing. Default is 1.
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @param ZI boolean flag that enables zero-inflated (ZI) Spearman correlations to be used. Default is TRUE. If FALSE, skip zero-inflation and calculate the normal Spearman.
#'
#' @return zi_spear_mat_tmp a data.table of tile pairs with associated correlations.
#'
#' @examples
#' runCoAccessibility(
#'   assays(SampleTileObj)[[1]],
#'   pairs = data.frame(
#'     Tile1 = c("chrX:1000-1499", "chr21:1000-1499"),
#'     Tile2 = c("chrX:1500-1999", "chr21:1500-1999")
#'   )
#' )
#' @noRd
#'
runCoAccessibility <- function(accMat, pairs, ZI = TRUE, verbose = TRUE, numCores = 1) {
  zero_inflated_spearman <- unlist(pbapply::pblapply(seq_len(dim(pairs)[1]),
    function(x) {
      tryCatch(
        {
          MOCHA:::weightedZISpearman(
            x = accMat[pairs[x, 1], ],
            y = accMat[pairs[x, 2], ],
            verbose = verbose,
            ZI = ZI
          )
        },
        error = function(e) {
          NA
        }
      )
    },
    cl = numCores
  ))

  # Create zero-inflated correlation matrix from correlation values
  zi_spear_mat_tmp <- data.table::data.table(
    Correlation = zero_inflated_spearman,
    Tile1 = pairs[, 1],
    Tile2 = pairs[, 2]
  )

  return(zi_spear_mat_tmp)
}


#' @title \code{getBackGroundObj}
#'
#' @description \code{getBackGroundObj} combines all celltypes in a SampleTileMatrix object into a SummarizedExperiment with one single matrix across all cell types and samples.
#'
#' @param STObj The SummarizedExperiment object output from getSampleTileMatrix containing your sample-tile matrices
#' @param NAToZero Set NA values in the sample-tile matrix to zero
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @return TileCorr A data.table correlation matrix
#'
#'
#' @export
getBackGroundObj <- function(STObj, NAtoZero = TRUE, verbose = FALSE) {

  Sample <- Freq <- . <- NULL
  # Extract all the Sample-Tile Matrices for each cell type
  assays <- SummarizedExperiment::assays(STObj)

  coldata <- SummarizedExperiment::colData(STObj)
  
  # Let's generate a new assay, that will contain the
  # the intensity for a given cell, as well as the
  # median intensity per sample-tile for all other cell types (i.e. the background)

  newAssays <- list(do.call("cbind", as(assays, "list")))
  newSamplesNames <- unlist(lapply(names(assays), function(x) {
    paste(x, colnames(STObj), sep = "__") %>% gsub(" ", "", .)
  }))

  names(newAssays) <- "counts"
  colnames(newAssays[[1]]) <- newSamplesNames

  if (NAtoZero) {
    newAssays[[1]][is.na(newAssays[[1]])] <- 0
  }


  allSampleData <- do.call("rbind", lapply(names(assays), function(x) {
    tmp_meta <- coldata
    tmp_meta$Sample <- paste(x, tmp_meta$Sample, sep = "__") %>% gsub(" ", "", .)
    tmp_meta$CellType <- rep(x, dim(tmp_meta)[1])
    rownames(tmp_meta) <- tmp_meta$Sample
    tmp_meta
  }))
  
  cellTypeLabelList <- Var1 <- NULL
  cellCounts <- as.data.frame(S4Vectors::metadata(STObj)$CellCounts) %>%
    dplyr::mutate(
      Sample = gsub(" ", "", paste(cellTypeLabelList, Var1, sep = "__"))) %>%
    dplyr::select(Sample, Freq)

  allSampleData <- dplyr::left_join(
    as.data.frame(allSampleData), cellCounts, by = "Sample")

  allRanges <- SummarizedExperiment::rowRanges(STObj)
  for (i in names(assays)) {
    GenomicRanges::mcols(allRanges)[, i] <- rep(TRUE, length(allRanges))
  }

  newObj <- SummarizedExperiment::SummarizedExperiment(
    assays = newAssays,
    colData = allSampleData,
    rowRanges = allRanges,
    metadata = S4Vectors::metadata(STObj)
  )
  
  tryCatch({
      newObj <- chromVAR::addGCBias(newObj, genome = S4Vectors::metadata(STObj)$Genome)
  }, error = function(e) {
    warning("The BSgenome for your organism is not installed")
    stop(e)
  })
  
  if (any(is.na(SummarizedExperiment::rowData(newObj)$bias))) {
    naList <- is.na(SummarizedExperiment::rowData(newObj)$bias)
    
    if (verbose) {
      message(paste(sum(naList), "NaNs found within GC Bias", sep = " "))
    }
    
    SummarizedExperiment::rowData(newObj)$bias[which(naList)] <- mean(rowData(newObj)$bias, na.rm = TRUE)
  }

  return(newObj)
}
