
#' @title \code{combinePseudobulkSE}
#'
#' @description \code{combinePseuduoulkSE} combines all celltypes in a
#'   PseudobulkSE objects (SampleTileMatrix or pseodobulked scRNA) into a SummarizedExperiment with one single matrix
#'   across all cell types and samples, annotating GC bias using
#'   chromVAR.
#'
#' @param SummarizedExp The SummarizedExperiment object output from
#'   getSampleTileMatrix or makePseudobulkSE contained cell-type-specific matrices.
#' @param NAToZero Set NA values in the sample-tile matrix to zero
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @return TileCorr A data.table correlation matrix
#'
#'
#' @export
combinePseudobulkSE <- function(SummarizedExp,
                                    NAtoZero = TRUE, 
                                    verbose = FALSE) {
  
  genome <- S4Vectors::metadata(SummarizedExp)$Genome
  genome <- BSgenome::getBSgenome(genome)

  Sample <- Freq <- . <- NULL
  # Extract all the Sample-Tile Matrices for each cell type
  assays <- SummarizedExperiment::assays(SummarizedExp)

  coldata <- SummarizedExperiment::colData(SummarizedExp)
  
  # Let's generate a new assay, that will contain the
  # the intensity for a given cell, as well as the
  # median intensity per sample-tile for all other cell types (i.e. the background)

  newAssays <- list(do.call("cbind", as(assays, "list")))
  newSamplesNames <- unlist(lapply(names(assays), function(x) {
    gsub(" ", "_", paste(x, colnames(SummarizedExp), sep = "__"))
  }))

  names(newAssays) <- "counts"
  colnames(newAssays[[1]]) <- newSamplesNames

  if (NAtoZero) {
    newAssays[[1]][is.na(newAssays[[1]])] <- 0
  }

  allSampleData <- do.call("rbind", lapply(names(assays), function(x) {
    tmp_meta <- coldata
    tmp_meta$Sample <- gsub(" ", "_", paste(x, tmp_meta$Sample, sep = "__"))
    tmp_meta$CellType <- rep(x, dim(tmp_meta)[1])
    rownames(tmp_meta) <- tmp_meta$Sample
    tmp_meta
  }))
  
  if(all(names(S4Vectors::metadata(SummarizedExp)) %in% FragmentCounts)){
  
    cellTypeLabelList <- Var1 <- NULL
    cellCounts <- as.data.frame(S4Vectors::metadata(SummarizedExp)$CellCounts)
    cellCounts <- dplyr::mutate(cellCounts,
        Sample = gsub(" ", "_", paste(cellTypeLabelList, Var1, sep = "__")))
    cellCounts <- dplyr::select(cellCounts, Sample, Freq)

    fragCounts <- as.data.frame(S4Vectors::metadata(SummarizedExp)$FragmentCounts)
    fragCounts <-  dplyr::select(fragCounts, dplyr::one_of(names(assays))) 
    fragCounts <-  dplyr::mutate(fragCounts, Sample = rownames(fragCounts))
    fragCounts <-  tidyr::pivot_longer(fragCounts, cols = names(assays), names_to = 'CellTypes', values_to = 'FragNumber') 
    fragCounts <-  dplyr::mutate(fragCounts, Sample = gsub(" ", "_", paste(CellTypes, Sample, sep = "__")))
    fragCounts <-  dplyr::select(fragCounts, Sample, FragNumber)
  }else{

    cellCounts <- as.data.frame(S4Vectors::metadata(SummarizedExp)$CellCounts)
    cellCounts <- dplyr::mutate(cellCounts,
        Sample = gsub(" ", "_", paste(cellTypeLabelList, Var1, sep = "__")))
    cellCounts <- dplyr::select(cellCounts, Sample, Freq)
    colnames(cellCounts) <- c('Sample','CellCounts')

  }

  allSampleData <- dplyr::left_join(
    as.data.frame(allSampleData), cellCounts, by = "Sample")

  allSampleData <- dplyr::left_join(
    allSampleData,  fragCounts, by = 'Sample'
  )

  allRanges <- SummarizedExperiment::rowRanges(SummarizedExp)
  for (i in names(assays)) {
    GenomicRanges::mcols(allRanges)[, i] <- rep(TRUE, length(allRanges))
  }

  newObj <- SummarizedExperiment::SummarizedExperiment(
      assays = newAssays,
      colData = allSampleData,
      rowRanges = allRanges,
      metadata = S4Vectors::metadata(SummarizedExp)
    )

  if(all(names(S4Vectors::metadata(SummarizedExp)) %in% FragmentCounts)){

    
    newObj <- chromVAR::addGCBias(newObj, genome = genome)
    
    if (any(is.na(SummarizedExperiment::rowData(newObj)$bias))) {
      naList <- is.na(SummarizedExperiment::rowData(newObj)$bias)
      
      if (verbose) {
        warning(paste(sum(naList), "NaNs found within GC Bias", sep = " "))
      }
      
      SummarizedExperiment::rowData(newObj)$bias[which(naList)] <- mean(rowData(newObj)$bias, na.rm = TRUE)
    }
  }else{
  
  }

  return(newObj)
}