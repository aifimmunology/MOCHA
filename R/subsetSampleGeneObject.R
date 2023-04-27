#' @title \code{subsetSampleGeneObject}
#'
#' @description \code{subsetSampleGeneObject} subsets SummarizedExperiment-type object that contains pseudobulk RNA expression, either by cell type or sample metadata.
#'
#' @param Object A SummarizedExperiment object
#' @param subsetBy the variable to subset by. Can either be 'celltype', or a
#'   column from the sample metadata (see colData(Object))
#' @param groupList the list of cell type names or sample-associated data that
#'   should be used to subset the Object
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#'
#' @return Object the input Object, filtered down to either the cell type or
#'   samples desired.
#'
#'
#' @export

subsetSampleGeneObject <- function(Object,
                              subsetBy,
                              groupList,
                              verbose = FALSE) {

 if (class(Object)[1] == "SummarizedExperiment") {
    sampleData <- SummarizedExperiment::colData(Object)

    # To subset by cell type, first we have to verify that all cell type names were found within the  object.
    # then we simply do a simple subsetting process, like you would with a list.
    if (grepl('celltype', tolower(subsetBy))) {
      if (!all(groupList %in% names(SummarizedExperiment::assays(Object)))) {
        stop("Error: groupList includes celltypes not found within Object.")
      }

      keep <- which(names(SummarizedExperiment::assays(Object)) %in% groupList)
      SummarizedExperiment::assays(Object) <- SummarizedExperiment::assays(Object)[keep]

      return(Object)
    }
  }


  if (subsetBy %in% colnames(sampleData)) {
    if (!all(groupList %in% unique(sampleData[[subsetBy]]))) {
      stop("Error: groupList includes names not found within the object sample data. Please check groupList.")
    }
  }
  keep <- rownames(sampleData)[which(sampleData[[subsetBy]] %in% groupList | is.na(sampleData[[subsetBy]]))]
  return(Object[, keep])
}

#' @title \code{combineSampleGeneMatrix}
#'
#' @description \code{combineSampleGeneMatrix} combines all celltypes in a
#'   SampleGeneMatrix object into a SummarizedExperiment with one single matrix
#'   across all cell types and samples, annotating GC bias using
#'   chromVAR.
#'
#' @param SampleGeneMatrix The SummarizedExperiment object output from
#'   getSampleGeneMatrix containing your sample-Gene matrices
#' @param NAToZero Set NA values in the sample-Gene matrix to zero
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @return SummarizedExperiment. Contains one assays for all cell types.
#'
#'
#' @export
combineSampleTileMatrix <- function(SampleTileMatrix,
                                    verbose = FALSE) {
  assays <- SummarizedExperiment::assays(SampleTileMatrix)

  coldata <- SummarizedExperiment::colData(SampleTileMatrix)
  
  # Let's generate a new assay, that will contain the
  # the intensity for a given cell, as well as the
  # median intensity per sample-tile for all other cell types (i.e. the background)

  newAssays <- list(do.call("cbind", as(assays, "list")))
  newSamplesNames <- unlist(lapply(names(assays), function(x) {
    paste(x, colnames(SampleTileMatrix), sep = "__") %>% gsub(" ", "", .)
  }))

  names(newAssays) <- "counts"
  colnames(newAssays[[1]]) <- newSamplesNames

  allSampleData <- do.call("rbind", lapply(names(assays), function(x) {
    tmp_meta <- coldata
    tmp_meta$Sample <- paste(x, tmp_meta$Sample, sep = "__") %>% gsub(" ", "_", .)
    tmp_meta$CellType <- rep(x, dim(tmp_meta)[1])
    rownames(tmp_meta) <- tmp_meta$Sample
    tmp_meta
  }))
  
  cellTypeLabelList <- Var1 <- NULL
  cellCounts <- as.data.frame(S4Vectors::metadata(SampleTileMatrix)$CellCounts) %>%
    dplyr::mutate(
      Sample = gsub(" ", "_", paste(cellTypeLabelList, Var1, sep = "__"))) %>%
    dplyr::select(Sample, Freq)

  allSampleData <- dplyr::left_join(
    as.data.frame(allSampleData), cellCounts, by = "Sample")

  newObj <- SummarizedExperiment::SummarizedExperiment(
    assays = newAssays,
    colData = allSampleData,
    metadata = S4Vectors::metadata(SampleTileMatrix)
  )
  
  return(newObj)
}