#' @title \code{subsetObject}
#'
#' @description \code{subsetObject} subsets a tileResults-type object (from
#'   callOpenTiles), or a SummarizedExperiment-type object (from
#'   getSampleTileMatrix), either by cell type or sample metadata.
#'
#' @param Object A MultiAssayExperiment or RangedSummarizedExperiment,
#' @param subsetBy the variable to subset by. Can either be 'celltype', or a
#'   column from the sample metadata (see colData(Object))
#' @param groupList the list of cell type names or sample-associated data that
#'   should be used to subset the Object
#' @param na.rm removes groups that are NA if set to true. If set to false, then
#'   you filter for everything in the groupList and also NA values.
#'
#' @return Object the input Object, filtered down to either the cell type or
#'   samples desired.
#'
#'
#' @export

subsetMOCHAObject <- function(Object,
                              subsetBy,
                              groupList,
                              na.rm = TRUE) {
  if (class(Object)[1] == "MultiAssayExperiment") {
    sampleData <- MultiAssayExperiment::colData(Object)

    if (!(subsetBy %in% colnames(sampleData)) & tolower(subsetBy) != "celltypes") {
      stop("Error: subsetBy must either be a column name within colData(Objects), or 'celltype'.")
    }

    if ((subsetBy %in% colnames(sampleData)) & tolower(subsetBy) == "celltypes") {
      warning("subsetBy is set to CellTypes, but that is also a column name within the Sample metadata. The object will be filtered by cell type annotation, not by sample metadata.")
    }

    # To subset by cell type, first we have to verify that all cell type names were found within the  object.
    # then we simply do a simple subsetting process, like you would with a list.
    if (tolower(subsetBy) == "celltypes") {
      if (!all(groupList %in% names(Object))) {
        stop("Error: groupList includes celltypes not found within Object.")
      }

      keep <- which(names(Object) %in% groupList)
      return(Object[keep])
    }
  } else if (class(Object)[1] == "RangedSummarizedExperiment") {
    sampleData <- SummarizedExperiment::colData(Object)

    # To subset by cell type, first we have to verify that all cell type names were found within the  object.
    # then we simply do a simple subsetting process, like you would with a list.
    if (tolower(subsetBy) == "celltypes") {
      if (!all(groupList %in% names(Object))) {
        stop("Error: groupList includes celltypes not found within Object.")
      }

      keep <- which(names(SummarizedExperiment::assays(Object)) %in% groupList)

      SummarizedExperiment::assays(object) <- SummarizedExperiment::assays(object)[keep]

      return(Object)
    }
  }


  if (subsetBy %in% colnames(sampleData)) {
    if (!all(groupList %in% unique(sampleData[[subsetBy]]))) {
      stop("Error: groupList includes names not found within the object sample data. Please check groupList.")
    }
  }

  if (na.rm) {
    keep <- rownames(sampleData)[which(sampleData[[subsetBy]] %in% groupList)]

    return(Object[, keep])
  } else {
    keep <- rownames(sampleData)[which(sampleData[[subsetBy]] %in% groupList | is.na(sampleData[[subsetBy]]))]
    return(Object[, keep])
  }
}
