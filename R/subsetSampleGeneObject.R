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