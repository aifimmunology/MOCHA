#' @title \code{subsetObject}
#'
#' @description \code{subsetObject} subsets a tileResults-type object (from
#'   callOpenTiles), or a SummarizedExperiment-type object (from
#'   getSampleTileMatrix), either by cell type or sample metadata.
#'
#' @param Object A MultiAssayExperiment or RangedSummarizedExperiment,
#' @param subsetBy The variable to subset by. Can either be 'celltype', or a
#'   column from the sample metadata (see `colData(Object)`).
#' @param groupList the list of cell type names or sample-associated data that
#'   should be used to subset the Object
#' @param removeNA If TRUE, removes groups in groupList that are NA. If FALSE, 
#'   keep groups that are NA. 
#' @param subsetPeaks If `subsetBy` = 'celltype', subset the tile 
#'   set to tiles only called in those cell types. Default is TRUE.
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#'
#' @return Object the input Object, filtered down to either the cell type or
#'   samples desired.
#'
#'
#' @export
subsetMOCHAObject <- function(Object,
                              subsetBy,
                              groupList,
                              removeNA = TRUE,
                              subsetPeaks = TRUE,
                              verbose = FALSE) {
  
  if (!(subsetBy %in% colnames(sampleData)) & !grepl('celltype', tolower(subsetBy))) {
    stop(
      "Variable given in subsetBy is not in the colData of the input Object.",
      "subsetBy must either be 'celltype', or a column name within colData(Object)."
    )
  }
  
  # Subset cell populations (assays)
  if (grepl('celltype', tolower(subsetBy))) {
    
    if (class(Object)[1] == "MultiAssayExperiment") {
      # Input is tileResults, output of callOpenTiles
      
      sampleData <- MultiAssayExperiment::colData(Object)
    
      if ((subsetBy %in% colnames(sampleData)) & grepl('celltype', tolower(subsetBy))) {
        if (verbose) {
          warning(
            "subsetBy is set to 'celltype', but that is also a column name ",
            "within the colData of the input Object. The object will be filtered",
            " by cell type annotation, not by colData of the input Object.")
        }
      }
    
      if (grepl('celltype', tolower(subsetBy))) {
        if (!all(groupList %in% names(Object))) {
          stop("groupList includes celltypes not found within Object.")
        }
    
        keep <- MultiAssayExperiment::subsetByAssay(Object, groupList)
        
        # TODO SUBSET SUMMARIZEDDATA
        return(keep)
      }
    
    } else if (class(Object)[1] == "RangedSummarizedExperiment") {
      # Input is a TSAM, output of getSampleTileMatrix
      
      sampleData <- SummarizedExperiment::colData(Object)
      
      # To subset by cell type, first we have to verify that all cell type 
      # names were found within the  object.
      # then we simply do a simple subsetting process, like you would with a list.
      if (!all(groupList %in% names(SummarizedExperiment::assays(Object)))) {
        stop("groupList includes celltypes not found within Object.")
      }
      
      keep <- which(names(SummarizedExperiment::assays(Object)) %in% groupList)
      SummarizedExperiment::assays(Object) <- SummarizedExperiment::assays(Object)[keep]
      
      # Subset peaks
      if(subsetPeaks){
      
        rowMeta <- GenomicRanges::mcols(SummarizedExperiment::rowRanges(Object))[,groupList]
      
        if(!is.null(dim(rowMeta))){
      
          rowMeta = rowSums(as.data.frame(rowMeta)) > 0
      
        }
        Object <- Object[rowMeta,]
      }
      
      # TODO SUBSET SUMMARIZEDDATA
      return(Object)
    }
  }

  # Subset samples by sample metadata
  if (subsetBy %in% colnames(sampleData)) {
    if (!all(groupList %in% unique(sampleData[[subsetBy]]))) {
      stop(str_interp(
        "groupList includes names not found within the column '${subsetBy}'"
        ),
        " in the sample metadata. (see `colData(Object)`). ")
    }
    
    if (removeNA) {
      keep <- rownames(sampleData)[which(sampleData[[subsetBy]] %in% groupList)]
    } else {
      keep <- rownames(sampleData)[which(
        sampleData[[subsetBy]] %in% groupList | is.na(sampleData[[subsetBy]])
      )]
    }
    
    # TODO SUBSET SUMMARIZEDDATA
    return(Object[, keep])
  }

  
}
