#' @title \code{exportCoverage}
#'
#' @description \code{exportCoverage} will export normalized coverage files to BigWigs, either as sample-specific or sample-averaged files.
#'
#' @param SampleTileObj The SummarizedExperiment object output from getSampleTileMatrix
#' @param dir string. Directory to save files to. 
#' @param cellPopulations vector of strings. Cell subsets for which to call
#'   peaks. This list of group names must be identical to names that appear in
#'   the SampleTileObj.  Optional, if cellPopulations='ALL', then peak
#'   calling is done on all cell populations. Default is 'ALL'.
#' @param type Boolean. Default is true, and exports Coverage. If set to FALSE, exports Insertions. 
#' @param groupColumn Optional, the column containing sample group labels for returning coverage within sample groups. Default is NULL, all samples will be used.
#' @param subGroups a list of subgroup(s) within the groupColumn from the metadata. Optional, default is NULL, all labels within groupColumn will be used.
#' @param sampleSpecific If TRUE, a BigWig will export for each sample-cell type combination.
#' @param saveFile Boolean. If TRUE, it will save to a BigWig. If FALSE, it will return the GRangesList. 
#' @param numCores integer. Number of cores to parallelize peak-calling across
#'   multiple cell populations
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#'
#' @return countSE a SummarizedExperiment containing coverage for the given input cell populations.
#'
#' @examples
#' \dontrun{
#'  MOCHA::exportBigWig(
#'   SampleTileObj = SampleTileMatrices,
#'   cellPopulations = "ALL",
#'   numCores = 30,
#'   sampleSpecific = FALSE
#' )
#' }
#' 
#' @export
#'

exportCoverage<- function(SampleTileObj,
                          dir = getwd(),
                          type = TRUE,
                          cellPopulations = "ALL",
                          groupColumn = NULL,
                          subGroups = NULL,
                          sampleSpecific = FALSE,
                          saveFile = TRUE,
                          numCores = 1,
                          verbose = FALSE) {
  . <- idx <- score <- NULL

  cellNames <- names(SummarizedExperiment::assays(SampleTileObj))
  metaFile <- SummarizedExperiment::colData(SampleTileObj)
  outDir <- SampleTileObj@metadata$Directory

  if (is.na(outDir)) {
    stop("Missing coverage file directory. SampleTileObj$metadata must contain 'Directory'.")
  }
  
  if (!file.exists(outDir)) {
    stop("Directory given by SampleTileObj@metadata$Directory does not exist.")
  }

  if (all(toupper(cellPopulations) == "ALL")) {
    cellPopulations <- cellNames
  }
  if (!all(cellPopulations %in% cellNames)) {
    stop("Some or all cell populations provided are not found.")
  }


  # Pull out a list of samples by group.

  if (!is.null(subGroups) & !is.null(groupColumn)) {
    # If the user defined a list of subgroup(s) within the groupColumn from the metadata, then it subsets to just those samples
    subSamples <- lapply(subGroups, function(x) metaFile[metaFile[, groupColumn] %in% x, "Sample"])
    names(subSamples) <- subGroups
  } else if (!is.null(groupColumn)) {

    # If no subGroup defined, then it'll form a list of samples across all labels within the groupColumn
    subGroups <- unique(metaFile[, groupColumn])

    subSamples <- lapply(subGroups, function(x) metaFile[metaFile[, groupColumn] %in% x, "Sample"])
  } else {

    # If neither groupColumn nor subGroup is defined, then it forms one list of all sample names
    subGroups <- "All"
    subSamples <- list("All" = metaFile[, "Sample"])
  }


  cl <- parallel::makeCluster(numCores)
  parallel::clusterEvalQ(cl, {
    library(GenomicRanges)
  })
  # Pull up the cell types of interest, and filter for samples and subset down to region of interest
  GRangesList1 <- NULL
  for(x in cellPopulations) {
    if(type){
      originalCovGRanges <- readRDS(paste(outDir, "/", x, "_CoverageFiles.RDS", sep = ""))
    }else{
      originalCovGRanges <- readRDS(paste(outDir, "/", x, "_InsertionFiles.RDS", sep = ""))
    }

    if (verbose) {
      message(stringr::str_interp("Extracting coverage from {x}."))
    }

    iterList <- lapply(subSamples, function(y){
        originalCovGRanges[y]
    })
    #If not sample specific, take the average coverage across samples. 
    # if it is sample specific, just subset down the coverage to the region of interest. 
    if (!sampleSpecific) {
        cellPopSubsampleCov <- pbapply::pblapply(cl = cl, X= iterList, averageCoverage)
        names(cellPopSubsampleCov) <- subGroups
    }else{
        names(iterList) <- subGroups
        cellPopSubsampleCov <- unlist(iterList, recursive = FALSE)
        names(cellPopSubsampleCov) <- gsub("\\.","__", names(cellPopSubsampleCov))
    }

    for(i in 1:length(cellPopSubsampleCov)){
        fileName <- gsub(" ","__", paste(groupColumn, names(cellPopSubsampleCov), sep = "__"))
        if(!type){
          fileName = paste(fileName, "__Insertions", sep = '')
        }
        plyranges::write_bigwig(cellPopSubsampleCov[[i]], paste(dir, "/", fileName, '.bw', sep =''))
    }
    if(saveFile){
      GRangesList1 <- append(GRangesList1, cellPopSubsampleCov)
      names(GRangesList1) <- c(names(GRangesList1), x)
    }
    
  }

  return(GRangesList1)
}


## helper function    
## Efficiently generates average single basepair coverage for a given region
averageCoverage <- function(coverageList){

    sampleCount <- length(coverageList)
    if(sampleCount > 1){
        mergedCounts <- IRanges::stack(methods::as(coverageList, "GRangesList"))
        mergedCounts <- plyranges::compute_coverage(mergedCounts, weight = mergedCounts$score / sampleCount)
    }else{
        mergedCounts = coverageList
    }                     


    return(mergedCounts)
}
