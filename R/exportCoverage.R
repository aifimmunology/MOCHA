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
        fileName <- gsub(" ","__", paste(x, groupColumn, names(cellPopSubsampleCov)[i], sep = "__"))
        if(!type){
          fileName = paste(fileName, "__Insertions", sep = '')
        }
        plyranges::write_bigwig(cellPopSubsampleCov[[i]],
                                paste(dir, "/", fileName, '.bw', sep =''))
    }
      
    if(saveFile){

      names(cellPopSubsampleCov) <- x
      GRangesList1 <- append(GRangesList1, cellPopSubsampleCov)
        
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


#' @title \code{exportDifferentials}
#'
#' @description \code{exportCoverage} will export normalized coverage files to BigWigs, either as sample-specific or sample-averaged files.
#'
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
exportDifferentials <- function(SampleTileObject, DifferentialsGRList, outDir, verbose = FALSE){
  
  genome <- BSgenome::getBSgenome(metadata(SampleTileObject)$Genome)
  
  for (i in seq_along(DifferentialsGRList)) {
    
    comparison_name <- names(DifferentialsGRList)[[i]]
    DiffPeaksGR <- DifferentialsGRList[[i]]
    
    # Set score and seqinfo for bigBed
    DiffPeaksGR$score <- 1
    seqinfo(DiffPeaksGR) <- seqinfo(genome)[seqnames(seqinfo(DiffPeaksGR))]
    
    outFile <- file.path(outDir, paste(comparison_name, sep="__"))
    if (verbose) { message("Exporting: ", outFile) } 
    
    # Output to bigbed
    rtracklayer::export.bb(DiffPeaksGR, paste0(outFile, ".bigBed"))
  }
  
}

#' @title \code{exportOpenTiles}
#'
#' @description \code{exportCoverage} will export normalized coverage files to BigWigs, either as sample-specific or sample-averaged files.
#'
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
exportOpenTiles <- function(SampleTileObject, outDir) {
  
  outDir <- './data/peaks_samplespecific/'
  genome <- BSgenome::getBSgenome(metadata(SampleTileObject)$Genome)
  
  for (cellPopulation in names(assays(SampleTileObject))) {
    cellPopMatrix <- MOCHA::getCellPopMatrix(
      SampleTileObject, 
      cellPopulation, 
      NAtoZero = FALSE
    )
    for (sample in colnames(SampleTileObject)){
      samplePeaks <- cellPopMatrix[, 
                                   colnames(cellPopMatrix) %in% c(sample),
                                   drop = FALSE]
      samplePeaksGR <- StringsToGRanges(rownames(samplePeaks))
      
      # Set score and seqinfo for bigBed
      samplePeaksGR$score <- 1
      seqinfo(samplePeaksGR) <- seqinfo(genome)[seqnames(seqinfo(samplePeaksGR))]
      
      sampleRow <- colData(SampleTileObject)[sample,]
      pbmc_sample_id <- sampleRow[['Sample']] # Enforced colname in callOpenTiles
      outFile <- file.path(outDir, paste(cellPopulation, pbmc_sample_id, sep="__"))
      if (verbose) { message("Exporting: ", outFile) }
      
      # Output to bigbed
      rtracklayer::export.bb(samplePeaksGR, paste0(outFile, ".bigBed"))
    }
  }
}

#' @title \code{exportMotifs}
#'
#' @description \code{exportCoverage} will export normalized coverage files to BigWigs, either as sample-specific or sample-averaged files.
#'
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
exportMotifs <- function(SampleTileObject, motifsGRanges, motifSetName, outDir, filterCellTypePeaks = FALSE, verbose = FALSE) {
  
  # Map over the seqinfo from our genome to the motifsGRanges
  # Required for bigBed export
  genome <- BSgenome::getBSgenome(metadata(SampleTileObject)$Genome)
  seqinfo(motifsGRanges) <- seqinfo(genome)[seqnames(seqinfo(motifsGRanges))]
  
  # # Truncate trailing zeroes from score https://www.biostars.org/p/235193/
  # # (Error : Trailing characters parsing integer in field 4 line 1 of text, got 10.3405262378482)
  motifsGRanges$score <- floor(motifsGRanges$score)
  
  # Filter out negative scores https://github.com/jhkorhonen/MOODS/issues/12
  motifsGRanges <- motifsGRanges[motifsGRanges$score > 0]
  
  
  if (filterCellTypePeaks) {
    # Split motifsGRangesFiltered by cell type peaks
    allPeaks <- rowRanges(SampleTileObject)
    samplePeakTable <- mcols(allPeaks)
    for (celltype in names(samplePeakTable)) {
      # Filter rows with boolean index
      samplePeaksGR <- allPeaks[samplePeakTable[[celltype]],]
      cellTypePeakMotifs <- plyranges::filter_by_overlaps(motifsGRanges, samplePeaksGR)
      
      outFile <- file.path(outDir, paste(celltype, motifSetName, "motifset.bigBed", sep="__"))
      if (verbose) { message("Exporting: ", outFile) }
      # Output to bigbed
      rtracklayer::export.bb(cellTypePeakMotifs, outFile)
    }
  } else {
    outFile <- file.path(outDir, paste(motifSetName, "motifset.bigBed", sep="__"))
    if (verbose) { message("Exporting: ", outFile) }
    # Output to bigbed
    rtracklayer::export.bb(motifsGRanges, outFile)
  }
}
