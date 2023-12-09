#' @title \code{exportCoverage}
#'
#' @description \code{exportCoverage} will export normalized coverage files to
#'   BigWig files, either as sample-specific or sample-averaged files, for
#'   visualization in genome browsers.
#'
#' @param SampleTileObj The SummarizedExperiment object output from
#'   getSampleTileMatrix
#' @param dir string. Directory to save files to.
#' @param cellPopulations vector of strings. Cell subsets for which to call
#'   peaks. This list of group names must be identical to names that appear in
#'   the SampleTileObj.  Optional, if cellPopulations='ALL', then peak calling
#'   is done on all cell populations. Default is 'ALL'.
#' @param type Boolean. Default is true, and exports Coverage. If set to FALSE,
#'   exports Insertions.
#' @param groupColumn Optional, the column containing sample group labels for
#'   returning coverage within sample groups. Default is NULL, all samples will
#'   be used.
#' @param subGroups a list of subgroup(s) within the groupColumn from the
#'   metadata. Optional, default is NULL, all labels within groupColumn will be
#'   used.
#' @param sampleSpecific If TRUE, a BigWig will export for each sample-cell type
#'   combination.
#' @param saveFile Boolean. If TRUE, it will save to a BigWig. If FALSE, it will
#'   return the GRangesList.
#' @param numCores integer. Number of cores to parallelize peak-calling across
#'   multiple cell populations
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#'
#' @return countSE a SummarizedExperiment containing coverage for the given
#'   input cell populations.
#'
#' @examples
#' \dontrun{
#' MOCHA::exportCoverage(
#'   SampleTileObj = SampleTileMatrices,
#'   cellPopulations = "ALL",
#'   numCores = 30,
#'   sampleSpecific = FALSE
#' )
#' }
#'
#' @export
#' 

exportCoverage <- function(SampleTileObj,
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
  for (x in cellPopulations) {
    # MOCHA::getCoverage outputs a single list with two named items: "Accessibility" 
    # and "Insertions". 
    # This is saved to *_CoverageFiles.RDS in MOCHA::callOpenTiles
    if (type) { # Accessibility
      originalCovGRanges <- readRDS(paste(outDir, "/", x, "_CoverageFiles.RDS", sep = ""))$Accessibility
    } else { # Insertions
      originalCovGRanges <- readRDS(paste(outDir, "/", x, "_CoverageFiles.RDS", sep = ""))$Insertions
    }

    if (verbose) {
      message(stringr::str_interp("Extracting coverage from {x}."))
    }

    iterList <- lapply(subSamples, function(y) {
      originalCovGRanges[y]
    })
    # If not sample specific, take the average coverage across samples.
    # if it is sample specific, just subset down the coverage to the region of interest.
    if (!sampleSpecific) {
      cellPopSubsampleCov <- pbapply::pblapply(cl = cl, X = iterList, averageCoverage)
      names(cellPopSubsampleCov) <- subGroups
    } else {
      names(iterList) <- subGroups
      cellPopSubsampleCov <- unlist(iterList, recursive = FALSE)
      names(cellPopSubsampleCov) <- gsub("\\.", "__", names(cellPopSubsampleCov))
    }

    for (i in 1:length(cellPopSubsampleCov)) {
      fileName <- gsub(" ", "__", paste(x, groupColumn, names(cellPopSubsampleCov)[i], sep = "__"))
      if (!type) {
        fileName <- paste(fileName, "__Insertions", sep = "")
      }
      fileToWrite <- plyranges::filter(cellPopSubsampleCov[[i]], score != 0)
      plyranges::write_bigwig(fileToWrite, paste(dir, "/", fileName, ".bw", sep = ""))
    }

    if (saveFile) {
      names(cellPopSubsampleCov) <- x
      GRangesList1 <- append(GRangesList1, cellPopSubsampleCov)
    }
      
  }

  return(GRangesList1)
}


## helper function
## Efficiently generates average single basepair coverage for a given region
averageCoverage <- function(coverageList) {
  sampleCount <- length(coverageList)
  if (sampleCount > 1) {
    mergedCounts <- IRanges::stack(methods::as(coverageList, "GRangesList"))
    mergedCounts <- plyranges::compute_coverage(mergedCounts, weight = mergedCounts$score / sampleCount)
  } else {
    mergedCounts <- coverageList
  }


  return(mergedCounts)
}


#' @title \code{exportDifferentials}
#'
#' @description \code{exportDifferentials} exports the differential peaks
#'  output GRangesList output from \code{getDifferentialAccessibleTiles} to
#'  bigBed format for visualization in genome browsers.
#'
#' @param SampleTileObject The SummarizedExperiment object output from
#'   \code{getSampleTileMatrix}
#' @param DifferentialsGRList The GRangesList output from
#'   \code{getDifferentialAccessibleTiles}
#' @param outDir Desired output directory where bigBed files will be saved
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#'
#' @examples
#' \dontrun{
#' MOCHA::exportDifferentials(
#'   SampleTileObject = SampleTileMatrices,
#'   DifferentialsGRList,
#'   outDir = tempdir(),
#'   verbose = TRUE
#' )
#' }
#'
#' @export
#'
exportDifferentials <- function(SampleTileObject,
                                DifferentialsGRList,
                                outDir,
                                verbose = FALSE) {
  genome <- BSgenome::getBSgenome(S4Vectors::metadata(SampleTileObject)$Genome)
  
  for (i in seq_along(DifferentialsGRList)) {
    comparison_name <- names(DifferentialsGRList)[[i]]
    DiffPeaksGR <- DifferentialsGRList[[i]]

    # Set score and seqinfo for bigBed
    DiffPeaksGR$score <- 1
    GenomeInfoDb::seqinfo(DiffPeaksGR) <- GenomeInfoDb::seqinfo(genome)[GenomicRanges::seqnames(GenomeInfoDb::seqinfo(DiffPeaksGR))]

    outFile <- file.path(outDir, paste(comparison_name, sep = "__"))
    if (verbose) {
      message("Exporting: ", outFile)
    }

    # Output to bigbed
    rtracklayer::export.bb(DiffPeaksGR, paste0(outFile, ".bigBed"))
  }
}

#' @title \code{exportOpenTiles}
#'
#' @description \code{exportOpenTiles} exports the open tiles of a given cell
#'  population to bigBed file for visualization in genome browsers.
#'
#' @param SampleTileObject The SummarizedExperiment object output from
#'   \code{getSampleTileMatrix}
#' @param cellPopulation The name of the cell population to export
#' @param outDir Desired output directory where bigBed files will be saved
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#'
#' @examples
#' \dontrun{
#' MOCHA::exportOpenTiles(
#'   SampleTileObj = SampleTileObject,
#'   cellPopulation,
#'   outDir = tempdir(),
#'   verbose = TRUE
#' )
#' }
#'
#' @export
#'
exportOpenTiles <- function(SampleTileObject,
                            cellPopulation,
                            outDir,
                            verbose = FALSE) {
    
  genome <- BSgenome::getBSgenome(S4Vectors::metadata(SampleTileObject)$Genome)

  for (cellPopulation in names(SummarizedExperiment::assays(SampleTileObject))) {
    cellPopMatrix <- MOCHA::getCellPopMatrix(
      SampleTileObject,
      cellPopulation,
      NAtoZero = FALSE
    )
    for (sample in colnames(SampleTileObject)) {
      samplePeaks <- cellPopMatrix[,
        colnames(cellPopMatrix) %in% c(sample),
        drop = FALSE
      ]
      samplePeaksGR <- StringsToGRanges(rownames(samplePeaks))

      # Set score and seqinfo for bigBed
      samplePeaksGR$score <- 1
      GenomeInfoDb::seqinfo(samplePeaksGR) <- GenomeInfoDb::seqinfo(genome)[GenomicRanges::seqnames(GenomeInfoDb::seqinfo(samplePeaksGR))]

      sampleRow <- SummarizedExperiment::colData(SampleTileObject)[sample, ]
      pbmc_sample_id <- sampleRow[["Sample"]] # Enforced colname in callOpenTiles
      outFile <- file.path(outDir, paste(cellPopulation, pbmc_sample_id, sep = "__"))
      if (verbose) {
        message("Exporting: ", outFile)
      }

      # Output to bigbed
      rtracklayer::export.bb(samplePeaksGR, paste0(outFile, ".bigBed"))
    }
  }
}

#' @title \code{exportMotifs}
#'
#' @description \code{exportMotifs} exports a motif set GRanges from running
#'    \code{addMotifSet(returnSTM=FALSE)} to bigBed file files for visualization
#'    in genome browsers.
#'
#' @param SampleTileObject The SummarizedExperiment object output from
#'   \code{getSampleTileMatrix}
#' @param motifsGRanges A GRanges containing motif annotations, typically from
#'   \code{addMotifSet(returnSTM=FALSE)}
#' @param filterByOpenTiles Boolean. If TRUE, a bigBed file will be exported
#'   for each cell population with motifs filtered to those occurring
#'   only in open tiles.
#' @param outDir Desired output directory where bigBed files will be saved
#' @param motifSetName Optional, a name indicating the motif set. Used to name
#'   files in the specified \code{outdir}. Default is "motifs".
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#'
#' @examples
#' \dontrun{
#' MOCHA::exportMotifs(
#'   SampleTileObj = SampleTileMatrices,
#'   motifsGRanges,
#'   motifSetName = "CISBP",
#'   filterByOpenTiles = FALSE,
#'   outDir = tempdir(),
#'   verbose = TRUE
#' )
#' }
#'
#' @export
#'
exportMotifs <- function(SampleTileObject,
                         motifsGRanges,
                         motifSetName = "motifs",
                         filterByOpenTiles = FALSE,
                         outDir,
                         verbose = FALSE) {

  # Map over the seqinfo from our genome to the motifsGRanges
  # Required for bigBed export
  genome <- BSgenome::getBSgenome(S4Vectors::metadata(SampleTileObject)$Genome)
  GenomeInfoDb::seqinfo(motifsGRanges) <- GenomeInfoDb::seqinfo(genome)[GenomicRanges::seqnames(GenomeInfoDb::seqinfo(motifsGRanges))]

  # # Truncate trailing zeroes from score https://www.biostars.org/p/235193/
  # # (Error : Trailing characters parsing integer in field 4 line 1 of text, got 10.3405262378482)
  motifsGRanges$score <- floor(motifsGRanges$score)

  # Filter out negative scores https://github.com/jhkorhonen/MOODS/issues/12
  motifsGRanges <- motifsGRanges[motifsGRanges$score > 0]


  if (filterByOpenTiles) {
    # Split motifsGRangesFiltered by cell population, filtered to open tiles in cell population
    allPeaks <- SummarizedExperiment::rowRanges(SampleTileObject)
    samplePeakTable <- GenomicRanges::mcols(allPeaks)
    for (celltype in names(samplePeakTable)) {
      # Filter rows with boolean index
      samplePeaksGR <- allPeaks[samplePeakTable[[celltype]], ]
      cellTypePeakMotifs <- plyranges::filter_by_overlaps(motifsGRanges, samplePeaksGR)

      outFile <- file.path(outDir, paste(celltype, motifSetName, "motifset.bigBed", sep = "__"))
      if (verbose) {
        message("Exporting: ", outFile)
      }
      # Output to bigbed
      rtracklayer::export.bb(cellTypePeakMotifs, outFile)
    }
  } else {
    outFile <- file.path(outDir, paste(motifSetName, "motifset.bigBed", sep = "__"))
    if (verbose) {
      message("Exporting: ", outFile)
    }
    # Output to bigbed
    rtracklayer::export.bb(motifsGRanges, outFile)
  }
}
