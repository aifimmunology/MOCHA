#' @title \code{callOpenTiles}
#'
#' @description \code{callOpenTiles} is the main peak-calling function in scMACS
#'   that serves as a wrapper function to call peaks provided a set of fragment
#'   files and an ArchR Project for meta-data purposes
#'
#'
#' @param ArchRProj an ArchR Project
#' @param cellPopLabel string indicating which column in the ArchRProject
#'   metadata contains the cell population label.
#' @param cellPopulations vector of strings. Cell subsets for which to call
#'   peaks. This list of group names must be identical to names that appear in
#'   the ArchRProject metadata.  Optional, if cellPopulations='ALL', then peak
#'   calling is done on all cell populations in the ArchR project metadata.
#'   Default is 'ALL'.
#' @param numCores integer. Number of cores to parallelize peak-calling across
#'   multiple cell populations.
#'
#' @return tileResults A MultiAssayExperiment object containing ranged data
#'   for each tile 
#'
#' @details The technical details of the algorithm are found in XX.
#'
#' @references XX
#'
#' @export

callOpenTiles <- function(ArchRProj,
                          cellPopLabel,
                          cellPopulations = "ALL",
                          numCores = 30) {

  # Get cell metadata and blacklisted regions from ArchR Project
  cellColData <- ArchR::getCellColData(ArchRProj)
  blackList <- ArchR::getBlacklist(ArchRProj)

  # Get frags grouped by cell population and sample
  # This will also validate the input cellPopulations
  frags <- getPopFrags(
    ArchRProj = ArchRProj,
    metaColumn = cellPopLabel,
    cellSubsets = cellPopulations,
    region = NULL,
    numCores = numCores,
    sampleSpecific = TRUE,
    NormMethod = "nfrags",
    blackList = NULL,
    overlapList = 50
  )

  # Check for and remove celltype-sample groups for which there are no fragments.
  fragsNoNull <- frags[lengths(frags) != 0]
  emptyFragsBool <- !(names(frags) %in% names(fragsNoNull[7:10]))
  emptyGroups <- names(frags)[emptyFragsBool]
  emptyGroups <- gsub("__.*", "", emptyGroups)
  
  if (length(emptyGroups)==0){
    warning("The following celltype#sample groupings have no fragments",
            "and will be ignored: ", emptyGroups)
  }
  
  prefilterCellPops <- unique(sapply(strsplit(names(frags),"#"), `[`, 1))
  if (any(prefilterCellPops != cellPopulations)){
    # Removed cell populations must be filtered from cellPopulations
    postfilterCellPops <- unique(sapply(strsplit(names(fragsNoNull),"#"), `[`, 1))
    # Match these labels by index to the original cellPopulatoins
    retainedPopsBool <- (prefilterCellPops %in% postfilterCellPops)
    finalCellPopulations <- cellPopulations[retainedPopsBool]
  } else {
    finalCellPopulations <- cellPopulations
  }
  
  # Split the fragments list into a list of lists per cell population
  splitFrags <- scMACS:::splitFragsByCellPop(fragsNoNull)
  
  # getPopFrags needs to replace spaces for underscores, so
  # here we rename the fragments with the original cell populations labels.
  # Fragments maintain their order by cell population.
  names(splitFrags) <- finalCellPopulations
  message("Cell populations for peak calling: ", paste(names(splitFrags), ""))

  # Main loop over all cell populations
  experimentList <- list()
  for (cellPop in names(splitFrags)) {
    message(str_interp("Calling open tiles for cell population ${cellPop}"))

    # Get our fragments for this cellPop
    popFrags <- unlist(splitFrags[[cellPop]])

    # Simplify sample names to remove everything before the first "."
    sampleNames <- gsub("^([^.]+).", "", names(popFrags))
    names(popFrags) <- sampleNames

    # Calculate normalization factors as the number of fragments for each celltype_samples
    normalization_factors <- as.integer(sapply(popFrags, length))

    # Add prefactor multiplier across datasets
    curr_frags_median <- median(cellColData$nFrags)
    study_prefactor <- 3668 / curr_frags_median # Training median

    # This mclapply will parallelize over each sample within a celltype.
    # Each arrow is a sample so this is allowed
    # (Arrow files are locked - one access at a time)
    tilesGRangesList <- parallel::mclapply(
      1:length(popFrags),
      function(x) {
        callTilesBySample(
          blackList = blackList,
          returnAllTiles = TRUE,
          numCores = numCores,
          totalFrags = normalization_factors[x],
          fragsList = popFrags[[x]],
          StudypreFactor = study_prefactor
        )
      },
      mc.cores = numCores
    )

    names(tilesGRangesList) <- sampleNames
    
    # Cannot make peak calls with < 5 cells (see make_prediction.R) 
    # so NULL will occur for those samples
    tilesGRangesListNoNull <- BiocGenerics::Filter(Negate(is.null), tilesGRangesList)
    # TODO: Add warning message about removed samples for this celltype.
    
    # Package rangeList into a RaggedExperiment
    ragExp <- RaggedExperiment::RaggedExperiment(
      tilesGRangesListNoNull
    )

    # And add it to the experimentList for this cell population
    experimentList <- append(experimentList, ragExp)
  }

  # Create sample metadata from cellColData using util function
  # "Sample" is the enforced col in ArchR containing the 
  # sample IDs which correspond to the arrow files.
  # We are assuming samples are synonymous to wells
  # (If hashed, samples were un-hashed during ArchR project generation)
  sampleData <- suppressWarnings(
    scMACS:::sampleDataFromCellColData(cellColData, sampleLabel = "Sample")
  )

  genome <- ArchR::validBSgenome(ArchR::getGenome(ArchRProj))
 
  # Add experimentList to MultiAssayExperiment
  names(experimentList) <- names(splitFrags)
  tileResults <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = experimentList,
    colData = sampleData,
    metadata = list('Genome' = genome)
  )

  return(tileResults)
}
