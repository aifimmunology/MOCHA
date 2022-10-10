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
#' @param TxDb is an AnnotationDbi object with transcript info for the organism.
#' @param Org is the genome-wide annotation package for your organism.
#' @param outDir is a string describing the output directory for coverage files per sample/celltype.
#'   Must be a complete directory string. Default is NULL, in which case it'll pull out a directory from ArchR
#'   and make a new fold named MOCHA for saving files.
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
                          TxDb,
                          Org,
                          outDir = NULL,
                          fast = FALSE,
                          numCores = 30,
                          verbose = TRUE,
                          force = FALSE,
                          studySignal=NULL,
                          sampleSpecific=TRUE
                         ) {
  
  if (is.null(outDir)) {
    ## Generate folder within ArchR for outputting results
    outDir <- paste(ArchR::getOutputDirectory(ArchRProj), "/MOCHA", sep = "")
  }

  if (!file.exists(outDir)) {
    message(stringr::str_interp("Creating directory for MOCHA at ${outDir}"))
    dir.create(outDir)
  }
  
  if (fast) {
    tileResults <- callOpenTilesFast(
      ArchRProj,
      cellPopLabel,
      cellPopulations,
      TxDb,
      Org,
      outDir,
      numCores,
      force
    )
    return(tileResults)
  }

  # Get cell metadata and blacklisted regions from ArchR Project
  cellColData <- ArchR::getCellColData(ArchRProj)
  blackList <- ArchR::getBlacklist(ArchRProj)

  # Get cell populations
  cellTypeLabelList <- cellColData[, cellPopLabel]
  cellCounts <- table(cellTypeLabelList[!is.na(cellTypeLabelList)])


  if (all(cellPopulations == "ALL")) {
    cellPopulations <- names(cellCounts)
  } else if (!all(cellPopulations %in% names(cellCounts))) {
    stop("Error: cellPopulations not all found in ArchR project.")
  }

  #Save the cell number per population-sample in the metadata
  allCellCounts <- table(cellColData[, "Sample"], cellTypeLabelList)
  allCellCounts <- allCellCounts[, cellPopulations]

  # Genome and TxDb annotation info is added to the metadata of
  # the final MultiAssayExperiment for downstream analysis
  genome <- ArchR::validBSgenome(ArchR::getGenome(ArchRProj))
  AnnotationDbi::saveDb(TxDb, paste(outDir, "/TxDb.sqlite", sep = ""))
  AnnotationDbi::saveDb(Org, paste(outDir, "/Org.sqlite", sep = ""))

  # Main loop over all cell populations
  experimentList <- list()
  for (cellPop in cellPopulations) {
    message(stringr::str_interp("Calling open tiles for cell population ${cellPop}"))

    # Get our fragments for this cellPop
    frags <- scMACS::getPopFrags(
      ArchRProj = ArchRProj,
      metaColumn = cellPopLabel,
      cellSubsets = cellPop,
      region = NULL,
      numCores = numCores,
      sampleSpecific = sampleSpecific,
      NormMethod = "nfrags",
      blackList = NULL,
	  verbose = verbose,
      overlapList = 50
    )

    # Check for and remove celltype-sample groups for which there are no fragments.
    fragsNoNull <- frags[lengths(frags) != 0]
    emptyFragsBool <- !(names(frags) %in% names(fragsNoNull))
    emptyGroups <- names(frags)[emptyFragsBool]
    emptyGroups <- gsub("__.*", "", emptyGroups)
	
	rm(fragsNoNull)
	
    if (length(emptyGroups) == 0) {
      warning(
        "The following celltype#sample groupings have no fragments",
        "and will be ignored: ", emptyGroups
      )
    }

    # Simplify sample names to remove everything before the first "#"
    sampleNames <- gsub("__.*", "", gsub(".*#", "", names(frags)))
    names(frags) <- sampleNames

    # Calculate normalization factors as the number of fragments for each celltype_samples
    normalization_factors <- as.integer(sapply(frags, length))

    # save coverage files to folder.
    if (!file.exists(paste(outDir, "/", cellPop, "_CoverageFiles.RDS", sep = "")) | force) {
      covFiles <- scMACS:::getCoverage(
        popFrags = frags,
        normFactor = normalization_factors / 10^6,
        filterEmpty = FALSE,
        numCores = numCores, TxDb = TxDb
      )
      saveRDS(covFiles, paste(outDir, "/", cellPop, "_CoverageFiles.RDS", sep = ""))
      rm(covFiles)
    }

	#rm(frags)
	
    # Add prefactor multiplier across datasets
    if(is.null(studySignal)){
        message('calculating study signal on ArchR project. Make sure data contains all cell populations')
        curr_frags_median <- median(cellColData$nFrags)
        
    }
    study_prefactor <- 3668 / studySignal # Training median

    # This mclapply will parallelize over each sample within a celltype.
    # Each arrow is a sample so this is allowed
    # (Arrow files are locked - one access at a time) 
    tilesGRangesList <- parallel::mclapply(
      1:length(frags),
      function(x) {
        scMACS:::callTilesBySample(
          blackList = blackList,
          returnAllTiles = TRUE,
          numCores = numCores,
          totalFrags = normalization_factors[x],
          fragsList = frags[[x]],
		  verbose = verbose,
          StudypreFactor = study_prefactor
        )
      },
      mc.cores = numCores
    )
      
    if(sampleSpecific){
        
        names(tilesGRangesList) <- sampleNames
        
    } else{
        names(tilesGRangesList) <- cellPop
        
    }

    # Cannot make peak calls with < 5 cells (see make_prediction.R)
    # so NULL will occur for those samples. We need to fill in dummy data so that we
    # preserve the existence of the sample, while also not including any information from it. 
    
    emptyGroups <- which(unlist(lapply(tilesGRangesList, is.null)))
    
    if(length(emptyGroups) > 0){
    warning(
      "The following celltype#sample groupings have too few celltypes (<5)",
      "and will be ignored: ", names(tilesGRangesList)[emptyGroups]
      ) 
    }
    for(i in emptyGroups){
      # This is an empty region placeholder that represents an empty sample
      # And is functionally ignored downstream but required for
      # the RaggedExperiment structure
      tilesGRangesList[[i]] <- data.table(
        tileID = 'chr1:1-499', seqnames = 'chr1', start = 1,
        end = 499, strand ='*', TotalIntensity = 0, maxIntensity = 0,
        numCells = 0, Prediction = 0, PredictionStrength = 0, peak = FALSE)
    }

    # Package rangeList into a RaggedExperiment
    ragExp <- RaggedExperiment::RaggedExperiment(
      tilesGRangesList
    )
    # And add it to the experimentList for this cell population
    experimentList <- append(experimentList, ragExp)
    # experimentList <- append(experimentList, tilesGRangesListNoNull)
  }

  # Create sample metadata from cellColData using util function
  # "Sample" is the enforced col in ArchR containing the
  # sample IDs which correspond to the arrow files.
  # We are assuming samples are synonymous to wells
  # (If hashed, samples were un-hashed during ArchR project generation)
  
  if(sampleSpecific){
      sampleData <- suppressWarnings(
        scMACS:::sampleDataFromCellColData(cellColData, 
                                           sampleLabel = "Sample")
      )
    } else{
      sampleData <- suppressWarnings(
        scMACS:::sampleDataFromCellColData(cellColData, 
                                           sampleLabel = cellPopLabel))
      
  }

  # Add experimentList to MultiAssayExperiment
  names(experimentList) <- cellPopulations
   
  tileResults <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = experimentList,
    colData = sampleData,
    metadata = list('CellCounts' = allCellCounts,
      "Genome" = genome, "TxDb" = paste(outDir, "/TxDb.sqlite", sep = ""),
      "Org" = paste(outDir, "/Org.sqlite", sep = ""), "Directory" = outDir
    )
  )

  return(tileResults)
}

##### Same function, but runs faster with much, much more RAM usage.
callOpenTilesFast <- function(ArchRProj,
                              cellPopLabel,
                              cellPopulations = "ALL",
                              TxDb = NULL,
                              Org = NULL,
                              outDir = NULL,
                              numCores = 30,
                              force = FALSE) {
  
  # Genome and TxDb annotation info is added to the metadata of 
  # the final MultiAssayExperiment for downstream analysis
  genome <- ArchR::validBSgenome(ArchR::getGenome(ArchRProj))
  AnnotationDbi::saveDb(TxDb, paste(outDir, "/TxDb.sqlite", sep = ""))
  AnnotationDbi::saveDb(Org, paste(outDir, "/Org.sqlite", sep = ""))
  
  # Get cell metadata and blacklisted regions from ArchR Project
  cellColData <- ArchR::getCellColData(ArchRProj)
  blackList <- ArchR::getBlacklist(ArchRProj)

  # Get frags grouped by cell population and sample
  # This will also validate the input cellPopulations
  frags <- scMACS::getPopFrags(
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

  if (length(emptyGroups) == 0) {
    warning(
      "The following celltype#sample groupings have no fragments",
      "and will be ignored: ", emptyGroups
    )
  }

  prefilterCellPops <- unique(sapply(strsplit(names(frags), "#"), `[`, 1))
  if (any(prefilterCellPops != cellPopulations)) {
    # Removed cell populations must be filtered from cellPopulations
    postfilterCellPops <- unique(sapply(strsplit(names(frags), "#"), `[`, 1))
    # Match these labels by index to the original cellPopulatoins
    retainedPopsBool <- (prefilterCellPops %in% postfilterCellPops)
    finalCellPopulations <- cellPopulations[retainedPopsBool]
  } else {
    finalCellPopulations <- cellPopulations
  }

  # Split the fragments list into a list of lists per cell population
  splitFrags <- scMACS:::splitFragsByCellPop(frags)

  # getPopFrags needs to replace spaces for underscores, so
  # here we rename the fragments with the original cell populations labels.
  # Fragments maintain their order by cell population.
  names(splitFrags) <- finalCellPopulations
  message("Cell populations for peak calling: ", paste(names(splitFrags), collapse = ", "))


  if (is.null(outDir)) {
    # Generate folder within ArchR project for outputting results
    outDir <- paste(ArchR::getOutputDirectory(ArchRProj), "/MOCHA", sep = "")
  }
  if (!dir.exists(outDir)) {
    dir.create(outDir)
  }

  if (!file.exists(outDir)) {
    message(stringr::str_interp("Creating directory for MOCHA at ${outDir}"))
    dir.create(outDir)
  }

  # Main loop over all cell populations
  experimentList <- list()
  for (cellPop in names(splitFrags)) {
    message(stringr::str_interp("Calling open tiles for cell population ${cellPop}"))

    # Get our fragments for this cellPop
    popFrags <- unlist(splitFrags[[cellPop]])

    # Simplify sample names to remove everything before the first "."
    sampleNames <- gsub("^([^.]+).", "", names(popFrags))
    names(popFrags) <- sampleNames

    # Calculate normalization factors as the number of fragments for each celltype_samples
    normalization_factors <- as.integer(sapply(popFrags, length))

    # save coverage files to folder.
    if (!file.exists(paste(outDir, "/", cellPop, "_CoverageFiles.RDS", sep = "")) | force) {
      covFiles <- scMACS:::getCoverage(
        popFrags = popFrags,
        normFactor = normalization_factors / 10^6,
        filterEmpty = FALSE,
        numCores = numCores, TxDb = TxDb
      )
      saveRDS(covFiles, paste(outDir, "/", cellPop, "_CoverageFiles.RDS", sep = ""))
      rm(covFiles)
    }


    # Add prefactor multiplier across datasets
    curr_frags_median <- median(cellColData$nFrags)
    study_prefactor <- 3668 / curr_frags_median # Training median

    # This mclapply will parallelize over each sample within a celltype.
    # Each arrow is a sample so this is allowed
    # (Arrow files are locked - one access at a time)
    tilesGRangesList <- parallel::mclapply(
      1:length(popFrags),
      function(x) {
        scMACS:::callTilesBySample(
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
  browser()
    
  # Add experimentList to MultiAssayExperiment
  names(experimentList) <- names(splitFrags)
  
  tileResults <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = experimentList,
    colData = sampleData,
    metadata = list(
      "Genome" = genome, "TxDb" = paste(outDir, "/TxDb.sqlite", sep = ""),
      "Org" = paste(outDir, "/Org.sqlite", sep = ""), "Directory" = outDir
    )
  )

  return(tileResults)
}
