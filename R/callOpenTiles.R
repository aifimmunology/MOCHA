#' @title \code{callOpenTiles} Perform peak-calling on a set of fragments or an 
#'   ArchR Project.
#'
#' @description \code{callOpenTiles} is the main peak-calling function in MOCHA
#'   that serves as a wrapper function to call peaks provided a set of fragment
#'   files and an ArchR Project for meta-data purposes
#'
#'
#' @param ATACFragments an ArchR Project, or a GRangesList of fragments
#' @param cellPopLabel string indicating which column in the ArchRProject
#'   metadata contains the cell population label.
#' @param cellPopulations vector of strings. Cell subsets for which to call
#'   peaks. This list of group names must be identical to names that appear in
#'   the ArchRProject metadata.  Optional, if cellPopulations='ALL', then peak
#'   calling is done on all cell populations in the ArchR project metadata.
#'   Default is 'ALL'.
#' @param cellColData A DataFrame containing cell-level metadata and a 'Sample' 
#'   column
#' @param additionalCellData Name of columns from the cellColData that the
#'   user wants summarized by pseodubulk sample and saved in Sample-Cell type metadata (in addition to Cell Counts). Data must numeric. 
#'   Default is 'nFrags', based on ArchR's metadata interpretation.
#' @param blackList A GRanges of blacklisted regions
#' @param genome A BSgenome object, or the full name of an installed 
#'   BSgenome data package, or a short string specifying the name of an NCBI 
#'   assembly (e.g. "GRCh38", "TAIR10.1", etc...) or UCSC genome (e.g. "hg38", 
#'   "bosTau9", "galGal6", "ce11", etc...). The supplied short string must refer
#'   unambiguously to an installed BSgenome data package. See 
#'   \link[BSgenome]{getBSgenome}.
#' @param studySignal The median signal (number of fragments) in your study. If 
#'   not set, this will be calculated using the input ArchR project but relies 
#'   on the assumption that the ArchR project encompasses your whole study (i.e. 
#'   is not a subset).
#' @param TxDb The exact package name of a TxDb-class transcript annotation 
#'   package for your organism (e.g. "TxDb.Hsapiens.UCSC.hg38.refGene"). This 
#'   must be installed. See 
#'   \href{https://bioconductor.org/packages/release/data/annotation/}{
#'   Bioconductor AnnotationData Packages}.
#' @param OrgDb The exact package name of a OrgDb-class genome wide annotation 
#'   package for your organism (e.g. "org.Hs.eg.db"). This must be installed. 
#'   See \href{https://bioconductor.org/packages/release/data/annotation/}{
#'   Bioconductor AnnotationData Packages}
#' @param outDir is a string describing the output directory for coverage files.
#'   Must be a complete directory string. With ArchR input, set outDir to NULL 
#'   to create a directory within the input ArchR project directory named MOCHA
#'   for saving files.
#' @param numCores integer. Number of cores to parallelize peak-calling across
#'   multiple cell populations.
#' @param force Optional, whether to force creation of coverage files if they
#'   already exist. Default is FALSE.
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#'
#' @return tileResults A MultiAssayExperiment object containing ranged data for
#'   each tile
#'
#' @examples
#' \dontrun{
#' # Starting from an ArchR Project:
#' tileResults <- MOCHA::callOpenTiles(
#'   ArchRProj = myArchRProj,
#'   cellPopLabel = "celltype_labeling",
#'   cellPopulations = "CD4",
#'   TxDb = "TxDb.Hsapiens.UCSC.hg38.refGene",
#'   OrgDb = "org.Hs.eg.db",
#'   numCores = 1
#' )
#' }
#' \donttest{
#' # Starting from GRangesList
#' if (
#'   requireNamespace(BSgenome.Hsapiens.UCSC.hg19) &&
#'     requireNamespace(TxDb.Hsapiens.UCSC.hg38.refGene) &&
#'     requireNamespace(org.Hs.eg.db)
#' ) {
#'   tiles <- MOCHA::callOpenTiles(
#'     ATACFragments = MOCHA::exampleFragments,
#'     cellColData = MOCHA::exampleCellColData,
#'     blackList = MOCHA::exampleBlackList,
#'     genome = "BSgenome.Hsapiens.UCSC.hg19",
#'     TxDb = "TxDb.Hsapiens.UCSC.hg38.refGene",
#'     OrgDb = "org.Hs.eg.db",
#'     outDir = tempdir(),
#'     cellPopLabel = "Clusters",
#'     cellPopulations = c("C2", "C5"),
#'     numCores = 1
#'   )
#' }
#' }
#'
#' @export
#' @docType methods
#' @rdname callOpenTiles-methods
setGeneric(
  "callOpenTiles",
  function(ATACFragments,
           cellColData,
           blackList,
           genome,
           cellPopLabel,
           cellPopulations = "ALL",
           additionalCellData = "nFrags",
           studySignal = NULL,
           TxDb,
           OrgDb,
           outDir,
           numCores = 30,
           verbose = FALSE,
           force = FALSE) {
    standardGeneric("callOpenTiles")
  },
  signature = "ATACFragments"
)


.callOpenTiles_default <- function(ATACFragments,
                                   cellColData,
                                   blackList,
                                   genome,
                                   cellPopLabel,
                                   cellPopulations = "ALL",
                                   additionalCellData = NULL,
                                   studySignal = NULL,
                                   TxDb,
                                   OrgDb,
                                   outDir,
                                   numCores = 30,
                                   verbose = FALSE,
                                   force = FALSE) {
  
  # Load Genome
  genome <- BSgenome::getBSgenome(genome)
  
  validGRanges <- sapply(ATACFragments, function(obj) {
    methods::is(obj, "GRanges")
  })
  if (!all(validGRanges)) {
    stop(
      "Invalid ATACFragments. ATACFragments must be a list of GRanges or ",
      "GRangesList."
    )
  }
  if (!("Sample" %in% colnames(cellColData))) {
    stop("Invalid cellColData. cellColData must contain column 'Sample'.")
  }
  if (!(cellPopLabel %in% colnames(cellColData))) {
    stop(
      stringr::str_interp("cellPopLabel {cellPopLabel} not found in "),
      "cellColData. cellColData must contain column cellPopLabel."
    )
  }

  # Save the fragment number per population-sample
  cellPopList <- sapply(names(ATACFragments), function(x) {
    unlist(stringr::str_split(x, "#"))[1]
  })
  sampleList <- sapply(names(ATACFragments), function(x) {
    sample <- unlist(stringr::str_split(x, "#"))[2]
    unlist(stringr::str_split(sample, "__"))[1]
  })

  

  if (!is.null(additionalCellData) & !all(additionalCellData %in% colnames(cellColData))) {
    stop(
      stringr::str_interp("additionalCellData {additionalCellData} not found in "),
      "cellColData. cellColData must contain columns for additionalCellData."
    )
  }

  # Save additional cell metrics. But first we need to test that these additional columns are numeric data. 
  test <- all(unlist(lapply(additionalCellData, function(x) is.numeric(cellColData[,x]))))
  if(!is.null(additionalCellData) & test){
      additionalMetaData <- lapply(additionalCellData, function(x){

        cellData_Summary <- as.data.frame(cellColData) %>%
          dplyr::group_by(!!as.name(cellPopLabel), Sample) %>%
          dplyr::summarize(tmp_value = mean({{x}}), .groups = "drop")
      
        tidyr::pivot_wider(
          cellData_Summary,
          names_from = tidyselect::all_of({{cellPopLabel}}),
          values_from = tmp_value
        )
      })
      names(additionalMetaData) = additionalCellData
  }else if(is.null(additionalCellData)){
    additionalMetaData = NULL
  }else if(!test){
     stop(
      stringr::str_interp("additionalCellData {additionalCellData} from "),
      "cellColData is not all numeric and so cannot be summarized. All columns for additionalCellData must be numeric."
    )
  }

  .callOpenTiles(
    ATACFragments,
    cellColData,
    blackList,
    genome,
    cellPopLabel,
    cellPopulations,
    additionalCellData = additionalMetaData,
    studySignal,
    TxDb,
    OrgDb,
    outDir,
    numCores,
    verbose,
    force,
    useArchR = FALSE
  )
}
#' @rdname callOpenTiles-methods
#' @aliases callOpenTiles,GRangesList-method
setMethod(
  "callOpenTiles",
  signature(ATACFragments = "GRangesList"),
  .callOpenTiles_default
)
#' @rdname callOpenTiles-methods
#' @aliases callOpenTiles,list-method
setMethod(
  "callOpenTiles",
  signature(ATACFragments = "list"),
  .callOpenTiles_default
)


#' @rdname callOpenTiles-methods
#' @aliases callOpenTiles,ArchRProject-method
.callOpenTiles_ArchR <- function(ATACFragments,
                                 cellPopLabel,
                                 cellPopulations = "ALL",
                                 additionalCellData = NULL,
                                 studySignal = NULL,
                                 TxDb,
                                 OrgDb,
                                 outDir = NULL,
                                 numCores = 30,
                                 verbose = FALSE,
                                 force = FALSE) {
  
  # Load Genome
  genome <- ArchR::validBSgenome(ArchR::getGenome(ATACFragments))
  
  if (is.null(outDir)) {
    outDir <- paste(
      ArchR::getOutputDirectory(ATACFragments),
      "/MOCHA",
      sep = ""
    )
  }

  if (!file.exists(outDir)) {
    if (verbose) {
      message(stringr::str_interp("Creating directory for MOCHA at ${outDir}"))
    }
    dir.create(outDir)
  }

  # Get cell metadata and blacklisted regions from ArchR Project
  cellColData <- ArchR::getCellColData(ATACFragments)
  blackList <- ArchR::getBlacklist(ATACFragments)
  
  if (!(cellPopLabel %in% colnames(cellColData))) {
    stop(
      stringr::str_interp("cellPopLabel {cellPopLabel} not found in "),
      "cellColData. cellColData must contain column cellPopLabel."
    )
  }

  if (!is.null(additionalCellData) & !all(additionalCellData %in% colnames(cellColData))) {
    stop(
      stringr::str_interp("additionalCellData {additionalCellData} not found in "),
      "cellColData. cellColData must contain columns for additionalCellData."
    )
  }

  # Save additional cell metrics. 
  test <- all(unlist(lapply(additionalCellData, function(x) is.numeric(cellColData[,x]))))
  if(!is.null(additionalCellData) & test){
      additionalMetaData <- lapply(additionalCellData, function(x){

        cellData_Summary <- as.data.frame(cellColData) %>%
          dplyr::group_by(!!as.name(cellPopLabel), Sample) %>%
          dplyr::summarize(tmp_value = mean({{x}}), .groups = "drop")
      
        tidyr::pivot_wider(
          cellData_Summary,
          names_from = tidyselect::all_of({{cellPopLabel}}),
          values_from = tmp_value
        )
      })
      names(additionalMetaData) = additionalCellData
  }else if(is.null(additionalCellData)){
    additionalMetaData = NULL
  }else if(!test){
     stop(
      stringr::str_interp("additionalCellData {additionalCellData} from "),
      "cellColData is not all numeric and so cannot be summarized. All columns for additionalCellData must be numeric."
    )
  }

  .callOpenTiles(
    ATACFragments,
    cellColData,
    additionalCellData = additionalMetaData,
    blackList,
    genome,
    cellPopLabel,
    cellPopulations,
    studySignal,
    TxDb,
    OrgDb,
    outDir,
    numCores,
    verbose,
    force,
    useArchR = TRUE
  )
}
setMethod(
  "callOpenTiles",
  signature(ATACFragments = "ArchRProject"),
  .callOpenTiles_ArchR
)


.callOpenTiles <- function(
    ATACFragments,
    cellColData,
    blackList,
    genome,
    cellPopLabel,
    cellPopulations,
    additionalCellData,
    studySignal,
    TxDb,
    OrgDb,
    outDir,
    numCores,
    verbose,
    force,
    useArchR) {
  
  # Load databases and save names for use in metadata
  TxDbName <- TxDb
  OrgDbName <- OrgDb
  TxDb <- getAnnotationDbFromInstalledPkgname(dbName=TxDb, type="TxDb")
  OrgDb <- getAnnotationDbFromInstalledPkgname(dbName=OrgDb, type="OrgDb")
 
  # Get cell populations
  cellTypeLabelList <- cellColData[, cellPopLabel]

  # Save the cell number per population-sample in the metadata
  cellCounts <- as.data.frame(table(cellColData[, "Sample"], cellTypeLabelList))
  names(cellCounts) <- c('Sample', 'CellType', 'CellCount')
  cellCounts <-  tidyr::pivot_wider(cellCounts, id_cols = 'CellType', names_from= 'Sample', values_from = 'CellCount') 
  allCellCounts <- as.data.frame(cellCounts[,-1])
  rownames(allCellCounts) <- cellCounts$CellType

  ## Create a dummy data.frame for storing fragment information with each iteration. 
  allFragmentCounts <- allCellCounts
  allFragmentCounts[!is.na(allFragmentCounts)] <- NA

  if (all(cellPopulations == "ALL")) {
    cellPopulations <- colnames(allCellCounts)
  } else {
    allCellCounts <- allCellCounts[rownames(allCellCounts) %in% cellPopulations,]
    allFragmentCounts <- allFragmentCounts[rownames(allCellCounts) %in% cellPopulations,]
  }

  # Add prefactor multiplier across datasets
  if (is.null(studySignal)) {
    if (verbose) {
      message(
        "studySignal was not provided. ",
        "Calculating study signal on cellColData as the median ",
        "nFrags with the assumption that all cell populations are ",
        "present in cellColData."
      )
    }
    if (!("nFrags" %in% colnames(cellColData))) {
      stop(
        "cellColData is missing fragment count information. ",
        "To calculate study signal, cellColData must contain a column",
        " 'nFrags' representing the number of fragments per cell. ",
        "Alternatively, provide a value for the parameter 'studySignal'."
      )
    }
    studySignal <- stats::median(cellColData$nFrags)
  }
  study_prefactor <- 3668 / studySignal # Training median

  # Main loop over all cell populations
  experimentList <- list()
  for (cellPop in cellPopulations) {
    if (verbose) {
      message(stringr::str_interp(
        "Calling open tiles for cell population ${cellPop}"
      ))
    }

    cl <- parallel::makeCluster(numCores)

    if (useArchR) {
      parallel::clusterEvalQ(cl, {
        library(ArchR)
        library(rhdf5)
      })
      # Get our fragments for this cellPop
      frags <- MOCHA::getPopFrags(
        ArchRProj = ATACFragments,
        cellPopLabel = cellPopLabel,
        cellSubsets = cellPop,
        numCores = cl,
        verbose = verbose
      )
    } else {
      frags <- ATACFragments[grep(cellPop, names(ATACFragments))]
      if (length(frags) == 0) {
        stop(
          "Provided ATACFragments does not contain any sample fragments ",
          stringr::str_interp("belonging to cellPopulations: ${cellPop}. "),
          "ATACFragments must have at lease one sample GRanges for each ",
          "cell population."
        )
      }
    }

    # Simplify sample names to remove celltype and normalization factor
    # Removes everything before the first "#" and after the first "__"
    # This is a convention from MOCHA::getPopFrags which is ArchR-dependent
    # E.g. C1#PBMCSmall__0.527239 becomes PBMCSmall
    sampleNames <- gsub("__.*", "", gsub(".*#", "", names(frags)))
    names(frags) <- sampleNames

    # Calculate normalization factors as the number of fragments for
    # each celltype_samples. Then save it in metadata. 
    normalization_factors <- as.integer(lengths(frags))

    #Save the number of fragments into the fragment count. 
    allFragmentCounts[cellPop,sampleNames] = normalization_factors

    ## Create cores to parallelize over. 
    cl <- parallel::makeCluster(numCores)
    # save coverage files to folder.
    # This doesn't include empty samples and might break. We may need to
    # reconsider how getCoverage works and add empty samples before this step.
    if (!file.exists(
      paste(outDir, "/", cellPop, "_CoverageFiles.RDS", sep = "")
    ) || force) {
      if (verbose) {
        message(stringr::str_interp(
          "Saving coverage files for cell population ${cellPop}"
        ))
      }
      covFiles <- getCoverage(
        popFrags = frags,
        normFactor = normalization_factors / 10^6,
        filterEmpty = FALSE,
        cl = cl, TxDb = TxDb
      )
      saveRDS(
        covFiles,
        paste(outDir, "/", cellPop, "_CoverageFiles.RDS", sep = "")
      )
      rm(covFiles)
    }
    # This pbapply will parallelize over each sample within a celltype.
    # Each arrow is a sample so this is allowed
    # (Arrow files are locked - one access at a time)
    iterList <- lapply(seq_along(frags), function(x) {
      list(blackList, frags[[x]], verbose, study_prefactor)
    })
    tilesGRangesList <- pbapply::pblapply(
      cl = cl,
      X = iterList,
      FUN = simplifiedTilesBySample
    )
    parallel::stopCluster(cl)
    names(tilesGRangesList) <- names(frags)

    # Where samples have no cells, add an empty GRanges placeholder
    if (!all(rownames(allCellCounts) %in% names(tilesGRangesList))) {
      emptySamples <- rownames(allCellCounts)[
        !rownames(allCellCounts) %in% names(tilesGRangesList)
      ]
      emptyGRanges <- lapply(emptySamples, function(x) {
        NULL
      })
      names(emptyGRanges) <- emptySamples
      tilesGRangesList <- append(tilesGRangesList, emptyGRanges)
    }

    tilesGRangesList <- tilesGRangesList[sort(names(tilesGRangesList))]

    # Cannot make peak calls with < 5 cells (see make_prediction.R)
    # so NULL will occur for those samples. We need to fill in
    # dummy data so that we preserve the existence of the sample, while
    # also not including any information from it.
    emptyGroups <- which(unlist(lapply(tilesGRangesList, is.null)))
    if (length(emptyGroups) > 0) {
      if (verbose) {
        warning(
          paste(
            "The following celltype#sample groupings have too few cells (<5)",
            "and will be ignored: ", names(tilesGRangesList)[emptyGroups]
          )
        )
      }
    }

    for (i in emptyGroups) {
      # This is an empty region placeholder that represents an empty sample
      # And is functionally ignored downstream but required for
      # the RaggedExperiment structure
      tilesGRangesList[[i]] <- data.table(
        tileID = "chr1:1-499", seqnames = "chr1", start = 1,
        end = 499, strand = "*", TotalIntensity = 0, maxIntensity = 0,
        numCells = 0, Prediction = 0, PredictionStrength = 0, peak = FALSE
      )
    }

    # Package rangeList into a RaggedExperiment
    ragExp <- RaggedExperiment::RaggedExperiment(
      tilesGRangesList
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
    sampleDataFromCellColData(cellColData, sampleLabel = "Sample")
  )

  ## Finish processing the fragment number metadata for each sample-cell type combination.
  allFragmentCounts[is.na(allFragmentCounts)] = 0

  sampleCellTypeMetaData = SummarizedExperiment::SummarizedExperiment(
    append(
      list("CellCounts" = as.data.frame(allCellCounts)),
      list("FragmentCounts" = allFragmentCounts),
      additionalCellData
      ),
      colData = sampleData)

  # Add experimentList to MultiAssayExperiment
  names(experimentList) <- cellPopulations
  tileResults <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = experimentList,
    colData = sampleData,
    metadata = list(
      "Genome" = metadata(genome)$genome,
      "TxDb" = list(pkgname = TxDbName, metadata = S4Vectors::metadata(TxDb)),
      "OrgDb" = list(pkgname = OrgDbName, metadata = S4Vectors::metadata(OrgDb)),
      "SampleCellTypeData" = sampleCellTypeMetaData,
      "Directory" = outDir
    )
  )
  return(tileResults)
}