#' @title \code{getSampleTileMatrix}
#'
#' @description \code{getSampleTileMatrix} takes the output of peak calling with
#'   callOpenTiles and creates sample-tile matrices containing the signal
#'   intensity at each tile.
#'
#' @param tileResults a MultiAssayExperiment returned by callOpenTiles
#'   containing containing peak calling results.
#' @param cellPopulations vector of strings. Cell subsets in TileResults for
#'   which to generate sample-tile matrices. This list of group names must be
#'   identical to names that appear in the ArchRProject metadata.  If
#'   cellPopulations='ALL', then peak calling is done on all cell populations in
#'   the ArchR project metadata. Default is 'ALL'.
#' @param groupColumn Optional, the column containing sample group labels for
#'   determining consensus tiles within sample groups. Default is NULL, all
#'   samples will be used for determining consensus tiles.
#' @param threshold Threshold for consensus tiles, the minimum \% of samples
#'   (within a sample group, if groupColumn is set) that a peak must be called
#'   in to be retained. If set to 0, retain the union of all samples' peaks
#'   (this is equivalent to a threshold of 1/numSamples). It is recommended to
#'   tune this parameter to omit potentially spurious peaks.
#' @param tiles Optional, the number of cores to use with multiprocessing.
#'   Default is 1.
#' @param numCores Optional, the number of cores to use with multiprocessing.
#'   Default is 1.
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#'
#' @return SampleTileMatrices a MultiAssayExperiment containing a sample-tile
#'   intensity matrix for each cell population
#'
#' @examples
#' \donttest{
#' }
#'
#' @export
mergeTileResults <- function(tileResultsList, numCores = 1, verbose = TRUE) {

  #Test for duplicate sample names
  sampleTest <- unlist(lapply(tileResultsList, function(x) rownames(colData(x))))
  if (any(duplicated(sampleTest))) {
    stop("Sample names are duplicated in your list of tileResults. Please provide tileResults from different samples, not the same one's duplicated. ")
  }

  #Test whether all the tileResultsList indices are MultiAssayExperiments. 
  classTest <- lapply(tileResultsList, function(x) class(x)[1])
  if (any(unlist(classTest) != "MultiAssayExperiment")) {
    stop("At least one index of the tileResultsList is not MultiAssayExperiment")
  }

  #Test whether all the tileResults objects have the same Transcript, organism, and genome databases involved before merging. 
  TxDbTest <- unique(unlist(lapply(tileResultsList, function(x) x@metadata$TxDb$pkgname)))
  if (length(TxDbTest) > 1) {
    stop("Different TxDb were used to generates these tile results. Please re-run with the same TxDb")
  }
  
  OrgDbTest <- unique(unlist(lapply(tileResultsList, function(x) x@metadata$OrgDb$pkgname)))
  if (length(OrgDbTest) > 1) {
    stop("These tileResults are from different organisms and should not be merged. Why are you doing this to us?")
  }

  GenTest <- unique(unlist(lapply(tileResultsList, function(x) x@metadata$Genome)))
  if (length(GenTest) > 1) {
    stop("These tileResults are from different genome assemblies and cannot be merged.")
  }

  #Find all celltypes across all tile results objects. 
  allCellTypes <- as.data.frame(table(unlist(lapply(tileResultsList, names))))

  #Find the subset of cell types that are in common across all tileResults objects. Merge those. 
  subCellTypes <- as.character(unlist(dplyr::filter(allCellTypes, Freq == length(tileResultsList))$Var))

  if(length(subCellTypes) == 0){
    stop('No cell types in common across tileResults objects.')
  }
  if(verbose){
    message(paste0(subCellTypes, sep = ', '), 'are in common between all tileResults objects. These will be merged.')
  }

  #Iterate across each tileResults object, extracting the RaggedExperiments and turning them back into GRangesList. 
  #These GRangesList can then be concatenated easily, and then turned back into RaggedExperiments, and joined into one final tile results object. 
  cl = parallel::makeCluster(numCores)
  allRaggedResults <- lapply(subCellTypes, function(y){
            ragRes <- lapply(tileResultsList, function(x) x[[as.character(y)]])
            ragExp <- do.call('c', pbapply::pblapply(cl = cl, X = ragRes, combineRagged))
            RaggedExperiment::RaggedExperiment(ragExp)
  })
  names(allRaggedResults) <- subCellTypes

  parallel::stopCluster(cl)

  #Merged sample and other metadata information. 
  allSampleData <- do.call(eval(parse(text="dplyr::bind_rows")), lapply(tileResultsList, function(x){ as.data.frame(SummarizedExperiment::colData(x))}))
 
  #This is complicated, but it lets us handle table or data.frame outputs for this info. Merges cell counts and fragment count info across tile results. 
  allCellCounts <- do.call('cbind',pbapply::pblapply(cl = NULL, X = subCellTypes, function(y){
     
     tmpCell <- data.frame(unlist(lapply(tileResultsList, function(z){ z@metadata$CellCounts[,y] })))
     rownames(tmpCell) <- unlist(lapply(tileResultsList, function(z){ rownames(z@metadata$CellCounts)}))
     colnames(tmpCell) = y
     tmpCell
  }))
  allFragmentCounts <- do.call('cbind',pbapply::pblapply(cl = NULL, X = subCellTypes, function(y){

     tmpFrag <- data.frame(y = unlist(lapply(tileResultsList, function(x){ x@metadata$FragmentCounts[,y] })))
     rownames(tmpFrag) <- unlist(lapply(tileResultsList, function(x){ rownames(x@metadata$FragmentCounts)}))
     tmpFrag
  }))
  colnames(allCellCounts) <- colnames(allFragmentCounts) <-  subCellTypes
 
  #Construct final tile results object. 
  tileResults <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = allRaggedResults,
    colData = allSampleData,
    metadata = list(
      "CellCounts" = allCellCounts,
      "FragmentCounts" = allFragmentCounts,
      "Genome" = GenTest,
      "TxDb" = list(pkgname = tileResultsList[[1]]@metadata$TxDb$pkgname, metadata = TxDbTest),
      "OrgDb" = list(pkgname = tileResultsList[[1]]@metadata$OrgDb$pkgname, metadata = OrgDbTest),
      "Directory" = NULL
    )
  )
  return(tileResults)
 
}


combineRagged <- function(x){

    as(x, 'GRangesList')

}