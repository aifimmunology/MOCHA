#' @title \code{customSampleTileMatrix}
#'
#' @description \code{customSampleTileMatrix} takes the output of peak calling with
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
#' @param tiles A GRanges object or vector of strings in the format Chr9:100-200. These will be used as tiles for counting. 
#' @param numCores Optional, the number of cores to use with multiprocessing.
#'   Default is 1.
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#'
#' @return SampleTileMatrices a MultiAssayExperiment containing a sample-tile
#'   intensity matrix for each cell population
#'
#' @examples
#' \donttest{
#' # Starting from GRangesList
#' if (
#'   require(BSgenome.Hsapiens.UCSC.hg19) && 
#'   require(TxDb.Hsapiens.UCSC.hg38.refGene) && 
#'   require(org.Hs.eg.db)
#' ) {
#' tiles <- MOCHA::callOpenTiles(
#'   ATACFragments = MOCHA::exampleFragments,
#'   cellColData = MOCHA::exampleCellColData,
#'   blackList = MOCHA::exampleBlackList,
#'   genome = BSgenome.Hsapiens.UCSC.hg19,
#'   TxDb = TxDb.Hsapiens.UCSC.hg38.refGene,
#'   Org = org.Hs.eg.db,
#'   outDir = tempdir(),
#'   cellPopLabel = "Clusters",
#'   cellPopulations = c("C2", "C5"),
#'   numCores = 1
#' )
#' 
#' SampleTileMatrices <- MOCHA::customSampleTileMatrix(
#'   tiles,
#'   cellPopulations = c('C2', 'C5'),
#'   tiles = regionGRange
#' )
#' }
#' }
#'
#' @export
customSampleTileMatrix <- function(tileResults,
                                tiles = NULL,
                                cellPopulations = "ALL",
                                numCores = 1,
                                verbose = FALSE) {
  if (class(tileResults)[1] != "MultiAssayExperiment") {
    stop("tileResults is not a MultiAssayExperiment")
  }

  # Any column can be used to group samples
  # Note that these are case-sensitive
  sampleData <- MultiAssayExperiment::colData(tileResults)
 
  if (length(cellPopulations) == 1 & any(tolower(cellPopulations) == "all")) {
    subTileResults <- tileResults
    cellPopulations <- names(tileResults)
  } else {
    if (all(cellPopulations %in% names(tileResults))) {
      subTileResults <- MultiAssayExperiment::subsetByAssay(tileResults, cellPopulations)
    } else {
      stop(paste(
        "All of `cellPopulations` must present in tileResults.",
        "Check `names(tileResults)` for possible cell populations."
      ))
    }
  }

  #Set up regions to count. 

  if(class(tiles)[[1]] == 'GRanges'){
    tile_str <- MOCHA::StringsToGRanges(tiles)
    tile_gr = tiles
  }else{
    tile_gr <- MOCHA::GRangesToStrings(tiles)
    tile_str = tiles
  }

  if(!all(GenomicRanges::width(tile_gr) == 499 & 
      GenomicRanges::start(tile_gr) %% 500 == 0)){
    allTiles <- plyranges::reduce_ranges(tile_gr)
    RangeBins <- plyranges::stretch(plyranges::anchor_end(allTiles),
          extend = GenomicRanges::start(allTiles) %% 500)
    tile_gr = plyranges::slide_ranges(RangeBins,
                 width = 500, step = 500)
    tile_str <- MOCHA::StringsToGRanges(tile_gr)
  }


  if (verbose) {
    message(stringr::str_interp("Generating sample-tile matrix across all populations."))
  }
  # consensusTiles is used to  extract rows (tiles) from this matrix
  iterList <- lapply(seq_along(MultiAssayExperiment::experiments(subTileResults)), function(x){
    list(MultiAssayExperiment::experiments(subTileResults)[[x]], tile_str)
  })
  #browser()
  cl <- parallel::makeCluster(numCores)
  tilesByCellPop <- pbapply::pblapply(cl = cl, X = iterList, FUN = anyCalled)
  names(tilesByCellPop) <- names(subTileResults)

  tilesToCount <- unique(unlist(lapply(tilesByCellPop, names)))
  tilePresence <- do.call('cbind', lapply(tilesByCellPop, function(x){
    
    tilesToCount %in% names(x)
    
    }))
  rownames(tilePresence) = tilesToCount

  sampleTileIntensityMatList <- pbapply::pblapply(cl = cl, X = iterList, FUN = simplifiedSampleTile)
  names(sampleTileIntensityMatList) <- names(subTileResults)

  parallel::stopCluster(cl)

  # Order sampleData rows to match the same order as the columns
  maxMat <- which.max(lapply(sampleTileIntensityMatList, ncol))
  colOrder <- colnames(sampleTileIntensityMatList[[maxMat]])
  sampleData <- sampleData[match(colOrder, rownames(sampleData)), ]

  . <- NULL
  GenomicRanges::mcols(tile_gr) <- tilePresence

  results <- SummarizedExperiment::SummarizedExperiment(
    sampleTileIntensityMatList,
    rowRanges = allTilesGR,
    colData = sampleData,
    metadata = MultiAssayExperiment::metadata(tileRtileesults)
  )
  return(results)
}

## Helper function for addressing whether any peaks had been called
anyCalled <- function(xList){
  mat1 <- RaggedExperiment::compactAssay(xList[[1]], 'peak')
  trueTiles <- which(apply(mat1,1, any))
   return(trueTiles)
}