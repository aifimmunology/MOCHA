#' @title \code{getSampleTileMatrix}
#'
#' @description \code{getSampleTileMatrix} is a function that can transform
#'   a set of sample-specific tile calls into tile X sample matrix containing
#'   lambda1 measurements for each sample.
#'
#'
#' @param tileResults output of callOpenTiles
#' @return sampleTileIntensityMat a sample X tile matrix containing observed
#'   measurements for each sample at each tile.
#'
#' @details The technical details of the algorithm are found in XX.
#'
#' @references XX
#'
#' @export

getSampleTileMatrix <- function(tileResults,
                                cellPopulations = "ALL",
                                groupColumn = NULL,
                                threshold = 0.2,
                                join = "union",
                                NAtoZero = FALSE,
                                log2Intensity = TRUE,
				numCores = 1,
				verbose = FALSE
                               ) {

  if (!(join %in% c("union", "intersect"))){
    stop("`join` must be either 'union' or 'intersect'")
  }
  
  if(class(tileResults)[1] != 'MultiAssayExperiment'){
    stop('tileResults is not a MultiAssayExperiment')  
  }
  
  # Any column can be used to group samples
  # Note that these are case-sensitive
  sampleData <- MultiAssayExperiment::colData(tileResults)
  validGroups <- colnames(sampleData)
  if (!is.null(groupColumn)){
    if (!(groupColumn %in% validGroups)){
      stop("`groupColumn` not found in the column data of tileResults.")
    }
  }
  
  if (!(all(cellPopulations %in% names(tileResults))) & tolower(cellPopulations) != 'all'){
       stop(paste(
          "All of `cellPopulations` must present in tileResults.",
          "Check `names(tileResults)` for possible cell populations."
        ))
  }
  
  if(length(cellPopulations) == 1){
    if (cellPopulations == 'ALL'){
      subTileResults <- tileResults
      cellPopulations = names(tileResults)
    }
  }else{

    subTileResults <- tileResults[names(tileResults) %in% cellPopulations]

  }
  

  if(verbose){message(str_interp("Extracting consensus tile set within each population:  ${cellPopulations} "))}

  tilesByCellPop <- mclapply(experiments(subTileResults), function(x){
    
    # Get consensus tiles for this cell population for filtering
    scMACS:::singlePopulationConsensusTiles(
      x,
      threshold,
      groupColumn,
      join
    )

  }, mc.cores = numCores)

  allTiles <- sort(unique(do.call('c',tilesByCellPop)))

  if(verbose){message(str_interp("Generating Sample Tile Matrix Across All Populations: ${cellPopulations} "))}

  # consensusTiles is used to  extract rows (tiles) from this matrix
  sampleTileIntensityMatList <- mclapply(experiments(tileResults), function(x){

		scMACS:::singlePopulationSampleTileMatrix(
      				x,
      				allTiles,
      				NAtoZero,
      				log2Intensity
    		)
   }, mc.cores = numCores)

  tilePresence <- lapply(tilesByCellPop, function(x) (allTiles %in% x)) %>% do.call('cbind', .) %>% as.data.frame()

  allTilesGR <- scMACS::StringsToGRanges(allTiles)
  mcols(allTilesGR) <- tilePresence

  results <- SummarizedExperiment::SummarizedExperiment(
  		sampleTileIntensityMatList,
    		rowRanges = allTilesGR,
		colData = sampleData,
		metadata = append(list('Log2Intensity' = log2Intensity ), MultiAssayExperiment::metadata(tileResults))
  )

  return(results)
}



##### function for annotating the tiles for the Sample Tile Matrix.
## Returns a sample tile matrix with the transcript Db annotated and tile annotations. 


annotateSampleTileMatrix <- function(SampleTileMatrix,
			  TxDb = TxDb.Hsapiens.UCSC.hg38.refGene,
			   Org = org.Hs.eg.db, 
                          promoterRegion = c(2000, 100)){

	tileList <- rowRanges(SampleTileMatrix)

	metadata(SampleTileMatrix) = append(SummarizedExperiment::metadata(SampleTileMatrix), 
					list('Transcripts' = TxDb,
						'Org' = Org))

	tss_list <- suppressWarnings(ensembldb::transcriptsBy(TxDb, by = ('gene')))

        names(tss_list) <- suppressWarnings(MariaDB::mapIds(Org, names(tss_list), "SYMBOL", "ENTREZID"))

    	promoterList <- stack(tss_list) %>% GenomicRanges::trim(.) %>% 
		GenomicRanges::promoters(., upstream = promoterRegion[1], downstream = promoterRegion[2]) 

	

    ttiles <- plyranges::join_overlap_intersect(allT, GRanges)

	


}

#################################### getCellTypeMatrix pulls out the SampleTileMatrix of tiles called in one given cell type
## @SampleTileObj - output from getSampleTileMatrix, a SummarizedExperiment of pseudobulk intensitys across all tiles & cell types
## @CellType - The cell type you want to pull out.

getCellTypeMatrix <- function(SampleTileObj, CellType){
    
        tilesCalled <- SummarizedExperiment::mcols(rowRanges(SampleTileObj))[,CellType]
    
        sampleTileMatrix <- SummarizedExperiment::assays(SampleTileObj)[[CellType]][tilesCalled,]
    
        return(sampleTileMatrix)
    }

